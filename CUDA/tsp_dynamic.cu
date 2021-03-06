#include <stdio.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <cuda.h>

#define NUMBER_RUNS 1
#define LCG_A 1103515245
#define LCG_C 12345
#define LCG_M 2147483646
#define MAX_TRIES 1000
#define N_LIMIT 20
#define MAX_TEMP_STEPS 100
#define TEMP_START 0.5
#define COOLING 0.99
#define THREADS 1024
#define BOLTZMANN_COEFF 0.01

using namespace std;

struct city {
	double x;
	double y;
};

struct permutation {
	int cost;
	int* order;
	int nSucc;
};

struct GlobalConstants {
	int CITY_N;
	city* cities;
	unsigned int* randSeeds;
};

//global variables
struct city *cities;
int CITY_N;

//global variables on GPU
__constant__ GlobalConstants cuTspParam;

/* rounding function, but at .5 rounds to the lower int. Due to the TSPLIB
 * standard library.
 */
__device__ __host__ __inline__ int nint(float x)
{
	return (int) (x + 0.5);
}
/* Randomisation is done by a simple linear congruential generator.
 * We use A and C values as done by glibc.
 */
__device__ __inline__ unsigned int rand(unsigned int *x)
{
	*x = ((LCG_A * (*x)) + LCG_C) & 0x7fffffff;
	return *x;
}

__device__ __inline__ float randomFloat(unsigned int *x)
{
	return (float) (rand(x) / (float) LCG_M);
}

__device__ __inline__ double randomDouble(unsigned int *x)
{
	return (double) (rand(x) / (double) LCG_M);
}

__device__ __inline__ unsigned int randomInt(unsigned int *x, unsigned int max)
{
	return rand(x) % max;
}

__device__ __inline__ bool randomBool(unsigned int *x)
{
	if ((randomInt(x, 256) >> 7) & 0x00000001)
		return true;
	else
		return false;
}

__device__ __host__ __inline__ int euclideanDistance(struct city *a, struct city *b)
{
	float dx = b->x - a->x;
	float dy = b->y - a->y;
	return nint((sqrt(dx * dx + dy * dy)));
}

/* Calcuates the delta of the costs given by a new order using reverse
 */
__device__ __inline__ int reverseCost(struct city *cities, int *order, int *n)
{
	int cost;

	cost = -euclideanDistance(&cities[order[n[0]]], &cities[order[n[2]]]);
	cost -= euclideanDistance(&cities[order[n[1]]], &cities[order[n[3]]]);
	cost += euclideanDistance(&cities[order[n[0]]], &cities[order[n[3]]]);
	cost += euclideanDistance(&cities[order[n[1]]], &cities[order[n[2]]]);

	return cost;
}

/* The order of the city is changed by swapping the
 * order between n[0] and n[1].
 * The swapping is done beginning from the outer end
 * going into the middle
 */
__device__ __inline__ void reverse(int *order, int *n)
{
	int CITY_N = cuTspParam.CITY_N;
	int swaps = (1 + ((n[1] - n[0] + CITY_N) % CITY_N)) / 2;	// this many elements have to be swapped to have a complete reversal
	for (int j = 0; j < swaps; ++j) {
		int k = (n[0] + j) % CITY_N;
		int l = (n[1] - j + CITY_N) % CITY_N;
		int tmp = order[k];
		order[k] = order[l];
		order[l] = tmp;
	}
}

/* Calculates the delta of the costs of the city order if
 * the transportation of this segments (given by n) are actually
 * done.
 */
__device__ __inline__  int transportCost(struct city *cities, int *order, int *n)
{
	int cost;

	cost = -euclideanDistance(&cities[order[n[1]]], &cities[order[n[5]]]);
	cost -= euclideanDistance(&cities[order[n[0]]], &cities[order[n[4]]]);
	cost -= euclideanDistance(&cities[order[n[2]]], &cities[order[n[3]]]);
	cost += euclideanDistance(&cities[order[n[0]]], &cities[order[n[2]]]);
	cost += euclideanDistance(&cities[order[n[1]]], &cities[order[n[3]]]);
	cost += euclideanDistance(&cities[order[n[4]]], &cities[order[n[5]]]);

	return cost;
}

/* Transport the path segment (consisting of the start n[0] and end at n[1]
 * to the path given by n[2] and n[3], which are adjacent and the segment is
 * to be placed in between. n[4] is the city preceding n[0] and n[5] succeeds
 * n[1].
 * Transportation should only be done if the metroplis algorithm agrees.
 *
 */
__device__ void transport(int *order, int *n)
{
	int CITY_N = cuTspParam.CITY_N;
	int *newOrder = &order[CITY_N];
	int m1 = (n[1] - n[0] + CITY_N) % CITY_N;
	int m2 = (n[4] - n[3] + CITY_N) % CITY_N;
	int m3 = (n[2] - n[5] + CITY_N) % CITY_N;
	int i = 0;
	for (int j = 0; j <= m1; ++j) {
		newOrder[i++] = order[(j + n[0]) % CITY_N];
	}
	for (int j = 0; j <= m2; ++j) {
		newOrder[i++] = order[(j + n[3]) % CITY_N];
	}
	for (int j = 0; j <= m3; ++j) {
		newOrder[i++] = order[(j + n[5]) % CITY_N];
	}
	for (int j = 0; j < CITY_N; ++j) {
		order[j] = newOrder[j];
	}
}

/* Metroplis algorithm: Always take the downhill path and
 * sometime take the uphill path to avoid local minima
 */
__device__ __inline__ bool metropolis(const int cost, const double t, unsigned int *x)
{
	return cost < 0 || randomDouble(x) < exp((double) (BOLTZMANN_COEFF * -cost / t));
}

__host__ __inline__ void copy_permutation(struct permutation* dest, const struct permutation* src) {
	dest->cost = src->cost;
	dest->nSucc = src->nSucc;
	for (int i = 0; i < CITY_N; ++i) {
		dest->order[i] = src->order[i];
	}
}

/* Main kernel function */
__global__ void solve( struct permutation *permutations, const float t)
{
	struct city *cities = cuTspParam.cities;
	int CITY_N = cuTspParam.CITY_N;
	int notSeg;						// number of cities not on the segment
	int maxChangeTries = MAX_TRIES * CITY_N;
	int succLimit = N_LIMIT * CITY_N;
	int dCost;
	bool ans;
	int n[6];
	int id = blockDim.x * blockIdx.x + threadIdx.x;
	struct permutation *perm = &(permutations[id]);
	unsigned int *x = cuTspParam.randSeeds;

	perm->nSucc = 0;
	for (int j = 0; j < maxChangeTries; ++j) {
		do {
			n[0] = randomInt(x, CITY_N);
			n[1] = randomInt(x, CITY_N - 1);
			if (n[1] >= n[0])
				++n[1];
			notSeg = (n[0] - n[1] + CITY_N - 1) % CITY_N;
		} while (notSeg < 2);

		/* It is randomly choosen whether a transportation or a reversion is done */
		if (randomBool(x)) {
			n[2] = (n[1] + randomInt(x, abs(notSeg - 1)) + 1) % CITY_N;
			n[3] = (n[2] + 1) % CITY_N;
			n[4] = (n[0] + CITY_N- 1) % CITY_N;
			n[5] = (n[1] + 1) % CITY_N;

			dCost = transportCost(cities, perm->order, n);
			ans = metropolis(dCost, t, x);
			if (ans) {
				++perm->nSucc;
				perm->cost += dCost;
				transport(perm->order, n);
			}
		} else {
			n[2] = (n[0] + CITY_N - 1) % CITY_N;
			n[3] = (n[1] + 1) % CITY_N;

			dCost = reverseCost(cities, perm->order, n);
			ans = metropolis(dCost, t, x);
			if (ans) {
				++perm->nSucc;
				perm->cost += dCost;
				reverse(perm->order, n);
			}
		}

		/* Finish early if there are enough successful changes */
		if (perm->nSucc > succLimit)
			break;
	}
}

class Anneal {
private:
	/* Calculates the length of the initial path, which is already given.
	 * This is in O(n)
	 */
	void initialPath(struct permutation *perm, struct city *cities)
	{
		int i, i1, i2;

		perm->cost= 0;
		for (i = 0; i < CITY_N - 1; i++) {
			i1 = perm->order[i];
			i2 = perm->order[i+1];
			perm->cost += euclideanDistance(&cities[i1], &cities[i2]);
		}
		i1 = perm->order[CITY_N - 1];
		i2 = perm->order[0];
		perm->cost += euclideanDistance(&cities[i1], &cities[i2]);
		cout << "Initial path length: " << perm->cost << endl;
	}

	void printInformation(struct permutation *currPerm, bool showOrder = true)
	{
		cout << "Path Length = " << currPerm->cost << endl;
		cout << "Successful Moves: " << currPerm->nSucc << endl;
		if (showOrder) {
			cout << "Order: ";
			for (int j = 0; j < CITY_N; j++) {
				cout << currPerm->order[j] << " ";
			}
		}
		cout << endl;
	}

public:
	double runtime;
	int resultCost;

	Anneal() {}

	void order(struct city *cities, int *order)
	{
		double t = TEMP_START;
		long seed = (long) (time(NULL));
		cudaError_t cudaStat;
		struct permutation* dPermutation;
		struct permutation* hPermutation = (struct permutation *) malloc(THREADS * sizeof(struct permutation));
		for (int i = 0; i < CITY_N; i++) {
			hPermutation[i].order = new int [CITY_N];
		}
		struct city *dCities;
		unsigned int *LCGX = (unsigned int *) malloc(THREADS * sizeof(unsigned int));
		unsigned int *dLCGX;
		struct permutation *currPerm = (struct permutation *) malloc(sizeof(struct permutation));
		currPerm->order = new int [CITY_N];
		struct permutation *allMinPerm= (struct permutation *) malloc(sizeof(struct permutation));
		allMinPerm->order = new int [CITY_N];
		int oldCost = 2147483647;
		int repeatCost = 0;
		clock_t startAll, endAll;			// timer to measure the overall run time
		double runtimeAll;
		clock_t startCuda, endCuda;			//timer to measure the run time of cuda
		double cudaRuntime = 0.0f;

		startAll = clock();

		//initialize RNG
		srand(seed);

		//initialize seeds for RNG on GPU
		for (int i = 0; i < THREADS; ++i) {
			LCGX[i] = rand();
		}

		// Kernel invocation
		int threadsPerBlock = 256;
		int blocksPerGrid = (THREADS + threadsPerBlock - 1) / threadsPerBlock;

		cout << "Threads: " << THREADS << ", Blocks: " << blocksPerGrid << endl;

		for (int i = 0; i < CITY_N; i++) {
			currPerm->order[i] = order[i];
		}

		initialPath(currPerm, cities);
		copy_permutation(allMinPerm, currPerm);

		//allocate and copy #threads permutations on the device
		cudaStat = cudaMalloc(&dPermutation, THREADS * sizeof(struct permutation));
		if (cudaStat != cudaSuccess) {
			cout << "couldn't allocate memory on the device. Exit." << endl;
			return;
		}
		for (int i = 0; i < THREADS; i++) {
			int* order;
			cudaStat = cudaMalloc(&order, 2 * CITY_N * sizeof(int));
			if (cudaStat != cudaSuccess) {
				cout << "couldn't allocate memory on the device. Exit." << endl;
				return;
			}
			cudaStat = cudaMemcpy(order, currPerm->order, CITY_N * sizeof(int), cudaMemcpyHostToDevice);
			if (cudaStat != cudaSuccess) {
				cout << "couldn't allocate memory on the device. Exit." << endl;
				return;
			}
			cudaStat = cudaMemcpy(&dPermutation[i].order, &order, sizeof(int*), cudaMemcpyHostToDevice);
			if (cudaStat != cudaSuccess) {
				cout << "couldn't allocate memory on the device. Exit." << endl;
				return;
			}
		}
		cout<<"After the first call\n";

		cudaStat = cudaMalloc(&dCities, CITY_N * sizeof(struct city));
		if (cudaStat != cudaSuccess) {
			cout << "couldn't allocate memory on the device. Exit." << endl;
			return;
		}
		cudaStat = cudaMemcpy(dCities, cities, CITY_N * sizeof(struct city), cudaMemcpyHostToDevice);
		if (cudaStat != cudaSuccess) {
			cout << "couldn't copy memory to global memory. Exit." << endl;
			return;
		}

		cudaStat = cudaMalloc(&dLCGX, THREADS * sizeof(unsigned int));
		if (cudaStat != cudaSuccess) {
			cout << "couldn't allocate memory on the device. Exit." << endl;
			return;
		}
		cudaStat = cudaMemcpy(dLCGX, LCGX, THREADS * sizeof(unsigned int), cudaMemcpyHostToDevice);
		if (cudaStat != cudaSuccess) {
			cout << "couldn't copy memory to global memory. Exit." << endl;
			return;
		}

		cout<<"Just after the allocation\n";
		GlobalConstants params;
		params.cities = dCities;
		params.randSeeds = dLCGX;
		params.CITY_N = CITY_N;
		cudaMemcpyToSymbol(cuTspParam, &params, sizeof(GlobalConstants));
		cout<<"Just before the for loop\n";

		/* Try up to MAX_TEMP_STEPS temperature steps. It could stop before if no kernel
		 * showed any succesful change or if the solution did not change 5 times
		 */
		for (int i = 0; i < MAX_TEMP_STEPS; ++i) {
			cudaThreadSynchronize();
			startCuda = clock();

			cout<<"This is the "<<i<<" th loop\n";
			//Copies the initial permutation to each result permutation
			for (int i = 0; i < THREADS; i++) {
				//cudaStat = cudaMemcpy(dPermutation[i].order, currPerm->order, CITY_N * sizeof(int), cudaMemcpyHostToDevice);
				if (cudaStat != cudaSuccess) {
					cout << "couldn't copy memory to global memory. Exit." << endl;
					return;
				}
				cudaStat = cudaMemcpy(&dPermutation[i].cost, &currPerm->cost, sizeof(int), cudaMemcpyHostToDevice);
				if (cudaStat != cudaSuccess) {
					cout << "couldn't copy memory to global memory. Exit." << endl;
					return;
				}
			}

			cout<<"Just before the kernel call\n";
			//invoke cuda
			solve<<<blocksPerGrid, threadsPerBlock>>>(dPermutation, t);

			cudaStat = cudaThreadSynchronize();
			if (cudaStat != cudaSuccess) {
				cout << "something went wrong during device execution. Exit." << endl;
				return;
			}

			endCuda = clock();
			cudaRuntime += (endCuda - startCuda) * 1000 / CLOCKS_PER_SEC;

			for (int i = 0; i < THREADS; i++) {
				cudaStat = cudaMemcpy(hPermutation[i].order, dPermutation[i].order, CITY_N * sizeof(int), cudaMemcpyDeviceToHost);
				if (cudaStat != cudaSuccess) {
					cout << "couldn't copy memory from global memory. Exit." << endl;
					return;
				}
				cudaStat = cudaMemcpy(&hPermutation[i].cost, &dPermutation[i].cost, sizeof(int), cudaMemcpyDeviceToHost);
				if (cudaStat != cudaSuccess) {
					cout << "couldn't copy memory from global memory. Exit." << endl;
					return;
				}
				cudaStat = cudaMemcpy(&hPermutation[i].nSucc, &dPermutation[i].nSucc, sizeof(int), cudaMemcpyDeviceToHost);
				if (cudaStat != cudaSuccess) {
					cout << "couldn't copy memory from global memory. Exit." << endl;
					return;
				}
			}
			/* Loops through all resulting permutations and store the one with minimal length but
			 * at least one swap.
			 * If all threads didn't swap, exit the program.
			 * Takes O(n) time.
			 */
			int minCost = 2147483647;
			bool swap = false;
			for (int j = 0; j < THREADS; ++j) {
				if (minCost >= hPermutation[j].cost && hPermutation[j].nSucc != 0) {
					currPerm = &(hPermutation[j]);
					minCost = currPerm->cost;
					swap = true;
					if (minCost < allMinPerm->cost)
						copy_permutation(allMinPerm, currPerm);
						//memcpy(allMinPerm, currPerm, sizeof(struct permutation));
				}
			}

			if (!swap) {
				cout << "No swaps occured. Exit" << endl;
				break;
			}

			if (oldCost == minCost) {
				if (++repeatCost == 5) {
					cout << "Cost did not change 5 times in a row. Exit" << endl;
					break;
				}
			} else
				repeatCost = 0;

			cout << endl << "T = " <<  t << endl;
			//cout << "repeat: " << repeatCost << ", old: " << oldCost << ", new: " << minCost << endl;
			printInformation(currPerm, false);
			//for (int j = 0; j < THREADS; ++j)
			//	printInformation(&(hPermutation[j]), false);

			oldCost = minCost;
			t *= COOLING;
		}

		endAll = clock();
		runtimeAll = (endAll - startAll) / (1.0f * CLOCKS_PER_SEC) * 1000;

		cout << endl << "Final Result:" << endl;
		cout << "=============" << endl;
		printInformation(allMinPerm);

		runtime = runtimeAll;
		resultCost = allMinPerm->cost;

		printf("\nThe program needed an overall time of %.2lf ms.\n", runtimeAll);
		printf("%.2lf ms were spent at the CUDA part.\n", cudaRuntime);
		printf("So %.2lf ms were spent at the host.", runtimeAll - cudaRuntime);

		cudaFree(dPermutation);
		cudaFree(dCities);
		cudaFree(dLCGX);

		free(allMinPerm);
		free(LCGX);
		free(hPermutation);
	}
};

void readFile(char* FILENAME)
{
	FILE *fp;
	char line[80];
	int i = 0;

	fp = fopen(FILENAME, "rt");
	for (int i = 0; i < 6; i++) {
		fgets(line, 80, fp);
		if (i == 3) {
			char* pch = strtok(line, ":");
			pch = strtok(NULL, "\t\n");
			CITY_N = atoi(pch);
		}
	}
	cout<<"Number of cities is "<<CITY_N<<endl;
	cities = (struct city *) malloc (CITY_N * sizeof(struct city));
	while (fgets(line, 80, fp) != NULL && i < CITY_N) {
		sscanf(line, "%*d %lf %lf", &(cities[i].x), &(cities[i].y));
		++i;
	}
}

void printCities(struct city *cities)
{
	cout << "Cities: " << endl;
	for (int i = 0; i < CITY_N; ++i)
		cout << i << ". x: " << cities[i].x << " y: " << cities[i].y << endl;
}

int main(int argc, char* argv[])
{
	if (argc < 2) {
		printf("Usage: ./tsp inputFiles\n");
		exit(-1);
	}
	readFile(argv[1]);

	int *order = (int *) malloc(CITY_N * sizeof(int));
	float avgResult = 0.0f;
	double avgRuntime = 0.0f;

	for (int runs = 0; runs < NUMBER_RUNS; ++runs) {
		for (int i = 0; i < CITY_N; ++i)
			order[i] = i;

		//printCities(cities);

		Anneal *a = new Anneal();
		a->order(cities, order);
		avgResult += a->resultCost / (NUMBER_RUNS * 1.0f);
		avgRuntime += a->runtime / (NUMBER_RUNS * 1.0f);
	}

	cout << endl << endl;
	cout << "Average Costs: " << avgResult << endl;
	cout << "Average Runtime: " << avgRuntime << endl;

	return 0;
}
