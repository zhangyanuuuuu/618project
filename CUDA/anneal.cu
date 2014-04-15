#include <stdio.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <cuda.h>

#define MAX_TRIES 10
#define N_LIMIT 20
#define MAX_TEMP_STEPS 500
#define TEMP_START 0.5
#define COOLING 0.99
#define THREADS 1024
#define MAX_CITY 1024
#define BOLTZMANN_COEFF 0.01

#define CUDA_CALL(x) do {if((x) != cudaSuccess) {\
	printf("Error at %s:%d\n",__FILE__,__LINE__); \
	return EXIT_FAILURE;}} while(0)

using namespace std;

struct city {
	double x;
	double y;
};

struct permutation {
	int cost;
	int order[MAX_CITY];
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
__device__ __host__ int nint(float x)
{
	return (int) (x + 0.5);
}
/* Randomisation is done by a simple linear congruential generator.
 * We use A and C values as done by glibc.
 */

__device__ unsigned int randomInt(curandState *state, unsigned int max) {
	return curand(state) % max;
}

__device__ double randomDouble(curandState *state)
{
	return (double) curand_uniform(state);
}

__device__ bool randomBool(curandState *state)
{
	if ((randomInt(state, 256) >> 7) & 0x00000001)
		return true;
	else
		return false;
}

__global__ void initCurand(curandState *state, unsigned long seed) {
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	curand_init(seed, idx, 0, &state[idx]);
}

__device__ __host__ int euclideanDistance(struct city *a, struct city *b)
{
	float dx = b->x - a->x;
	float dy = b->y - a->y;
	return nint((sqrt(dx * dx + dy * dy)));
}

/* Calcuates the delta of the costs given by a new order using reverse
 */
__device__ int reverseCost(struct city *cities, int *order, int *n)
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
__device__ void reverse(int *order, int *n)
{
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
__device__ int transportCost(struct city *cities, int *order, int *n)
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
	int newOrder[MAX_CITY];
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
__device__ bool metropolis(const int cost, const double t, unsigned int *x)
{
	return cost < 0 || randomDouble(x) < exp((double) (BOLTZMANN_COEFF * -cost / t));
}

/* Main kernel function */
__global__ void solve(struct permutation *permutations, const float t)
{
	struct city* cities = cuTspParam.cities;
	int CITY_N = cuTspParam.CITY_N;
	int notSeg;						// number of cities not on the segment
	int maxChangeTries = MAX_TRIES * CITY_N;
	int succLimit = N_LIMIT * CITY_N;
	int dCost;
	bool ans;
	int n[6];
	int id = blockDim.x * blockIdx.x + threadIdx.x;
	struct permutation *perm = &(permutations[id]);
	unsigned int *x = &(lcg_x[id]);

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
		cudaError_t cudaStat;
		struct permutation *dPermutation;
		struct permutation *hPermutation = (struct permutation *) malloc(THREADS * sizeof(struct permutation));
		struct city *dCities;
		struct permutation *currPerm = (struct permutation *) malloc(sizeof(struct permutation));
		struct permutation *allMinPerm= (struct permutation *) malloc(sizeof(struct permutation));
		int oldCost = 2147483647;
		int repeatCost = 0;
		clock_t startAll, endAll;			// timer to measure the overall run time
		double runtimeAll;
		clock_t startCuda, endCuda;			//timer to measure the run time of cuda
		double cudaRuntime = 0.0f;
		curandState *devStates;

		startAll = clock();

		// Kernel invocation
		int threadsPerBlock = 256;
		int blocksPerGrid = (THREADS + threadsPerBlock - 1) / threadsPerBlock;

		cout << "Threads: " << THREADS << ", Blocks: " << blocksPerGrid << endl;

		memcpy(currPerm->order, order, CITY_N * sizeof(int));
		initialPath(currPerm, cities);
		memcpy(allMinPerm, currPerm, sizeof(struct permutation));

		CUDA_CALL(cudaMalloc(&dPermutation, THREADS * sizeof(struct permutation)));
		CUDA_CALL(cudaMalloc(&dCities, CITY_N * sizeof(struct city)));
		CUDA_CALL(cudaMemcpy(dCities, cities, CITY_N * sizeof(struct city), cudaMemcpyHostToDevice));

		// for generate random numbers directly on the device
		CUDA_CALL(cudaMalloc(void **)&devStates, THREADS * sizeof(curandState));
		initCurand<<<blocksPerGrid, threadsPerBlock>>>(devStates, 1234);
		
		/* Try up to MAX_TEMP_STEPS temperature steps. It could stop before if no kernel
		 * showed any succesful change or if the solution did not change 5 times
		 */
		for (int i = 0; i < MAX_TEMP_STEPS; ++i) {
			cudaDeviceSynchronize();
			startCuda = clock();

			//Copies the initial permutation to each result permutation
			for (int i = 0; i < THREADS; ++i) {
				memcpy(hPermutation[i].order, currPerm->order, CITY_N * sizeof(int));
				hPermutation[i].cost = currPerm->cost;
			}
			cudaStat = cudaMemcpy(dPermutation, hPermutation, THREADS * sizeof(struct permutation), cudaMemcpyHostToDevice);
			if (cudaStat != cudaSuccess) {
				cout << "couldn't copy memory to global memory. Exit." << endl;
				return;
			}

			//invoke cuda
			solve<<<blocksPerGrid, threadsPerBlock>>>(dPermutation, t);

			cudaStat = cudaThreadSynchronize();
			if (cudaStat != cudaSuccess) {
				cout << "something went wrong during device execution. Exit." << endl;
				return;
			}

			endCuda = clock();
			cudaRuntime += (endCuda - startCuda) * 1000 / CLOCKS_PER_SEC;

			cudaStat = cudaMemcpy(hPermutation, dPermutation, THREADS * sizeof(struct permutation), cudaMemcpyDeviceToHost);
			if (cudaStat != cudaSuccess) {
				cout << "couldn't copy memory from global memory. Exit." << endl;
				return;
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
						memcpy(allMinPerm, currPerm, sizeof(struct permutation));
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

		free(allMinPerm);
		free(hPermutation);
	}
};
