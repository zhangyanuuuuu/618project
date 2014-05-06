#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <mpi.h>
#include <assert.h>
#include <algorithm>

#define ROOT 0
#define LCG_A 1103515245
#define LCG_C 12345
#define LCG_M 2147483646
#define MAX_TRIES 10240
#define N_LIMIT 20
#define MAX_TEMP_STEPS 500
#define TEMP_START 3
#define COOLING 0.95
#define BOLTZMANN_COEFF 0.1
const float initialMutationRate = 0.3;
const float mutationRate = 0.4;
const int generationNum = 700;

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

//global variable
struct city *cities;
int CITY_N;
int size;
int selectionNum;

/* rounding function, but at .5 rounds to the lower int. Due to the TSPLIB
 * standard library.
 */
int nint(float x)
{
	return (int) (x + 0.5);
}
/* Randomisation is done by a simple linear congruential generator.
 * We use A and C values as done by glibc.
 */
unsigned int myrand(unsigned int *x)
{
	*x = ((LCG_A * (*x)) + LCG_C) & 0x7fffffff;
	return *x;
}

float myrandomFloat(unsigned int *x)
{
	return (float) (myrand(x) / (float) LCG_M);
}

double myrandomDouble(unsigned int *x)
{
	return (double) (myrand(x) / (double) LCG_M);
}

unsigned int myrandomInt(unsigned int *x, unsigned int max)
{
	return myrand(x) % max;
}

bool myrandomBool(unsigned int *x)
{
	if ((myrandomInt(x, 256) >> 7) & 0x00000001)
		return true;
	else
		return false;
}

int euclideanDistance(struct city *a, struct city *b)
{
	float dx = b->x - a->x;
	float dy = b->y - a->y;
	return nint((sqrt(dx * dx + dy * dy)));
}

/* Calcuates the delta of the costs given by a new order using reverse
 */
int reverseCost(struct city *cities, int *order, int *n)
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
void reverse(int *order, int *n)
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
int transportCost(struct city *cities, int *order, int *n)
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
void transport(int *order, int *n)
{
	int* newOrder = new int [CITY_N];
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
	delete [] newOrder;
}

/* Metroplis algorithm: Always take the downhill path and
 * sometime take the uphill path to avoid local minima
 */
bool metropolis(const int cost, const double t, unsigned int *x)
{
	return cost < 0 || myrandomDouble(x) < exp((double) (BOLTZMANN_COEFF * -cost / t));
}


//functions for GA

void mutate(int* tour, float rate) {
	int pos2, temp;
	for (int pos1 = 0; pos1 < CITY_N; pos1++) {
		if ((static_cast<float> (rand()) / static_cast<float> (RAND_MAX)) < rate) {
			pos2 = rand() % CITY_N;
			temp = tour[pos1];
			tour[pos1] = tour[pos2];
			tour[pos2] = temp;
		}
	}
}

void optMutate(int* tour) {
	int n[4];
	int notSeg;
	int dCost;
	for (int i = 0; i < (int) ((float)CITY_N * mutationRate); i++) {
		do {
			n[0] = rand() % CITY_N;
			n[1] = rand() % (CITY_N - 1);
			if (n[1] >= n[0]) 
				++n[1];
			notSeg = (n[0] - n[1] + CITY_N - 1) % CITY_N;
		} while (notSeg < 2);
		n[2] = (n[0] + CITY_N - 1) % CITY_N;
		n[3] = (n[1] + 1) % CITY_N;

		dCost = reverseCost(cities, tour, n);
		if (dCost < 0) {
			reverse(tour, n);
		}
	}
}

int* crossOver(int* parent1, int* parent2) {
	bool* child_contains = new bool [CITY_N];
	for (int i = 0; i < CITY_N; i++) {
		child_contains[i] = false;
	}
	int* child = new int [CITY_N];
	int startPos = rand() % CITY_N;
	int endPos = rand() % CITY_N;

	if (startPos > endPos) {
		int temp = startPos;
		startPos = endPos;
		endPos = temp;
	}
#ifdef DEBUG
	cout<<"start position is "<<startPos<<endl;
	cout<<"end position is "<<endPos<<endl;
#endif
	
	// Loop and add the sub tour from parent1 to child
	for (int i = startPos; i <= endPos; i++) {
		child[i] = parent1[i];
		child_contains[child[i]] = true;
	}

	int k = endPos + 1;
	//first loop through the end of the tour
	for (int i = endPos + 1; i < CITY_N; i++) {
		if (!child_contains[parent2[i]]) {
			child[k++ % CITY_N] = parent2[i];
			child_contains[parent2[i]] = true;
		}
	}

	//then go through the start of the tour
	for (int i = 0; i <= endPos; i++) {
		if (!child_contains[parent2[i]]) {
			child[k++ % CITY_N] = parent2[i];
			child_contains[parent2[i]] = true;
		}
	}

#ifdef DEBUG
	assert(k >= startPos);
#endif
	delete [] child_contains;
	return child;
}

bool checkTour(int* tour, int CITY_N) {
	bool* test = new bool [CITY_N];
	for (int i = 0; i < CITY_N; i++) {
		test[i] = false;
	}
	for (int i = 0; i < CITY_N; i++) {
		test[tour[i]] = true;
	}
	for (int i = 0; i < CITY_N; i++) {
		if (!test[i]) {
			return false;
		}
	}
	return true;
}

void printTour(int* tour) {
	for (int i = 0; i < CITY_N; i++) {
		cout<<tour[i]<<' ';
	}
	cout<<endl;
}

int calculate_cost(int* order, struct city *cities) {
	int cost = 0;
	int i1, i2;
	for (int i = 0; i < CITY_N - 1; i++) {
		i1 = order[i];
		i2 = order[i+1];
		cost += euclideanDistance(&cities[i1], &cities[i2]);
	}
	i1 = order[CITY_N - 1];
	i2 = order[0];
	cost += euclideanDistance(&cities[i1], &cities[i2]);
#ifdef DEBUG
	cout<<"New path cost: "<<cost<<endl;
#endif
	return cost;
}

int tournamentSelection(int* cost) {
	int randomIdx = rand() % size;
	int minIdx = randomIdx;
	int minCost = cost[minIdx];
	for (int i = 0; i <selectionNum; i++) {
		randomIdx = rand() % size;
		if (cost[randomIdx] < minCost) {
			minCost = cost[randomIdx];
			minIdx = randomIdx;
		}
	}
	return minIdx;
}

int** evolveGeneration(int** generation, int*& cost) {
	int **new_generation = new int* [size];
	int *new_cost = new int [size];
	int minIdx = 0;
	int minCost = cost[0];
	for (int i = 1; i < size; i++) {
		if (cost[i] < minCost) {
			minIdx = i;
			minCost = cost[i];
		}
	}
	new_generation[0] = generation[minIdx];
	new_cost[0] = minCost;
	for (int i = 1; i < size; i++) {
		int idx1 = tournamentSelection(cost);
		int idx2 = tournamentSelection(cost);
		while (idx2 == idx1) {
			idx2 = tournamentSelection(cost);
		}
		int *parent1 = generation[idx1];
		int *parent2 = generation[idx2];
#ifdef DEBUG
		cout<<"parent 1 is "<<idx1<<"; parent 2 is "<<idx2<<endl;
		printTour(parent1);
		printTour(parent2);
#endif
		new_generation[i] = crossOver(parent1, parent2);
		optMutate(new_generation[i]);
		new_cost[i] = calculate_cost(new_generation[i], cities);
	}

	for (int i = 0; i < size; i++) {
		if (i != minIdx) {
			delete [] generation[i];
		}
	}
	delete [] generation;
	delete [] cost;
	cost = new_cost;
	return new_generation;
}

void printGeneration(int** generation) {
	for (int i = 0; i < size; i++) {
		printTour(generation[i]);
	}
}

void printAvgCost(int* cost) {
	int totalCost = 0;
	for (int i = 0; i < size; i++) {
		totalCost += cost[i];
	}
	cout<<"average cost for this generation is "<< totalCost/size <<endl;
}

void printCost(int* cost) {
	int totalCost = 0;
	for (int i = 0; i < size; i++) {
		totalCost += cost[i];
		cout<<"cost for the "<<i<<"th element is "<<cost[i]<<endl;
	}
	cout<<"average cost for this generation is "<< totalCost/size <<endl;
}
	
void GA(int* global_order, int*& cost) {
	int** generation = new int* [size];
	for (int i = 0; i < size; i++) {
		generation[i] = new int [CITY_N];
		memcpy(generation[i], global_order+i*CITY_N, CITY_N*sizeof(int));
	}
	for (int i = 0; i < generationNum; i++) {
		generation = evolveGeneration(generation, cost);
	}
	for (int i = 0; i < size; i++) {
		memcpy(global_order+i*CITY_N, generation[i],CITY_N*sizeof(int));
		delete [] generation[i];
	}
	delete [] generation;
}		

void fullGA(int* global_order, int*& cost) {
	int** generation = new int* [size];
	for (int i = 0; i < size; i++) {
		generation[i] = new int [CITY_N];
		memcpy(generation[i], global_order+i*CITY_N, CITY_N*sizeof(int));
	}
	int bestCost = cost[0];
	int repeatCost = 0;
	while(1) {
		generation = evolveGeneration(generation, cost);
		//printCost(cost);
		if (cost[0] < bestCost) {
			bestCost = cost[0];
			repeatCost = 0;
		} else if (cost[0] == bestCost) {
			if (++repeatCost == 100) {
				cout<<"Cost did not change 100 times in a row. Exit"<<endl;
				break;
			}
		} else {
			repeatCost = 0;
		}
	}
	for (int i = 0; i < size; i++) {
		memcpy(global_order+i*CITY_N, generation[i],CITY_N*sizeof(int));
		delete [] generation[i];
	}
	delete [] generation;
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
	int rank;

	Anneal(int rank) {
		this->rank = rank;
	}

	void order(struct city *cities, int *order)
	{
		double t = TEMP_START;
		long seed = (long) (time(NULL));
		struct permutation *currPerm = new struct permutation;
		currPerm->order = new int [CITY_N];
		int oldCost = 2147483647;
		int repeatCost = 0;
		clock_t startAll, endAll;			// timer to measure the overall run time
		double runtimeAll;
		int notSeg;						// number of cities not on the segment
		int maxChangeTries = MAX_TRIES * CITY_N;
		int succLimit = N_LIMIT * CITY_N;
		int dCost;
		bool ans;
		int n[6];
		int *global_cost;
		int *global_order;
		MPI_Status status;
		int minIdx;
		char out = 0;

		if (rank == ROOT) {
			global_cost = new int [size];
			global_order = new int [size * CITY_N];
		}

		startAll = clock();

		//initialize RNG
		unsigned int xval = (unsigned int) seed * rank;
		unsigned int *x = &xval;

		for (int j = 0; j < CITY_N; ++j)
			currPerm->order[j] = order[j];

		initialPath(currPerm, cities);
		
		struct permutation *allMinPerm;
		
		if (rank == ROOT) { //only root keeps the best permutation so far
			allMinPerm = new struct permutation;
			allMinPerm->order = new int [CITY_N];
		}

		//full GA at the beginning of the loop
		assert(MPI_Gather(&currPerm->cost, 1, MPI_INT, global_cost, 1, MPI_INT, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS);
		assert(MPI_Gather(currPerm->order, CITY_N, MPI_INT, global_order, CITY_N, MPI_INT, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS);
			
		if (rank == ROOT) {
			fullGA(global_order, global_cost);
			memcpy(currPerm->order, global_order, CITY_N*sizeof(int));
			currPerm->cost = global_cost[0];
			memcpy(allMinPerm->order, currPerm->order, CITY_N*sizeof(int));
			allMinPerm->cost = currPerm->cost;
		}

		assert(MPI_Bcast(currPerm->order, CITY_N, MPI_INT, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS);
		assert(MPI_Bcast(&currPerm->cost, 1, MPI_INT, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS);
				
		/* Try up to MAX_TEMP_STEPS temperature steps. It could stop before if no kernel
		 * showed any succesful change or if the solution did not change 5 times
		 */
		for (int i = 0; i < MAX_TEMP_STEPS; ++i) {
			currPerm->nSucc = 0;
			for (int j = 0; j < maxChangeTries; ++j) {
				do {
					n[0] = myrandomInt(x, CITY_N);
					n[1] = myrandomInt(x, CITY_N - 1);
					if (n[1] >= n[0]) 
						++n[1];
					notSeg = (n[0] - n[1] + CITY_N - 1) % CITY_N;
				} while (notSeg < 2);

				/* It is myrandomly choosen whether a transportation or a reversion is done */
				if (myrandomBool(x)) {
					n[2] = (n[1] + myrandomInt(x, abs(notSeg - 1)) + 1) % CITY_N;
					n[3] = (n[2] + 1) % CITY_N;
					n[4] = (n[0] + CITY_N- 1) % CITY_N;
					n[5] = (n[1] + 1) % CITY_N;

					dCost = transportCost(cities, currPerm->order, n);
					ans = metropolis(dCost, t, x);
					if (ans) {
						++currPerm->nSucc;
						currPerm->cost += dCost;
						transport(currPerm->order, n);
					}
				} else {
					n[2] = (n[0] + CITY_N - 1) % CITY_N;
					n[3] = (n[1] + 1) % CITY_N;

					dCost = reverseCost(cities, currPerm->order, n);
					ans = metropolis(dCost, t, x);
					if (ans) {
						++currPerm->nSucc;
						currPerm->cost += dCost;
						reverse(currPerm->order, n);
					}
				}

				/* Finish early if there are enough successful changes */
				if (currPerm->nSucc > succLimit)
					break;
			}
			
			//Gather all the cost and order from all the processes
			assert(MPI_Gather(&currPerm->cost, 1, MPI_INT, global_cost, 1, MPI_INT, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS);
			assert(MPI_Gather(currPerm->order, CITY_N, MPI_INT, global_order, CITY_N, MPI_INT, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS);
			
			if (rank == ROOT) {
				//find the best in all of the tours
        minIdx = 0;
        int minCost = currPerm->cost;
        for (int i = 1; i < size; i++) {
        	if (minCost > global_cost[i]) {
          	minCost = global_cost[i];
            minIdx = i;
          }
        }
				if (minCost < allMinPerm->cost) {
	 				allMinPerm->cost = minCost;	
					memcpy(allMinPerm->order, global_order+CITY_N*minIdx, CITY_N*sizeof(int));
				}
				//apply GA to get next xth generation
				fullGA(global_order,global_cost);
				memcpy(currPerm->order, global_order, CITY_N*sizeof(int));
				currPerm->cost = global_cost[0];
			}

			//scatter to all the processes
			//assert(MPI_Scatter(global_order, CITY_N, MPI_INT, currPerm->order, CITY_N, MPI_INT, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS);
			//assert(MPI_Scatter(global_cost, 1, MPI_INT, &currPerm->cost, 1, MPI_INT, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS);
			assert(MPI_Bcast(currPerm->order, CITY_N, MPI_INT, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS);
			assert(MPI_Bcast(&currPerm->cost, 1, MPI_INT, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS);
			//assert(MPI_Bcast(&currPerm->nSucc, 1, MPI_INT, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS);
						
			if (rank == ROOT) {
				if (currPerm->nSucc == 0) {
					cout << "No swaps occured. Exit" << endl;
					out = 1;
				}

				if (abs(oldCost -  currPerm->cost) <= 3) {
					if (++repeatCost == 5) {
						cout << "Cost did not change 5 times in a row. Exit" << endl;
						out = 1;
					}
				} else
					repeatCost = 0;

				cout << endl << "T = " <<  t << endl;
				printInformation(currPerm, false);
			}
			
			assert(MPI_Bcast(&out, 1, MPI_CHAR, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS);
			if (out == 1) {
				break;
			}

			oldCost = currPerm->cost;
			t *= COOLING;
		}

		resultCost = oldCost;
		endAll = clock();
		if (rank == ROOT) {
			delete [] global_cost;
			runtimeAll = (endAll - startAll) / (1.0f * CLOCKS_PER_SEC) * 1000;

                cout << endl << "Final Result:" << endl;
                cout << "=============" << endl;
                printInformation(allMinPerm);

                runtime = runtimeAll;
                resultCost = allMinPerm->cost;

			printf("The program needed an overall time of %.2lf ms.\n", runtimeAll);
		}
	
		delete [] currPerm->order;
		delete currPerm;
		
		if (rank == ROOT) {
			delete [] allMinPerm->order;
			delete allMinPerm;
		}
	}
};
