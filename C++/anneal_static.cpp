#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <cstring>

#define LCG_A 1103515245
#define LCG_C 12345
#define LCG_M 2147483646
#define MAX_TRIES 10240
#define N_LIMIT 20
#define MAX_TEMP_STEPS 500
#define TEMP_START 20
#define COOLING 0.95
#define CITY_N 280
#define BOLTZMANN_COEFF 0.1

using namespace std;

struct city {
	double x;
	double y;
};

struct permutation {
	int cost;
	int order[CITY_N];
	int nSucc;
};

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
	int newOrder[CITY_N];
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
bool metropolis(const int cost, const double t, unsigned int *x)
{
	return cost < 0 || myrandomDouble(x) < exp((double) (BOLTZMANN_COEFF * -cost / t));
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
		struct city *dCities;
		struct permutation *currPerm = (struct permutation *) malloc(sizeof(struct permutation));
		struct permutation *allMinPerm = (struct permutation *) malloc(sizeof(struct permutation));
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

		startAll = clock();

		//initialize RNG
		unsigned int xval = (unsigned int) seed;
		unsigned int *x = &xval;

		for (int j = 0; j < CITY_N; ++j)
			currPerm->order[j] = order[j];

		initialPath(currPerm, cities);
		memcpy(allMinPerm, currPerm, sizeof(struct permutation));

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

			if (currPerm->cost < allMinPerm->cost)
				memcpy(allMinPerm, currPerm, sizeof(struct permutation));

			if (currPerm->nSucc == 0) {
				cout << "No swaps occured. Exit" << endl;
				break;
			}

			if (oldCost == currPerm->cost) {
				if (++repeatCost == 5) {
					cout << "Cost did not change 5 times in a row. Exit" << endl;
					break;
				}
			} else
				repeatCost = 0;

			cout << endl << "T = " <<  t << endl;
			//cout << "repeat: " << repeatCost << ", old: " << oldCost << ", new: " << minCost << endl;
			printInformation(currPerm, false);

			oldCost = currPerm->cost;
			t *= COOLING;
		}

		resultCost = oldCost;
		endAll = clock();
		runtimeAll = (endAll - startAll) / (1.0f * CLOCKS_PER_SEC) * 1000;

                cout << endl << "Final Result:" << endl;
                cout << "=============" << endl;
                printInformation(allMinPerm);

                runtime = runtimeAll;
                resultCost = allMinPerm->cost;

		printf("The program needed an overall time of %.2lf ms.\n", runtimeAll);

		free(currPerm);
		free(allMinPerm);
	}
};
