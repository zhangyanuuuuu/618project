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
#define N_LIMIT 200
#define MAX_TEMP_STEPS 500
#define TEMP_START 20
#define COOLING 0.95
#define BOLTZMANN_COEFF 1

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

/* rounding function, but at .5 rounds to the lower int. Due to the TSPLIB
 * standard library.
 */
int nint(double x)
{
	return (int) (x + 0.5);
}
/* Randomisation is done by a simple linear congruential generator.
 * We use A and C values as done by glibc.
 */

double myrandomDouble()
{
	return (double) rand() / RAND_MAX;
}

unsigned int myrandomInt(unsigned int max)
{
	return rand() % max;
}

int euclideanDistance(struct city *a, struct city *b)
{
	double dx = b->x - a->x;
	double dy = b->y - a->y;
	return nint((sqrt(dx * dx + dy * dy)));
}


/* Metroplis algorithm: Always take the downhill path and
 * sometime take the uphill path to avoid local minima
 */
bool metropolis(const int cost, const double t)
{
	return cost < 0 || myrandomDouble() < exp((double) (BOLTZMANN_COEFF * -cost / t));
}

double prevLinkCost(int i, int* order) {
   
    double delX, delY;
   int otherEnd; 
   // this fussy business handles the case when city i
   // is the very FIRST city in our order, and the "previous"
   // city is the LAST city
   if( i==0) 
     otherEnd = CITY_N-1;
   else otherEnd = i-1;

   delX = cities[order[i]].x - cities[order[otherEnd]].x;
   delX = delX*delX;
   delY = cities[order[i]].y - cities[order[otherEnd]].y;
   delY = delY*delY;

   return sqrt( delX + delY );
} 

// compute length of the TSP link the "starts" at city i
double nextLinkCost(int i, int* order) {
   
   double delX, delY;
   int otherEnd;
   // this fussy business handles the case when city i
   // is the very LAST city in our order, and the "previous"
   // city is the FIRST city 
   if( i==CITY_N-1) 
     otherEnd = 0;
   else otherEnd = i+1;

   delX = cities[order[i]].x - cities[order[otherEnd]].x;
   delX = delX*delX;
   delY = cities[order[i]].y - cities[order[otherEnd]].y;
   delY = delY*delY;

   return sqrt( delX + delY );
}  

int evalSwap(int i,  int j, int* order)
{
double afterCost, beforeCost;

    // case 1:  i, j not adjacent
    if ( abs(i-j) != 1
         || ( (i==0) && (j!=(CITY_N-1)) )
         || ( (j==0) && (i!=(CITY_N-1)) ) ) {

        // eval the cost of the 4 critical links
        beforeCost = prevLinkCost(i, order) + nextLinkCost(i, order) 
                   + prevLinkCost(j, order) + nextLinkCost(j, order);

        // now, swap the cities :
        //    city[i] <==> city[j] in order order
        int tmp= order[i];
        order[i] = order[j];
        order[j] = tmp;

        // eval the changed order length.  Look only at 4 critical links
        afterCost = prevLinkCost(i, order) + nextLinkCost(i, order) 
                   + prevLinkCost(j, order) + nextLinkCost(j, order);

        // now, UN-swap the cities  back
        //    city[j] <==> city[i]
        tmp= order[i];
        order[i] = order[j];
        order[j] = tmp;

        // return the delta in the order cost
        return  nint(afterCost - beforeCost);
    }

    if ( abs(i-j) == 1
         || ( (i==0) && (j==(CITY_N-1)) )
         || ( (j==0) && (i==(CITY_N-1)) ) ) {

        // we need to know which one is first in the order
        int minCity = min(i,j);
        int maxCity = max(i,j);

        // eval the cost of the 3 critical links
        if( i!=0  && j!= 0) {
            // its the sub-case illustrated in comment above
            beforeCost = prevLinkCost(minCity, order) + nextLinkCost(minCity, order) 
                       + nextLinkCost(maxCity, order);
        }
        else  {
            // cities on the "ends" like:    i 1 2 3 ... n-3 n-2 j
            // we will swap to this:         j 1 2 3 ... n-3 n-2 i
            // eval the proper link cost, note its different
            beforeCost = prevLinkCost(minCity, order) + nextLinkCost(minCity, order) 
                       + prevLinkCost(maxCity, order);
        }

        // now, swap the cities :
        //    city[i] <==> city[j] in order order
        int tmp= order[i];
        order[i] = order[j];
        order[j] = tmp;

        // eval ONLY the changed order lengths.  Look only at 3 critical links
        minCity = min(i,j);
        maxCity = max(i,j);
        // eval the cost of the 3 critical links
        if( i!=0  && j!= 0) {
            // its the sub-case illustrated in first comment in this routine
            afterCost =  prevLinkCost(minCity, order) + nextLinkCost(minCity, order) 
                       + nextLinkCost(maxCity, order);
        }
        else  {
            // cities on the "ends" like:    i 1 2 3 ... n-3 n-2 j
            // we will swap to this:         j 1 2 3 ... n-3 n-2 i
            // eval the proper link cost, note its different
            afterCost =  prevLinkCost(minCity, order) + nextLinkCost(minCity, order) 
                       + prevLinkCost(maxCity, order);
        }

        // now, UN-swap the cities  back
        //    city[j] <==> city[i]
        tmp= order[i];
        order[i] = order[j];
        order[j] = tmp;

        // return the delta in the order cost
        return  nint(afterCost - beforeCost);
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
		srand(time(NULL));
		struct city *dCities;
		struct permutation *currPerm = (struct permutation *) malloc(sizeof(struct permutation));
		currPerm->order = new int [CITY_N];
		struct permutation *allMinPerm = (struct permutation *) malloc(sizeof(struct permutation));
		allMinPerm->order = new int [CITY_N];
		int oldCost = 2147483647;
		int repeatCost = 0;
		clock_t startAll, endAll;			// timer to measure the overall run time
		double runtimeAll;
		int notSeg;						// number of cities not on the segment
		int maxChangeTries = MAX_TRIES * CITY_N;
		int succLimit = N_LIMIT * CITY_N;
		int dCost;
		bool ans;
		int n[2];

		startAll = clock();

		//initialize RNG

		for (int j = 0; j < CITY_N; ++j)
			currPerm->order[j] = order[j];

		initialPath(currPerm, cities);
		allMinPerm->nSucc = currPerm->nSucc;
	 	allMinPerm->cost = currPerm->cost;	
		for (int i = 0; i < CITY_N; i++) {
			allMinPerm->order[i] = currPerm->order[i];
		}

		/* Try up to MAX_TEMP_STEPS temperature steps. It could stop before if no kernel
		 * showed any succesful change or if the solution did not change 5 times
		 */
		for (int i = 0; i < MAX_TEMP_STEPS; ++i) {
			currPerm->nSucc = 0;
			for (int j = 0; j < maxChangeTries; ++j) {
				n[0] = myrandomInt(CITY_N);
				n[1] = myrandomInt(CITY_N);
				if (n[0] == n[1]) {
					if (n[0] - 1 > 0)
						n[1] = n[0] - 1;
					else 
						n[1] = n[0] + 1;
				}
				dCost = evalSwap(n[0], n[1], currPerm->order);
				ans = metropolis(dCost, t);
				if (ans) {
					++currPerm->nSucc;
					currPerm->cost += dCost;
					cout<<"The change of cost of node "<<n[0]<<" "<<n[1]<<" is "<<dCost<<endl;
					cout<<"The cost now is "<<currPerm->cost<<endl;
					//perform the swap

					int tmp = currPerm->order[n[0]];
					currPerm->order[n[0]] = currPerm->order[n[1]];
					currPerm->order[n[1]] = tmp;
				}
				/* Finish early if there are enough successful changes */
				if (currPerm->nSucc > succLimit)
					break;
			}

			if (currPerm->cost < allMinPerm->cost) {
				allMinPerm->nSucc = currPerm->nSucc;
	 			allMinPerm->cost = currPerm->cost;	
				for (int i = 0; i < CITY_N; i++) {
					allMinPerm->order[i] = currPerm->order[i];
				}
			}

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
