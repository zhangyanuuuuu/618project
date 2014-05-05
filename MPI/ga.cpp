#include <iostream>
#include <assert.h>
#include <time.h>
#include <algorithm>
#include <math.h>
using namespace std;

struct city {
	double x;
	double y;
};

struct city *cities;
int CITY_N;
const int size = 10;
const float initialMutationRate = 0.3;
const float mutationRate = 0.1;
const float selectionFactor = 0.4;
const int generationNum = 500;

int nint(float x)
{
	return (int) (x + 0.5);
}

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

int euclideanDistance(struct city *a, struct city *b)
{
	float dx = b->x - a->x;
	float dy = b->y - a->y;
	return nint((sqrt(dx * dx + dy * dy)));
}

int reverseCost(struct city *cities, int *order, int *n)
{
	int cost;

	cost = -euclideanDistance(&cities[order[n[0]]], &cities[order[n[2]]]);
	cost -= euclideanDistance(&cities[order[n[1]]], &cities[order[n[3]]]);
	cost += euclideanDistance(&cities[order[n[0]]], &cities[order[n[3]]]);
	cost += euclideanDistance(&cities[order[n[1]]], &cities[order[n[2]]]);

	return cost;
}

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

	cout<<"start position is "<<startPos<<endl;
	cout<<"end position is "<<endPos<<endl;
	
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
	cout<<"New path cost: "<<cost<<endl;
	return cost;
}

int tournamentSelection(int* cost) {

	int randomIdx = rand() % size;
	int minIdx = randomIdx;
	int minCost = cost[minIdx];
	for (int i = 0; i <(int) ((float)size * selectionFactor); i++) {
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
		if (idx2 == idx1) {
			idx2 = tournamentSelection(cost);
		}
		cout<<"parent 1 is "<<idx1<<"; parent 2 is "<<idx2<<endl;
		int *parent1 = generation[idx1];
		int *parent2 = generation[idx2];
		printTour(parent1);
		printTour(parent2);
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
#ifdef DEBUG
			printf("# of cities : %d\n", CITY_N);
#endif
		}
	}
	cities = new city [CITY_N];
	while (fgets(line, 80, fp) != NULL && i < CITY_N) {
		sscanf(line, "%*d %lf %lf", &(cities[i].x), &(cities[i].y));
		++i;
	}
}


void printGeneration(int** generation) {
	for (int i = 0; i < size; i++) {
		printTour(generation[i]);
	}
}

void printCost(int* cost) {
	int totalCost = 0;
	for (int i = 0; i < size; i++) {
		totalCost += cost[i];
		cout<<"cost for the "<<i<<"th element is "<<cost[i]<<endl;
	}
	cout<<"average cost for this generation is "<< totalCost/size <<endl;
}
			
int main (int argc, char* argv[]) {
	srand(time(NULL));

	if (argc < 2) {
		printf("Usage: ./tsp inputFile\n");
		exit(-1);
	}

	readFile(argv[1]);

	int **generation = new int *[size];
	int *cost = new int [size];
	for (int i = 0; i < size; i++) {
		generation[i] = new int [CITY_N];
		for (int j = 0; j < CITY_N; j++) {
			generation[i][j] = j;
		}
		mutate(generation[i],initialMutationRate);
		cost[i] = calculate_cost(generation[i], cities);
	}

	printGeneration(generation);
	printCost(cost);

	for (int i = 0; i < generationNum; i++) {
		generation = evolveGeneration(generation, cost);
		printGeneration(generation);
		printCost(cost);
	}

/*
	int a[8] = {1, 2, 5, 6, 4, 3, 8, 7};
	int b[8] = {1, 4, 2, 3, 6, 5, 7, 8};
	int *c =  crossOver(a,b);
	for (int i = 0; i < 8; i++) {
		cout<<c[i]<<' ';
	}
	cout<<endl;
	if (!checkTour(c,8)) {
		cout<<"wrong tour"<<endl;
		exit(-1);
	}
*/
}
