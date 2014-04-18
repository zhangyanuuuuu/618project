#include "anneal.cpp"

//#define FILENAME "../benchmark/kroA100.tsp"
#define NUMBER_RUNS 1
#define DEBUG

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
	cities = (struct city *) malloc(CITY_N * sizeof(struct city));
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
		printf("Usage: ./tsp inputFile\n");
		exit(-1);
	}
	readFile(argv[1]);

	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	int *order = (int *) malloc(CITY_N * sizeof(int));
	float avgResult = 0.0f;
	double avgRuntime = 0.0f;

	for (int runs = 0; runs < NUMBER_RUNS; ++runs) {
		for (int i = 0; i < CITY_N; ++i) {
				order[i] = i;
		}
		//broad cast
		assert(MPI_Bcast(order, CITY_N, MPI_INT, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS);

		//printCities(cities);

		Anneal *a = new Anneal(rank, size);
		a->order(cities, order);
		avgResult += a->resultCost / (NUMBER_RUNS * 1.0f);
		avgRuntime += a->runtime / (NUMBER_RUNS * 1.0f);
	}

	if (rank == ROOT) {
		cout << endl << endl;
		cout << "Average Costs: " << avgResult << endl;
		cout << "Average Runtime: " << avgRuntime << endl;
	}
	MPI_Finalize();
	return 0;
}
