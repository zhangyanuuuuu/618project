#include "anneal.cpp"
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
	cities = new city [CITY_N];
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
	
	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	readFile(argv[1]); //read file seperately
		
	int *order = new int [CITY_N];
	float avgResult = 0.0f;
	double avgRuntime = 0.0f;
	int result[NUMBER_RUNS];

	for (int i = 0; i < CITY_N; ++i) {
		order[i] = i;
	}
	srand(time(NULL)+rank);

	for (int runs = 0; runs < NUMBER_RUNS; ++runs) {
		random_shuffle(order, order+CITY_N);

		//printCities(cities);

		Anneal a(rank, size);
		a.order(cities, order);
		avgResult += a.resultCost / (NUMBER_RUNS * 1.0f);
		avgRuntime += a.runtime / (NUMBER_RUNS * 1.0f);
		result[runs] = a.resultCost;
	}

	if (rank == ROOT) {
		cout << endl << endl;
		cout << "Average Costs: " << avgResult << endl;
		cout << "Average Runtime: " << avgRuntime << endl;
		for (int i = 0; i < NUMBER_RUNS; i++) {
			cout << "Run result for "<<i<<" is "<<result[i]<<endl;
		}
	}
	MPI_Finalize();
	
	delete [] order;
	delete [] cities;
	return 0;
}
