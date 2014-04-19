#include "anneal.cu"

#define NUMBER_RUNS 10

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
		printf("Usage: ./tsp <inputFile>\n");
		exit(-1);
	}

	readFile(argv[1]);

	int *order = (int *) malloc(CITY_N * sizeof(int));
	float avgResult = 0.0f;
	double avgRuntime = 0.0f;
	int result[NUMBER_RUNS];

	for (int runs = 0; runs < NUMBER_RUNS; ++runs) {
		for (int i = 0; i < CITY_N; ++i)
			order[i] = i;

		//printCities(cities);

		Anneal *a = new Anneal();
		a->order(cities, order);
		avgResult += a->resultCost / (NUMBER_RUNS * 1.0f);
		avgRuntime += a->runtime / (NUMBER_RUNS * 1.0f);
		result[runs] = a->resultCost;
	}

	cout << endl << endl;
	cout << "Average Costs: " << avgResult << endl;
	cout << "Average Runtime: " << avgRuntime << endl;
	for (int i = 0; i < NUMBER_RUNS; i++) {
		cout<<"CUDA run result for "<<i<<" is "<<result[i]<<endl;
	}
	return 0;
}
