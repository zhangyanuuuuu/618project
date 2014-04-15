#include "anneal_static.cpp"

#define FILENAME "../benchmark/a280.tsp"
#define NUMBER_RUNS 1

void readFile(struct city *cities)
{
	FILE *fp;
	char line[80];
	int i = 0;

	fp = fopen(FILENAME, "rt");
	while (fgets(line, 80, fp) != NULL) {
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

int main()
{
	struct city *cities = (struct city *) malloc(CITY_N * sizeof(struct city));
	int *order = (int *) malloc(CITY_N * sizeof(int));
	float avgResult = 0.0f;
	double avgRuntime = 0.0f;

	readFile(cities);

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
