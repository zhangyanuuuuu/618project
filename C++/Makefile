all: tsp tsp_static tsp_swap

tsp: anneal.cpp tsp.cpp
	g++ tsp.cpp -O3 -o tsp

tsp_swap: anneal_swap.cpp tsp_swap.cpp
	g++ tsp_swap.cpp -O3 -o tsp_swap

tsp_static: anneal_static.cpp tsp_static.cpp
	g++ tsp_static.cpp -O3 -o tsp_static

clean:
	rm tsp tsp_static tsp_swap
