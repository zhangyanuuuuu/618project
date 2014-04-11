#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>

//#define DEBUG
//#define PRINT

//---------------------------------------------------------------
// A simple annealer for doing Traveling Salesman problem
//
// To compile:  gcc tsp.cc -lm -o tsp
//
// UNIX command line args are:
//
// CITIES:     int, how many cities in TSP problem
// CLUSTERS:   int, group the cities into this many clusters
// SEED: int,  wants an integer (like 1234567) to
//             SEED the random number generation.
//             You get a different answer each time you
//             use a different seed to start the random
//             num generation.
// HOT: double,  an initial annealing HOT temperature, eg  20.0
// COOL: double, an annealing cooling rate, eg  0.90
// MOVESPER: int, how many moves per temperature, which is
//              MOVESPER * CITIES
//--------------------------------------------------------------


//--------------------------------------------------------------
//
// Simple changeable TSP control parameters 
//
//--------------------------------------------------------------

//  You can CHANGE  this stuff from the command line...
int        CITIES;                   // how any cities in TSP
int        SEED;                     // how we restart the random num gen
                                     // comes from UNIX command line   
double     HOT;                      // initial hot temperature 
double     COOL;                     // multiplier for how fast we cool Tnew = COOL*Told  
                                     // it's a fraction < 1
int        MOVESPER;                 // we do MOVESPER * nummods moves 
                                     
/* these params we assume just stay fixed */                                     

#define    TEMPS          40         // annealer must do at least this many temps 
#define    MINACCEPT      0.005      // must anneal until accept rate goes
                                     //  below this fraction 
#define    TOLERANCE      0.005      // tells how close the 3 consecutive
                                     // cost numbers seen at 3 temps
                                     // must be in order to say we're
                                     // frozen and can quit 


// annealing global variables that change their values  

int        tempcount;                // count how many temps we run


//--------------------------------------------------------------
//
// useful macros	
//
//--------------------------------------------------------------
#define abs(x)	(((x)>0)?(x):(-(x)))
#define max(x,y)  ( (x)>(y)?(x):(y) )
#define min(x,y)  ( (x)<(y)?(x):(y) )


//--------------------------------------------------------------
//
// ways to gen Random Vars with specific distributions 
//
//--------------------------------------------------------------


// returns a real number uniformly distributed on [0,1]. 
double uniformRandom( void )
{
    // get a random int from rand()
    // This is a number between 0 and RAND_MAX, inclusive
    int r = rand();   

    // turn it into a fraction in [0,1]  
    return (((double) r) / ((double) RAND_MAX));
}


// produces a random INTEGER in [imin, imax],
//  ie, including the endpoints imin and imax
int uniformRandomIntRange( int imin, int imax)
{
    double u = uniformRandom();
    int m = imin + (int)floor((double)(imax + 1 - imin)*u ) ;
    return  m ;
}

// produces a random Real in [rmin, rmax]
double uniformRandomRealRange( double rmin, double rmax)
{
    double u = uniformRandom();
    double m = rmin + (rmax - rmin)*u  ;
    return  m ;
}



//--------------------------------------------------------------
//
// Basic data structures & methods for the Trav Salesman 
//
//--------------------------------------------------------------

// Remember, there are CITIES cities, read from UNIX cmd line
// A city_struct is just a pair of coords, x and y
struct city_struct {
    double x, y;          
};   //should be dynamic--but this is easier...

city_struct* city;

// The solution is just an ordering of the cities, which
// means an ordering of the array index of the cities,
// a number in [0, CITIES-1].  So, tour[] holds this
// solution ordering
int *tour;

// the cost of the current solution, the sum of the lengths of
// all the distances between each pair of cities in the tour
double tourcost;


// compute length of the TSP link the "ends" at city i
double prevLinkCost(int i) {
   
    double delX, delY;
   int otherEnd; 
   // this fussy business handles the case when city i
   // is the very FIRST city in our tour, and the "previous"
   // city is the LAST city
   if( i==0) 
     otherEnd = CITIES-1;
   else otherEnd = i-1;

   delX = city[tour[i]].x - city[tour[otherEnd]].x;
   delX = delX*delX;
   delY = city[tour[i]].y - city[tour[otherEnd]].y;
   delY = delY*delY;

   return sqrt( delX + delY );
} 

// compute length of the TSP link the "starts" at city i
double nextLinkCost(int i) {
   
   double delX, delY;
   int otherEnd;
   // this fussy business handles the case when city i
   // is the very LAST city in our tour, and the "previous"
   // city is the FIRST city 
   if( i==CITIES-1) 
     otherEnd = 0;
   else otherEnd = i+1;

   delX = city[tour[i]].x - city[tour[otherEnd]].x;
   delX = delX*delX;
   delY = city[tour[i]].y - city[tour[otherEnd]].y;
   delY = delY*delY;

   return sqrt( delX + delY );
}  
 

// compute the cost of any given solution, FLAT, by looking
// at ALL links in the TSP tour. We only do this once
// per temperature during annealing, as a sort of sanity check.
double costTSP() {

    double cost = 0.0;
    for(int c=0; c<CITIES; c++) {
        // the tour length is the length of the link  "after" each city
        cost += nextLinkCost(c);
    }
    return cost;
}

// readin the file for location of cities, and then setup initial tour solution
void parse(char* fileName) 
{
	FILE* inFile;
	inFile = fopen(fileName, "r");

	ssize_t read;
	char* line = NULL;
	size_t len;
	//dump the first 6 lines
	for (int i = 0; i < 6; i++) {
		if ((read = getline(&line, &len, inFile)) != -1) {
#ifdef DEBUG
			printf("%s", line);
#endif
		} else {
			printf("input file format wrong\n");
			exit(-1);
		}
		if (i == 3) {
			char* pch = strtok(line, ":");
			pch = strtok(NULL, "\t\n");
			CITIES = atoi(pch);	
#ifdef DEBUG
			printf("# of cities : %d\n", CITIES);
#endif
		}
	}
	city = (city_struct*) malloc (CITIES*sizeof(city_struct));
	tour = (int*) malloc (CITIES*sizeof(int));

	// readin all the cities locations
	int c_idx;
	for (int i = 0; i < CITIES; i++) {
		fscanf(inFile, "%d %lf %lf", &c_idx, &city[i].x, &city[i].y);
#ifdef DEBUG
		assert(c_idx == i + 1);
#endif
	}
	if (line)
		free(line);
 
	// Now, do a stupid initial tour: the city visited
	// in the Nth place of the tour order is just
	// city N
	
	for(int n=0; n<CITIES; n++) {
		tour[n] = n;
	}

	// compute the cost of this stupid initial solution
  tourcost = costTSP();
}


//-------------------------------------------------------
//
// Annealing stuff
//
//-------------------------------------------------------


// ACCEPT CRITERION:
// This is the usual METROPOLIS accept criterion,
// returns 1 for accept, 0 reject
int accept(double deltac, double temperature)
{
    // annealing accept criterion 
    if( deltac < 0.0 )
        return( 1 );
    else {
        double pa = exp( (double)(-deltac)/temperature);
        if( uniformRandom() <= pa )
            return  1 ;
        else return  0 ;
    }
}

// MOVES:
// THis tries to swap the cities in the tour at
// locations  i, j and then
// returns the delta cost (distance) associated 
// with this swap.  This routine does NOT
// actually swap them; it just evaluates the
// change in distance due to the swap
double evalSwapFlat(int i,  int j)
{

    // first, swap the cities :
    //    city[i] <==> city[j] in tour order
    int tmp= tour[i];
    tour[i] = tour[j];
    tour[j] = tmp;

    // eval the new tour length.  This is the DUMB way
    // to do it -- only 2 cities moved, but we eval
    // the WHOLE tour cost, since we're being lazy
    double cost = costTSP();

    // now, UN-swap the cities  back
    //    city[j] <==> city[i]
    tmp= tour[i];
    tour[i] = tour[j];
    tour[j] = tmp;

    // return the delta in the tour cost
    return  (cost - tourcost);

}

// SMARTER version of MOVE
// This one tries to swap the cities in the tour at
// locations  i, j and then
// returns the delta cost (distance) associated 
// with this swap.  Again, routine does NOT
// actually swap them; it just evaluates the
// change in distance due to the swap
// BUT -- we do it incrementally, fast, this time
double evalSwap(int i,  int j)
{
double afterCost, beforeCost;

    // basic idea is simple:  just eval cost of the
    // links that we KNOW have changed

    // if i, j are not adjacent cities in the tour, then situation is
    // like this:
    //   0 1 2 3 ... i-1 i i+1 ... j-1 j j+1 ... n-1
    //               ---- ----     ---- ----
    // when we swap, we get
    //   0 1 2 3 ... i-1 j i+1 ... j-1 i j+1 ... n-1
    //               ---- ----     ---- ----
    // and these dashes are the only 4 links we have to worry
    // about computing and recomputing length of, before/after

    // case 1:  i, j not adjacent
    if ( abs(i-j) != 1
         || ( (i==0) && (j!=(CITIES-1)) )
         || ( (j==0) && (i!=(CITIES-1)) ) ) {

        // eval the cost of the 4 critical links
        beforeCost = prevLinkCost(i) + nextLinkCost(i) 
                   + prevLinkCost(j) + nextLinkCost(j);

        // now, swap the cities :
        //    city[i] <==> city[j] in tour order
        int tmp= tour[i];
        tour[i] = tour[j];
        tour[j] = tmp;

        // eval the changed tour length.  Look only at 4 critical links
        afterCost = prevLinkCost(i) + nextLinkCost(i) 
                   + prevLinkCost(j) + nextLinkCost(j);

        // now, UN-swap the cities  back
        //    city[j] <==> city[i]
        tmp= tour[i];
        tour[i] = tour[j];
        tour[j] = tmp;

        // return the delta in the tour cost
        return  (afterCost - beforeCost);
    }

    // case 2:  i, j are in fact adjacent
    // if i, j are ADJACENT cities in tour, then situation is
    // like this:
    //   0 1 2 3 ... i-1  i  j  j+1 ... n-1
    //               ----- -- ---- 
    // when we swap, we get
    //   0 1 2 3 ... i-1  j  i  j+1 ... n-1
    //               ----- -- ----- 
    // and these dashes are the only 3 links we have to worry
    // about computing and recomputing length of, before/after
    // Note also that we have to be careful to check if i,j are 1st or last
    // elements in the tour.  In that case, the "shared" link between
    // i, j wraps around from city 0 to to city (CITIES-1).
    if ( abs(i-j) == 1
         || ( (i==0) && (j==(CITIES-1)) )
         || ( (j==0) && (i==(CITIES-1)) ) ) {

        // we need to know which one is first in the tour
        int minCity = min(i,j);
        int maxCity = max(i,j);

        // eval the cost of the 3 critical links
        if( i!=0  && j!= 0) {
            // its the sub-case illustrated in comment above
            beforeCost = prevLinkCost(minCity) + nextLinkCost(minCity) 
                       + nextLinkCost(maxCity);
        }
        else  {
            // cities on the "ends" like:    i 1 2 3 ... n-3 n-2 j
            // we will swap to this:         j 1 2 3 ... n-3 n-2 i
            // eval the proper link cost, note its different
            beforeCost = prevLinkCost(minCity) + nextLinkCost(minCity) 
                       + prevLinkCost(maxCity);
        }

        // now, swap the cities :
        //    city[i] <==> city[j] in tour order
        int tmp= tour[i];
        tour[i] = tour[j];
        tour[j] = tmp;

        // eval ONLY the changed tour lengths.  Look only at 3 critical links
        minCity = min(i,j);
        maxCity = max(i,j);
        // eval the cost of the 3 critical links
        if( i!=0  && j!= 0) {
            // its the sub-case illustrated in first comment in this routine
            afterCost =  prevLinkCost(minCity) + nextLinkCost(minCity) 
                       + nextLinkCost(maxCity);
        }
        else  {
            // cities on the "ends" like:    i 1 2 3 ... n-3 n-2 j
            // we will swap to this:         j 1 2 3 ... n-3 n-2 i
            // eval the proper link cost, note its different
            afterCost =  prevLinkCost(minCity) + nextLinkCost(minCity) 
                       + prevLinkCost(maxCity);
        }

        // now, UN-swap the cities  back
        //    city[j] <==> city[i]
        tmp= tour[i];
        tour[i] = tour[j];
        tour[j] = tmp;

        // return the delta in the tour cost
        return  (afterCost - beforeCost);
    }
}


// THERMAL EQUILIBRIUM LOOP:
// this does the actual evolution of the TSP
// by annealing at a fixed tempt t passed in
// as an input.  startcost is also passed in
// so we can compute some statistics.
// acceptRatio (fraction of moves tried that were accepted)
// is passed back since the caller wants to use it
// to help decide if we are frozen yet
double annealAtTemperature( double t )
{
    // vars for  computing statistics of annealing
    // at the temperature, and some control parameters 
    double meanCost = 0.0;
    double varCost = 0.0;
    double currCost = costTSP();  // a little bit inefficient, but ok

    int attempts = MOVESPER*CITIES;
    double sample;
    int acceptCount = 0;
    double deltaCost;
    double totalDeltaCost = 0.0;

    // this is the main loop doing moves.
    // we do  'attempts' moves in all, then quit
    // at this temperature
    for(int m=1; m<=attempts; m++) {

        //generate a TSP swap move 
        int i = uniformRandomIntRange(0, CITIES-1);
        int j = uniformRandomIntRange(0, CITIES-1);
        if (j == i) {
           // just swap i with one of its neighbors
           if (i-1 >= 0 )
             j = i-1;
           else j = i+1;
        }
        deltaCost = evalSwap(i, j);

#ifdef DEBUG
 double debugCost = evalSwapFlat(i, j);
if( abs(deltaCost - debugCost) > 1e-9) {
  fprintf(stdout, "error T=%g  inc=%g  flat=%g  dif=%g\n",
    t, deltaCost, debugCost, abs(deltaCost-debugCost));
}
#endif

        // do we take this move?  run the
        // delta cost and temperature thru
        // metropolis criterion
        if( accept(deltaCost, t) ) {
            acceptCount++ ;
            totalDeltaCost += deltaCost;
            tourcost += deltaCost;

            // actually perform the swap     
           int tmp= tour[i];
           tour[i] = tour[j];
           tour[j] = tmp;  
        } 
        else {  // reject 
            deltaCost = 0;
        }

        // compute running sample statistics; these
        // stats are updated incrementally, with data taken
        // after EACH and every move.

        currCost = currCost + deltaCost;
        sample = (double)m ;

        // compute variance^2 now, since we need old mean value 
        if( sample <= 1.01 ) { // ie == 1 integer 
            varCost = 0.0;
	  }
	  else {
	      varCost = ( (sample - 2.0)/(sample - 1.0) ) * varCost;
	      varCost += ( (currCost - meanCost)
                         *(currCost - meanCost) )
                         /sample;
	  }

        // now update mean 
	  meanCost = meanCost +  (currCost - meanCost)/sample;
	
    } // for m=1 to attempts 

    // sanity checking:  lets make sure the cost we got using
    // the INCREMENTAL fast cost eval, evalSwap, gives
    // us the correct, current TOTAL cost.  Of course, since
    // these are floats, computed 2 different ways, we CANNOT
    // expect perfect equality.  They will different way down
    // in the significant digits, so we include this in the check.
    if ( abs( (currCost - costTSP()) / currCost ) > 1e-9 ) {
       fprintf(stderr, "OOPS:  T=%g currCost=%g  != flat costTSP=%g\n",
               t, currCost, costTSP() );
       exit(-1);
    }

    // compute accept rate at this temp, since we want to print it
    double accRate = (double)acceptCount / (double)attempts;

    // OK, we finished annealing all the moves at this
    // temperature.  Before we return, we print
    // some useful stuff so we can see the annealing
    // progress in the output file.
#ifdef PRINT
    fprintf(stdout, "run %d ", ++tempcount);
    fprintf(stdout, "T %g ", t);
    fprintf(stdout, "att %d ", attempts);
    fprintf(stdout, "acc %d ", acceptCount);
    fprintf(stdout, "accrate %g ", accRate);
    fprintf(stdout, "Cost %g ", currCost);
    fprintf(stdout, "meanC %g ", meanCost);
    fprintf(stdout, "varC %g ", varCost);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    //and, draw this current TSP solution

    // we return the accept rate for this temp
    return  accRate;
}

void print_tour(void)
{
	for (int i = 0; i < CITIES; i++) {
		printf("%d\n", tour[i]);
	}
	printf("The total cost is %lf\n", tourcost);
}

// ANNEAL:
// main annealing routine.  This controls
// the cooling procedure, calling the
// routine to anneal at each temperature
// and actually do the work.
void anneal(void)
{
double   cost3, cost2, costCurrent;
int      done, tempCount;
double   tol3, tol2, temp;
double   acceptRate;

    // set up the global control parameters for this
    // annealing run
    temp = HOT;
    tempCount = 0;
    cost3 = 999999999.0;
    cost2 = tourcost;
    costCurrent = cost2;

    // here is the temperature cooling loop of the annealer 
    done = 0;
    do {
        acceptRate = annealAtTemperature(temp);

        // we count the temps 
        tempCount++;

        // update the current cost after this temperature.
        costCurrent = tourcost; 

        // Now look at whether we are
        // frozen, or whether we should cool some more.
        // We basically just look at the last 2
        // temperatures, and see if the cost is not
        // changing much (thats the TOLERANCE test)
        // and if the we have done enough temperatures
        // (thats the TEMPS test), and if the accept
        // ratio fraction is small enough (that is the
        // MINACCEPT test).  If all are satisfied,
        // we quit.    
        tol3 = abs( (cost3 - cost2)/cost3 );
        tol2 = abs( (cost2 - costCurrent)/cost2 );

        if( tol3 < TOLERANCE
            && tol2 < TOLERANCE 
            && tempCount > TEMPS
            && acceptRate < MINACCEPT ){
            done = 1;
        }
        else {  // no, not frozen 

            // save the relevant info to test for frozen
            // after the NEXT temperature.
            cost3 = cost2;
            cost2 = costCurrent;

            // lower the temperature 
            temp = COOL * temp;
        }
    }while( ! done);   // go back and do annealing at next cooler temp 

    // OK, we're done lowering the temperature, its FROZEN.
    // BUT -- its common to do a little bit more work here.
    // We set to temp = 0 exactly (which will NEVER happen
    // if we do temp = COOL*temp.  This means we ONLY accept
    // downhill moves that improve the cost function.
    // That's OK, we are just trying to squeeze the last bits
    //  of improvement out of this process.  lets just do
    //  5 more of these "temperature=0" calls to the annealing loop
    for(int i=0 ; i<5; i++) {
        acceptRate = annealAtTemperature(0.0);
        tempCount++;
    }
}


//-----------------------------------------------------------------
//
// Main for the TSP */
//
//-----------------------------------------------------------------
int main(int argc, char *argv[]) 
{
 
    // To run:  tsp graphicsFile CITIES CLUSTERS SEED HOT COOL MOVESPER > outfile

    // go get command line arguments; complain if not right number 
    if(argc != 6) {
        fprintf(stderr, "usage: ./tsp inputFiles SEED HOT COOL MOVESPER \n");
        exit(-1);
    }

		// the first one is the file name for solving
		parse(argv[1]);

    // the second one is the seed for the random number
    // generator.  get it and convert it to int and save it
    SEED = atoi( argv[2] );
    
    // the fifth one is the initial value for the temperature. it's a double
    //   so you'd best type it with a decimal point.
    //   get it and convert it to int and save it
    HOT = atof( argv[3] );
    
    // the sixth one is the cooling rate (fraction < 1). it's a double
    //   so you'd best type it with a decimal point.
    //   get it and convert it to int and save it
    COOL = atof( argv[4] );
    
    // the seventh one is the moves-per-temp multiplier.  multiply this num
    // by the number of cities to be solved to get moves per temperature. 
    // it's an integer.
    //   get it and convert it to int and save it
    MOVESPER = atoi( argv[5] );

    // set up the random num stream
    srand(SEED);
    
    // do the TSP solution 
		clock_t t;
		t = clock();
    anneal();
		t = clock() - t;

		// print the tour
		print_tour();
		printf("Program took %d clocks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
    // that's it. 
		return 0;
}

