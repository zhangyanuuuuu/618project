#include    <stdio.h>
#include    <stdlib.h>
#include    <math.h>


/* ----------------------------------------
** Simple annealer for points in a grid, with
** fixed pads placed around the top, bottom,
** left, right edges of this grid.  Simple
** cooling schedule, no optimizations,
** only rather simple range limiting (which
** you can turn off by setting RANGERESTRICT = 1.0 
**
** To compile   gcc  760f96anneal.c  -lm -o anneal
** will do it.  Remember you need the
** math library for exp( ), for floor( ), etc.,
** so it will complain with out the -lm
**
** To run it   anneal netlistfile SEED HOT COOL MOVESPER RANGE > outfile
** will do it.  
**
** netlist (string):  It wants a file with a netlist.
**
** SEED (int) It wants an integer (like 1234567) to
** SEED the random number generation.
** You get a different answer each time you
** use a different seed to start the random
** num generation.
**
** HOT (float): It wants an initial annealing HOT temperature, eg  20.0
**
** COOL (float): It wants an annealing cooling rate, eg  0.90
**
** MOVESPER (int): It wants to know how many moves per temperature, which is
**  MOVESPER * number of objects.
**
** RANGE (float):  It wants to know how fast to shrink the range limit
**  window,  eg   0.98   per temperature.
** 
*/



/* simple changeable annealing control parameters 
**
**  You can CHANGE  this stuff from the command line...
**
*/
int        SEED;                     /* how we restart the random num gen
                                        comes from UNIX command line   */
double     HOT;                      /* initial hot temperature */ 
double     COOL;                     /* multiplier for how fast we cool Tnew = COOL*Told  
                                     ** it's a fraction < 1
                                     */
int        MOVESPER;                 /* we do MOVESPER * nummods moves */
double     RANGERESTRICT;            /* allowable distance of a swap
                                     ** reduces by this much after every
                                     ** temp.  It's a fraction < 1
                                     */  
                                     
/* these params we assume just stay fixed */                                     

#define    TEMPS          40         /* annealer must do at least this many temps */
#define    MINACCEPT      0.02       /* must anneal until accept rate goes
                                        below this fraction */
#define    TOLERANCE      0.01       /* tells how close the 3 consecutive
                                        cost numbers seen at 3 temps
                                        must be in order to say we're
                                        frozen and can quit */
#define    RANGEMIN       4          /* Smallest range limit we will
                                     ** shrink to using above shrink */


/* simple annealing glopbal variables that change their values  */

double     GlobalRangeLimit;         /* the actual current range limit,
                                     ** updated each temperature
                                     */
int        tempcount;                /* count how many temps we run*/




/* useful defs for status of nets and modules (a pad, or a gate) */
#define    UNDEFINED -1
#define    DEFINED    1
#define    PAD        0
#define    GATE       1

/*---------------------------------------------------*/
/* useful macros				     */


#define abs(x)	(((x)>0)?(x):(-(x)))
#define max(x,y)  ( (x)>(y)?(x):(y) )
#define min(x,y)  ( (x)<(y)?(x):(y) )


/* the netlist file we open to read out netlist */
FILE *nlfp;

/* nets data structure, with a simple limit of 10000 nets */
#define    NETS    10000
struct {
    int type;          /* is it defined or not, connected to pad or gate */
    int nummods;       /* how many moduled (gates, pads) on this net */
    int mods[1000];     /* IDs of the modules on this net */
    int len, templen;  /* how long is this net, permanent and temp vars */
    char touched;      /* have we seen this net already in a computation
                          of the length of this net? */
}  net[NETS];
int    numnets;        /* how many nets do we have in the whole design */

/*  modules (means a gate or a pad) data structure; 10000 max */
#define    MODS    10000
struct {
    int type;           /* is it defined or not, is it a pad or a gate */
    int numnets;        /* how many nets are attached to this module */
    int nets[50];        /* ID of those nets; max is 4 in this simple code*/
    int px, py;         /* x, y location of this placed object */
} mod[MODS];  
int    nummods;         /* how many total mods (pads + gates) we have */


/* chip surface modeled as a simple grid */
int        Xsize;              /* it's this many slots wide */
int        Ysize;              /* it's this many slots tall */
struct {
    int    occupied;           /* says if there is indeed an object here */
    int    mod;                /* tells ID of the object here */
} grid[100][100];              /* 100 x 100 default grid */

/* we reserve the edges of the grid for the pads.
** these defs just tell the X, Y numbers of the
** edges of the grid.  L means LEFT, R means RIGHT
** T means TOP, B means BOTTOM.   Note we reference
** the grid as   grid[x][y]
*/
#define    LPAD    (0)
#define    RPAD    (Xsize-1)
#define    TPAD    (Ysize-1)
#define    BPAD    (0)




/* Random number generator initialization */
#define A   147453245
#define C   226908347
#define M   1073741824
int     last_x;


/*----------------------------------------------------*/
/* ways to gen Random Vars with specific distributions */


/*
 * Simple random number generator based on the linear-congruential method
 * using parameters from example D, p 40, Knuth Vol 2. 
 *
 * uniform_rv() returns a real number uniformly distributed on [0,1]. This
 * version has the advantae that it should behave the same on different
 * machines, since the generator and starting point are explicitly
 * specified. 
 */
double
uniform_rv()
{
    /*
     * one small problem: the sequence we use can produce integers
     * larger than the word size used,  ie they can wrap around
     * negative.  We wimp out  on this matter and just make them
     * positive again. 
     */
    last_x = ((A * last_x) + C) % M;
     if (last_x < 0)
	last_x = -last_x; 
    return (((double)last_x) / ((double) M));
}

/* Initialize the random number stream */
void
init_rand (int seed)
{
    last_x = seed;
}

/* produces a random INTEGER in [imin, imax],
**  ie, including the endpoints imin and imax
*/
uniform_int_rv(imin, imax)
int	imin, imax;
{
double	u;
int	m;

	u = uniform_rv();
	m = imin + floor((double)(imax + 1 - imin)*u ) ;
	return( m );
}
	

/* go read the input netlist and build the
** data structures for the nets and the modules
*/
read_input()
{
int    done;
int    thismod, thisnet;
int    x, y;
char   type[3];
int    nnets;
int    n;

    /* init all nets, mods */
    for(thismod=0; thismod < MODS; thismod++) {
        mod[thismod].type = UNDEFINED;
    }
    nummods = 0;
    for(thisnet=0; thisnet<NETS; thisnet++) {
        net[thisnet].type = UNDEFINED;
        net[thisnet].touched = 0;
    }
    numnets = 0;

    /* read the size of the grid, Xsize then Ysize */
    fscanf(nlfp, "%d %d",  &Xsize, &Ysize);


    /* loop to read all the netlist data.
    *  each line look like either
    *    idnum  p  netid x y    (which says its a pad at x,y 
    *                            attached to net  netid)
    *  or
    *    idnum  g  numofnets  netid1  netid2, ... netidlast
    *                           (which says its a gate,
    *                            and how many nets, and the id of each)
    */
    done = 0;
    while( !done ) {
        fscanf(nlfp, "%d", &thismod);
        if( thismod < 0 ) {
            done = 1;
            break;
        }
        nummods++;

        fscanf(nlfp, "%1s", type);
        if(type[0]=='p') {

            /* its a pad */
            fscanf(nlfp, "%d %d %d", &x, &y, &thisnet);

            /* install it as a module */
            mod[thismod].type = PAD;
            mod[thismod].px = x;
            mod[thismod].py = y;
            mod[thismod].numnets = 1;
            mod[thismod].nets[0] = thisnet;
            
            /* install its net */
            if( net[thisnet].type == UNDEFINED ) {
                net[thisnet].type = DEFINED;
                net[thisnet].nummods = 1;
                net[thisnet].mods[0] = thismod;
    
                if( thisnet > numnets )
                    numnets = thisnet;
            }
            else {
                net[thisnet].mods[net[thisnet].nummods] = thismod;
                net[thisnet].nummods++;
            }
        }
        else if(type[0] == 'g') {
            /* its a regular gate */
            fscanf(nlfp, "%d", &nnets);

            /* install it as a module */
            mod[thismod].type = GATE;
            mod[thismod].numnets = nnets;
            
            /* get and install each of its nets  */
            for(n=0; n<nnets; n++) {
            
                fscanf(nlfp, "%d", &thisnet);

                /* add net to mod structure */
                mod[thismod].nets[n] = thisnet;
                
                /* install this net if its new */
                if( net[thisnet].type == UNDEFINED ) {
                    net[thisnet].type = DEFINED;
                    net[thisnet].nummods = 1;
                    net[thisnet].mods[0] = thismod;

                    if( thisnet > numnets )
                        numnets = thisnet;

                }
                else {
                    /* thisnet is not new, so just add thismod to net */

                    net[thisnet].mods[net[thisnet].nummods] = thismod;
                    net[thisnet].nummods++;
                }
            }
        }
        else {
            fprintf(stderr, "unknown module type: %s\n", type);
            exit(-1);
        }
    }

    fprintf(stdout, "init nummods %d  numnets %d\n", nummods, numnets);
    fflush(stdout);
        
}


/* a debugging routine to dump info on all the nets
** in the nets data structure
*/
dump_nets()
{
int    thisnet, m;

    for(thisnet=1; thisnet<=numnets; thisnet++) {
        printf("    net %d  len %d templen %d touched %d mods %d: ",
                thisnet, net[thisnet].len, net[thisnet].templen, 
                net[thisnet].touched, net[thisnet].nummods);
        for(m=0; m<net[thisnet].nummods; m++) {
            printf("%d ", net[thisnet].mods[m]);
        }
        printf("\n");
    }
}

/* a debugging routine to dump all infor the mod
** (gate or pad) data structure
*/
dump_mods()
{
int    thismod, m;

    for(thismod=1; thismod<=nummods; thismod++) {
        printf("    mod %d  type %d  x %d  y %d nets %d: ",
                thismod, mod[thismod].type, mod[thismod].px, mod[thismod].py,
                mod[thismod].numnets );
        for(m=0; m<mod[thismod].numnets; m++) {
            printf("%d ", mod[thismod].nets[m]);
        }
        printf("\n");
    }
}


/* start up the annealing by setting up
** a dumb random placemenet, computing
** an initial wirelength for all the nets
** in the netlist, and setting up the annealing
** parameters
*/
void
init_anneal()
{
int    x, y, n;
int    thisnet, thismod, m, left, right, top, bot;
int    len;

    /* init the random number generator;
    **  every time you want a different sequence of
    ** random nums, you have te re-seed this generator
    ** with a different starting SEED
     */
    init_rand(SEED);

    /* init the grid with empty's everywhere */
    for(x=0; x<Xsize; x++) {
        for(y=0; y<Ysize; y++) {
            grid[x][y].occupied = 0;
        }
    }

    /* do a stupid random placement */
    x = 1;
    y = 1;
    for(n=1; n<=nummods; n++) {

        if( mod[n].type == PAD ) {
            /* note that this location in grid is occupied
            ** but don't place pads -- they are already placed
            */
            grid[mod[n].px][mod[n].py].occupied = 1;
        }
        else {  /* it's a gate */
            grid[x][y].mod = n;
            grid[x][y].occupied = 1;
            mod[n].px = x;
            mod[n].py = y;

#ifdef DEBUG
           fprintf(stdout, "init place mod %d at (%d,%d)\n", n, x, y);
           fflush(stdout);
#endif

            x = x + 1;
            if( x >= (Xsize-1) ) {
                x = 1;
                y = y + 1;
                if( y > (Ysize-1) ) {
                    fprintf(stdout, "At (x,y) = (%d,%d), too many modules, not enough grid slots\n",
                                x, y);
                    fflush(stdout);
                    exit(-1);
                }
            }
        }
    } /* for each module */

#ifdef DEBUG
    fprintf(stdout, "OCCUPIED grid after random place \n");
    for(y=Ysize-1; y>=0; y--) {
        for(x=0; x<Xsize; x++) {
            fprintf(stdout, "%1d", grid[x][y].occupied);       
        }
        fprintf(stdout, "\n");
    }
#endif
        

    
    /* now that we have some placement, we can, 
    ** for each net, compute init wirelen
    */
    for(thisnet=1; thisnet<=numnets; thisnet++) {

        if(net[thisnet].type == UNDEFINED) {
            fprintf(stdout, "undefined net %d in init_anneal\n", thisnet);
            fflush(stdout);
            exit(-1);
        }

        /* do bounding box */
        left = RPAD+1;
        right = LPAD -1 ;
        top = BPAD - 1;
        bot = TPAD + 1;

        for(m=0; m<net[thisnet].nummods; m++) {
            thismod = net[thisnet].mods[m];

            /* update the bounding box for this net */
            if( mod[thismod].px < left )
                left = mod[thismod].px;
            if( mod[thismod].px > right )
                right = mod[thismod].px;
                
            if( mod[thismod].py < bot )
                bot = mod[thismod].py;
            if( mod[thismod].py > top )
                top = mod[thismod].py;

        } /* for each mod on thisnet */

    /* compute bounding box for this net */
    net[thisnet].len = ((top - bot) + (right - left) );
    
    } /* for each net */

#ifdef DEBUG
    printf("init_anneal: the total netlist\n");
    dump_mods();
    dump_nets();
#endif

    /* print basic initialization stuff for this run */
    fprintf(stdout, "init ");
    fprintf(stdout, "SEED %d ", SEED);
    fprintf(stdout, "HOT %f ", HOT);
    fprintf(stdout, "COOL %f ", COOL);
    fprintf(stdout, "TEMPS %d ", TEMPS);
    fprintf(stdout, "TOLERANCE %f ", TOLERANCE);
    fprintf(stdout, "MINACCEPT %f ", MINACCEPT);
    fprintf(stdout, "MOVESPER %d ", MOVESPER);
    fprintf(stdout, "RANGERESTRICT %g ", RANGERESTRICT);
    fprintf(stdout, "RANGEMIN %d ", RANGEMIN);
    fprintf(stdout, "\n");
    fflush(stdout);

}

/* the usual metropolis accept criterion,
** returns 1 for accept, 0 reject
*/
int
accept(deltac, temperature)
int	deltac;
double	temperature;
{
double	pa;

	/* annealing accept criterion */
	if( deltac <= 0 )
		return( 1 );
	else {
                pa = exp( (double)(-deltac)/temperature);
		if( uniform_rv() <= pa )
			return( 1 );
		else	return( 0 );
	}
}

/* try to swap the modules in the grid at
** locations  (x1 y1) and (x2 y2) and
** return the delta cost (wirelen) associated 
** with this swap.  This routine does NOT
** actually swap them; it just evaluates the
** change in wirelen due to the swap
*/
int
eval_swap(x1,y1,x2,y2)
int    x1,y1,x2,y2;
{
int    tempmod, tempocc;
int    delta_net_len;
int    left, right, top, bot;
int    mod1, mod2;
int    n, m, thisnet, thismod;

#ifdef DEBUG
        if(debug) {
        printf("\n\nstart eval_swap %d@(%d,%d)<> %d@(%d,%d)\n",
                grid[x1][y1].mod, x1, y1, grid[x2][y2].mod, x2,y2);
        }
#endif

    /* first, swap the modules :
    **    mod1@(x1,y1) <==> mod2@(x2,y2)
    */
    tempmod = grid[x1][y1].mod;
    tempocc = grid[x1][y1].occupied;
    grid[x1][y1].mod = grid[x2][y2].mod;
    grid[x1][y1].occupied = grid[x2][y2].occupied;
    grid[x2][y2].mod = tempmod;
    grid[x2][y2].occupied = tempocc;
    if(grid[x1][y1].occupied) {
        mod[grid[x1][y1].mod].px = x1;
        mod[grid[x1][y1].mod].py = y1;
    }
    if(grid[x2][y2].occupied) {
        mod[grid[x2][y2].mod].px = x2;
        mod[grid[x2][y2].mod].py = y2;
    }

    /* initialize the var for the delta in netlength */
    delta_net_len = 0;

    /* stick with original mod names, even though they just moved */
    mod1 = grid[x2][y2].mod;
    mod2 = grid[x1][y1].mod;

    /* look at every wire on module mod1 where mod1 is now (ie at x2,y2) */
    /* recall that we assume mod1 is not null, ie, mod had to be nonempty */
    for(n=0; n<mod[mod1].numnets; n++) {
    
        thisnet = mod[mod1].nets[n];
        net[thisnet].templen = 0;

#ifdef DEBUG
        if(debug) {
        printf("    net %d.....\n", thisnet);
        }
#endif
    
        /* set up to eval change in bounding box for nets */
        left = RPAD+1;
        right = LPAD -1 ;
        top = BPAD - 1;
        bot = TPAD + 1;

        for(m=0; m<net[thisnet].nummods; m++) {
            thismod = net[thisnet].mods[m];

            /* update the bounding box for this net */
            if( mod[thismod].px < left )
                left = mod[thismod].px;
            if( mod[thismod].px > right )
                right = mod[thismod].px;
                
            if( mod[thismod].py < bot )
                bot = mod[thismod].py;
            if( mod[thismod].py > top )
                top = mod[thismod].py;

#ifdef DEBUG
            if(debug) {
            printf("       mod %d at (%d,%d), h-p so far = %d\n",
                    thismod, mod[thismod].px, mod[thismod].py,
                    (top - bot) + (right - left)  );
            }
#endif


        } /* for each mod on thisnet */

        /* compute bounding box for this net */
        net[thisnet].templen = ((top - bot) + (right - left) ) ;
        net[thisnet].touched = 1;

        /* compute contribution to delta */
        delta_net_len += net[thisnet].templen - net[thisnet].len;

#ifdef DEBUG
        if(debug) {
        printf("    done: net %d newlen %d oldlen %d del %d\n", 
                thisnet, net[thisnet].templen, net[thisnet].len,
                net[thisnet].templen - net[thisnet].len);
        }
#endif

    }  /* for each net on mod 1 */

    /* look at every wire on module mod2 where mod2 is now at, ie x1,y1
    ** but NOT at wires we already touched,
    ** and NOT if mod2 is actually empty (which it can be!)
    */
    if(grid[x1][y1].occupied) {

        for(n=0; n<mod[mod2].numnets; n++) {
    
            thisnet = mod[mod2].nets[n];

#ifdef DEBUG
        if(debug) {
            printf("    net %d.....\n", thisnet);
        }
#endif

            /* have we seen it before? if so, next wire */
            if( net[thisnet].touched ) {

#ifdef DEBUG
                if(debug) {
                printf("    Touched before -- quit\n");
                }
#endif
                continue;    /* so we don't double count the length */
            }

            net[thisnet].templen = 0;

            /* set up to eval new bounding box for net */
            left = RPAD+1;
            right = LPAD -1 ;
            top = BPAD - 1;
            bot = TPAD + 1;

            for(m=0; m<net[thisnet].nummods; m++) {
                thismod = net[thisnet].mods[m];

                /* update the bounding box for this net */
                if( mod[thismod].px < left )
                    left = mod[thismod].px;
                if( mod[thismod].px > right )
                    right = mod[thismod].px;
                
                if( mod[thismod].py < bot )
                    bot = mod[thismod].py;
                if( mod[thismod].py > top )
                    top = mod[thismod].py;

#ifdef DEBUG
                if(debug) {
                printf("       mod %d at (%d,%d), h-p so far = %d\n",
                    thismod, mod[thismod].px, mod[thismod].py,
                    (top - bot) + (right - left)  );
                }
#endif

            } /* for each mod on thisnet */

            /* compute bounding box for this net */
            net[thisnet].templen = ((top - bot) + (right - left) ) ;
            net[thisnet].touched = 1;

            /* compute contribution to delta */
            delta_net_len += net[thisnet].templen - net[thisnet].len;

#ifdef DEBUG
        if(debug) {
        printf("    done: net %d newlen %d oldlen %d del %d touched %d\n", 
                thisnet, net[thisnet].templen, net[thisnet].len,
                net[thisnet].templen - net[thisnet].len, 
                net[thisnet].touched);
        }
#endif


        }  /* for each net on mod 2 */
    } /* if mod2 site is occupied */                


    /* now, swap the modules BACK:
    **    mod2@(x1,y1) <==> mod1@(x2,y2)
    ** and untouch all the nets whose length we updated
    */
    tempmod = grid[x1][y1].mod;
    tempocc = grid[x1][y1].occupied;
    grid[x1][y1].mod = grid[x2][y2].mod;
    grid[x1][y1].occupied = grid[x2][y2].occupied;
    grid[x2][y2].mod = tempmod;
    grid[x2][y2].occupied = tempocc;
    if(grid[x1][y1].occupied) {
        thismod = grid[x1][y1].mod;
        mod[thismod].px = x1;
        mod[thismod].py = y1;

        for(n=0; n<mod[thismod].numnets; n++) {
            thisnet = mod[thismod].nets[n];
            net[thisnet].touched = 0;
        }
    }
    if(grid[x2][y2].occupied) {
        thismod = grid[x2][y2].mod;
        mod[thismod].px = x2;
        mod[thismod].py = y2;

        for(n=0; n<mod[thismod].numnets; n++) {
            thisnet = mod[thismod].nets[n];
            net[thisnet].touched = 0;
        }
    }

    /* that's it.  We have the modules in their
    **  original (unswapped locations) and we
    **  have the change in netlength due to this 
    **  propsed swap
    */


    return( delta_net_len);

}

/* this is both a useful debugging routine, since
** you can call
** it to walk the entire netlist (hence the name "flat",
** since it's not incremental) and compute the total
** wirelength cost any time you need it to check
** against an evolving incremental cost 
** It's also an essential part of the annealing loop,
** since you need to call it at the beginning to get
** the starting cost at each new temperature.
*/
int
eval_flat_cost()
{
int    thisnet, thismod, m, left, right, top, bot;
int    len, total_len;
    
    total_len = 0;

    /* for each net, compute wirelen */
    for(thisnet=1; thisnet<=numnets; thisnet++) {

        if(net[thisnet].type == UNDEFINED) {
            fprintf(stdout, "undefined net %d in eval_flat_cost\n", thisnet);
            fflush(stdout);
            exit(-1);
        }

#ifdef DEBUG
        printf("evalflat: start net %d\n", thisnet);
        fflush(stdout);
#endif

        /* do bounding box */
        left = RPAD+1;
        right = LPAD -1 ;
        top = BPAD - 1;
        bot = TPAD + 1;

        for(m=0; m<net[thisnet].nummods; m++) {
            thismod = net[thisnet].mods[m];

            /* update the bounding box for this net */
            if( mod[thismod].px < left )
                left = mod[thismod].px;
            if( mod[thismod].px > right )
                right = mod[thismod].px;
                
            if( mod[thismod].py < bot )
                bot = mod[thismod].py;
            if( mod[thismod].py > top )
                top = mod[thismod].py;

#ifdef DEBUG
            printf("    mod %d at (%d,%d) has half-perim = %d\n",
                   thismod, mod[thismod].px, mod[thismod].py, 
                   (top - bot) + (right - left)   );
            fflush(stdout);
#endif

        } /* for each mod on thisnet */

        /* compute bounding box for this net */
        len = ((top - bot) + (right - left) );
        total_len += len;

#ifdef DEBUG
        printf("  total len is now: %d\n", total_len);
        fflush(stdout);
#endif
        
    }
    
    return(total_len);
}


/* this does the actual evolution of the placement
** by annealing at a fixed tempt t passed in
** as an input.  startcost is also passed in
** so we can compute some statistics.
** acceptratio (fraction of moves tried that were accepted)
** is passed back since the caller wants to use it
** to help decide if we are frozen yet
*/
int
anneal_at_temp(t, startcost, acceptratio)
double    t;
int       startcost;
double    *acceptratio;
{
int    m, got;
int    acceptcount, attempts;
double  sample, mean_cost, var_cost, current_cost;
int    delta_cost, total_delta_cost;
int    x1, y1, x2, y2, xmin, xmax, ymin, ymax;
int    tempmod, tempocc;
int    n, thismod, thisnet;
int    accstatus, save_delta, currentcost;

#ifdef DEBUG
    printf("entering anneal_at_temp\n");
    fflush(stdout);
#endif

    /* vars for  computing statistics of annealing
    ** at the temperature, and some control parameters 
    */
    mean_cost = 0.0;
    var_cost = 0.0;
    current_cost = (double) startcost;
    attempts = MOVESPER*nummods;
    acceptcount = 0;
    total_delta_cost = 0;

    /* this is the main loop doing moves.
    ** we do  'attempts' moves in all, then quit
    ** at this temperature
    */
    for(m=1; m<attempts; m++) {

        /* generate a legit swap move */

        /* get a source gate.  
        ** NOTE  we never pick a pad to swap, since they
        ** are fixed, so we make sure not to grab any
        ** placed objects on top, bottom, left or right
        ** rows of the gris, since these are the pads
        */
        got = 0;
        while( !got ){
            x1 = uniform_int_rv(1, RPAD-1);
            y1 = uniform_int_rv(1, TPAD-1);

            /* cannot be empty slot */
            if(grid[x1][y1].occupied == 1)
                got = 1;
        }

#ifdef DEBUG
        printf("move %d just got source (%d, %d)\n", m, x1, y1);
        fflush(stdout);
#endif

        /* figure out what our simple global range limiter
        ** says for limits for where we can swap to
        ** in the grid array.  GlobalRangeLimit is a number
        ** for how far we are allowed to look in X, Y
        ** for another module to swap with
        */
        xmin = max(1, x1 - (int)(GlobalRangeLimit));
        xmax = min(RPAD-1, x1 + (int)(GlobalRangeLimit));
        ymin = max(1, y1 - (int)GlobalRangeLimit);
        ymax = min(TPAD-1, y1 + (int)GlobalRangeLimit);
        
        /* find a target cell--can be empty now
        ** but it has to be different that the
        ** source cell we just got above
        */
        got = 0;
        while( !got ) {
            x2 = uniform_int_rv(xmin, xmax);
            y2 = uniform_int_rv(ymin, ymax);
            if( (x2!=x1 || y2!=y1) )
                got = 1;
        }
#ifdef DEBUG
        printf("move %d just got target (%d, %d)\n", m, x2, y2);
        fflush(stdout);
#endif

        /* OK, we have 2 different slots in the grid,
        ** let's see what the delta cost is if we
        ** swap their contents
        */
        delta_cost = eval_swap(x1,y1,x2,y2);
        save_delta = delta_cost;

#ifdef DEBUG
       printf("swap (%d %d) %d with  (%d %d) %d, accept(delta = %d, T = %g) is %d\n",
       x1, y1, (grid[x1][y1].occupied)?grid[x1][y1].mod:(-1),
       x2, y2, (grid[x2][y2].occupied)?grid[x2][y2].mod:(-1),
       delta_cost, t,
       accept(delta_cost, t));
#endif

        /* do we take this move?  run the
        ** delta cost and temperature thru
        ** metropolis criterion
        */
        if( accept(delta_cost, t) ) {
            acceptcount++ ;
            accstatus = 1;
            total_delta_cost += delta_cost;

#ifdef DEBUG
            printf("t %f move %d %d@(%d,%d)->%d@(%d,%d) del %d  ACCEPT\n",
                    t, m, grid[x1][y1].mod, x1, y1, 
                    grid[x2][y2].mod, x2, y2, delta_cost);

            fflush(stdout);
#endif

            /* actually perform the swap */
            tempmod = grid[x1][y1].mod;
            tempocc = grid[x1][y1].occupied;
            grid[x1][y1].mod = grid[x2][y2].mod;
            grid[x1][y1].occupied = grid[x2][y2].occupied;
            grid[x2][y2].mod = tempmod;
            grid[x2][y2].occupied = tempocc;
            if(grid[x1][y1].occupied) {
                mod[grid[x1][y1].mod].px = x1;
                mod[grid[x1][y1].mod].py = y1;
            }
            if(grid[x2][y2].occupied) {
                mod[grid[x2][y2].mod].px = x2;
                mod[grid[x2][y2].mod].py = y2;
            }

            /* When we did the eval of this swap, that eval
            ** eval routine tagged each net with a temporary
            ** netlength, just in case we accepted the move.
            ** So, now we go back and we make that templen 
            ** of each net affected by this swap the REAL len 
            */
            if(grid[x1][y1].occupied) {
                thismod = grid[x1][y1].mod;
                for(n=0; n<mod[thismod].numnets; n++) {
                    thisnet = mod[thismod].nets[n];
                    net[thisnet].len = net[thisnet].templen;
                }
            }
            if(grid[x2][y2].occupied) {
                thismod = grid[x2][y2].mod;
                for(n=0; n<mod[thismod].numnets; n++) {
                    thisnet = mod[thismod].nets[n];
                    net[thisnet].len = net[thisnet].templen;
                }
            }

#ifdef DEBUG
        printf("After ACCEPTED swap...\n");
        dump_nets();
        dump_mods();
        printf("\n\n\n");
#endif


        } /* if accept */
        else {  /* reject ! */

#ifdef DEBUG
            printf("t %f move %d %d@(%d,%d)->%d@(%d,%d) del %d  REJECT\n",
                    t, m, grid[x1][y1].mod, x1,y1,
                    grid[x2][y2].mod, x2, y2, delta_cost);
            fflush(stdout);
#endif

            accstatus = 0;
            delta_cost = 0;

        }

#ifdef DEBUG
        /* ultra-sanity check --
        ** go thru EVERY net in the design and compute
        ** the total new cost.  Expensive to do 
        ** but occasionally helpful when your
        ** cost function or moves just are not
        ** working right
        */
        currentcost = eval_flat_cost();
        if(startcost + total_delta_cost != currentcost) {

            /* error somewhere */
            printf("ERROR: temp %f  move %d  accstatus %d\n", t, m, accstatus);
            printf("       startcost %d + total_delta_cost %d != currentcost %d\n",
                    startcost, total_delta_cost, currentcost);
            printf("       move %d: %d@(%d,%d)->%d@(%d,%d) del was %d \n",
                    m, grid[x2][y2].mod, x1,y1,
                    grid[x1][y1].mod, x2, y2, save_delta);
            printf("netlist status NOW:\n");
            dump_mods();
            dump_nets();
            fflush(stdout);
            exit(-1);
        }
#endif

        /* compute running sample statistics */
	/* stats updated incrementally, with data taken
	** after each and every move.
	** Formulas from Kobayashi book "Modeling and Analysis"
	** page 325 prob 5.1
	*/

        current_cost = (double)(current_cost) + (double)(delta_cost);
        sample = (double)m ;


        /* compute variance^2 now, since we need old mean value */
        if( sample <= 1.01 ) { /* ie == 1 integer */
            var_cost = 0.0;
	}
	else {
	    var_cost = ( (sample - 2.0)/(sample - 1.0) ) * var_cost;
	    var_cost += ( (current_cost - mean_cost)
                         *(current_cost - mean_cost) )
                         /sample;
	}

        /* now update mean */
	mean_cost = mean_cost +  (current_cost - mean_cost)/sample;

	
    }/* for m=1 to attempts */

    /* OK, we finished annealing all the moves at this
    ** temperature.  Before we return, we print
    ** some useful stuff so we can see the annealing
    ** progress int he output file.
    */
    fprintf(stdout, "run %d ", ++tempcount);
    fprintf(stdout, "T %g ", t);
    fprintf(stdout, "att %d ", attempts);
    fprintf(stdout, "acc %d ", acceptcount);
    fprintf(stdout, "accratio %g ", (double)acceptcount / (double)attempts);
    fprintf(stdout, "cost %g ", current_cost);
    fprintf(stdout, "meancost %g ", mean_cost);
    fprintf(stdout, "GlobalRangeLim %g ", GlobalRangeLimit);
    fprintf(stdout, "varcost %g ", var_cost);
    fprintf(stdout, "\n");
    fflush(stdout);

    /* a hack: shouldn't have to pass this back thru ptr */
    *acceptratio =  (double)acceptcount / (double)attempts;

    return(total_delta_cost);
}


/* main annealing routine.  This controls
** the cooling procedure, calling the
** routine to anneal at each temperature
** and actually do the work.
*/
void
anneal()
{
int   cost3, cost2, costcurrent;
int   done;
double tol3, tol2, temp;
double acceptratio;
int   delta_cost;

    /* set up the global control parameters for this
    ** annealing run
    */
    temp = HOT;
    tempcount = 0;
    cost3 = 999999999;
    cost2 = eval_flat_cost();
    costcurrent = cost2;
   
    /* initial range limit on swaps is ~ 1/2 the distance across
    ** the larger dimension of the grid representing the chip
    */
    GlobalRangeLimit =  ( max(RPAD-LPAD+1, TPAD-BPAD+1) / 2 );


    fprintf(stdout, "init StartingCost %d \n", costcurrent);
    fflush(stdout);

    /* here is the temperature cooling loop of the annealer */
    done = 0;
    do {
        delta_cost = anneal_at_temp(temp, costcurrent, &acceptratio);

        /* get the current cost after this temperature.
        ** we do this the simple but expensive way, by
        ** using a full eval of all nets
        ** Could also just do this as cost2 + delta_cost,
        ** which is the incremental cost change added
        ** to cost at last temp.
        */
        costcurrent = eval_flat_cost();   


        /* Sanity check: see if the real current cost
        ** is the same as the one the annealing at
        ** this temperature computed incrementally,
        ** one move at a time
        */ 
        if( (cost2 + delta_cost) != costcurrent ) {
            fprintf(stdout, "ERROR: cost2 %d + delta %d != costcurrent %d\n",
                    cost2, delta_cost, costcurrent);
            dump_mods();
            dump_nets();
            fflush(stdout);
            exit(-1);
        } 


        /* OK, if we got here the cost function is working
        ** fine. We can now look at whether we are
        ** frozen, or whether we should cool some more.
        ** We basically just look at the last 2
        ** temperatures, and see if the cost is not
        ** changing much (thats the TOLERANCE test)
        ** and if the we have done enough temperatures
        ** (thats the TEMPS test), and if the accept
        ** ratio fraction is small enough (that is the
        ** MINACCEPT test).  If all are satisfied,
        ** we quit.
        */    
        tol3 = ((double)cost3 - (double)cost2)/ (double)cost3;
        if(tol3 < 0) tol3 = -tol3;
        tol2 = ((double)cost2 - (double)costcurrent)/ (double)cost2;
        if(tol2 < 0) tol2 = -tol2;

        if( tol3 < TOLERANCE
            && tol2 < TOLERANCE 
            && tempcount > TEMPS
            && acceptratio < MINACCEPT ){
            done = 1;
        }
        else {  /* no, not frozen */

            /* save the relevant info to test for frozen
            ** after the NEXT temperature.
            */
            cost3 = cost2;
            cost2 = costcurrent;

            /* lower the temperature */
            temp = COOL * temp;

            /* shrink the range limit for swaps 
            ** make sure we don't shrink it TOO much
            */
            GlobalRangeLimit = RANGERESTRICT *  GlobalRangeLimit;
            GlobalRangeLimit = max( RANGEMIN, GlobalRangeLimit );
        }
    }while( ! done);   /* go back and do annealing at next cooler temp */

}

/* routine to be used at the very end of the placement
** to dump some useful info:  the actual final placemet
** x, y of each  module, including both gates and pads
*/
print_results()
{
int    thismod;

    printf("\n\n\n");
    for(thismod=1; thismod<=nummods; thismod++) {
        printf("%d %c %d %d\n", thismod, 
                (mod[thismod].type == GATE)?('g'):('p'),
                mod[thismod].px, mod[thismod].py);
    }

}

/* main for the placer */
void
main(argc, argv) 
int    argc;
char    *argv[];
{
 
    /* go get command line arguments; complain if not right number */
    if(argc != 7) {
        fprintf(stderr, "usage: anneal netlistfile SEED HOT COOL MOVESPER RANGE > outfile\n");
        exit(-1);
    }

    /* the first one is supposed to be the file with the
    ** netlist.  Go open it, complain and die if we can't
    */
    if( (nlfp=fopen(argv[1], "r"))==NULL) {
        fprintf(stderr, "cannot open %s\n", argv[1]);
        exit(-1);
    }
   
    /* the second one is the seed for the random number
    ** generator.  get it and convert it to int and save it
    */
    SEED = atoi( argv[2] );
    
    /* the third one is the initial value for the temperature. it's a double
    ** so you'd best type it with a decimal point.
    **   get it and convert it to int and save it
    */
    HOT = atof( argv[3] );
    
    /* the fourth one is the cooling rate (fraction < 1). it's a double
    ** so you'd best type it with a decimal point.
    **   get it and convert it to int and save it
    */
    COOL = atof( argv[4] );
    
    /* the fifth one is the moves-per-temp multiplier.  multiply this num
    ** by the number of gates to be placed to get moves per temperature. 
    ** it's an integer.
    **   get it and convert it to int and save it
    */
    MOVESPER = atoi( argv[5] );
    
    /* the sixth one is the range limiter multiplier (fraction<1).  
    ** window in which we will try swaps shrinks by this much each temp. 
    ** it's a double -- best type it with a decimal point.
    **   get it and convert it to int and save it
    */
    RANGERESTRICT = atof( argv[6] );


    /* go read the netlist */
    read_input();

    /* go set up annealing */
    init_anneal();

    /* do the placement */
    anneal();

    /* print the final placement */
    print_results();

    /* that's it. */
}
