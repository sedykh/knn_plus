//Feature Selection algorithms: 
//Genetic Algorithm (GA), Simulated Annealing (SA), Ant Colony Optimization (ACO), Particple Swarm Optimization (PSO)

#if !defined(FSEL_ALG_CLASS)
#define FSEL_ALG_CLASS

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "core.h"			//basic types


//gene flags
#define	GENE_INTEGER		true
#define	GENE_REAL			false

#define UNDEFINED			ZERO

//modes for selection of parents
#define RANK_BASED			10				//some thresshold is used to pick up small elite group of invidiuals which then randomly intercourse :)
#define ROULETTE			20				//probablity for the individual to be selected as a parent is proportional to its "fitness" factor
#define TOURNAMENT			30				//best individual of a small, randomlly picked subset of strings from the current population is selected as a parent

//modes for the crossover
#define ONE_POINT			40
#define	TWO_POINT			50
#define UNIFORM				60
#define TWO_OR_ONE_POINT	70
#define CRS_ALIGN_HEADS		1				//only for varying-size GA

//-------------------- Default Settings --------------------//
#define		DEFAULT_POOL						100			//100 genes
#define		DEFAULT_IND_SIZE					20			//20 genes per individual
#define		DEFAULT_MAX_IND_SIZE				40			//used only if size of individuals can be different
#define		DEFAULT_MIN_IND_SIZE				5			//used only if size of individuals can be different
#define		DEFAULT_SIZE						200			//default size of the population
#define		DEFAULT_MAXGEN						500			//default maximum number of generations in optimization procedure
#define		DEFAULT_MUTATION_RATE				0.05		//%
#define		DEFAULT_CROSSOVER_CHANCE			0.8			//%
#define		DEFAULT_TEMPERATURE					300			//K
#define		DEFAULT_MIN_DIVERSION				0.001
#define		DEFAULT_STABILITY_RANGE				10
#define		DEFAULT_IDEAL_FITNESS				INVALID		//means it's not used
#define		DEFAULT_VARIATIONOFSIZE_MODE		false			
#define		DEFAULT_TEMPERATURE_CHANGE_MODE		false
#define		DEFAULT_ELITISM_MODE				false
#define		DEFAULT_GROUP_SIZE					5
#define		DEFAULT_ELITEPART_SIZE				0.1			//%
#define		DEFAULT_CROSSOVER_MODE				TWO_OR_ONE_POINT
#define		DEFAULT_PARENTSELECTION_MODE		RANK_BASED
//--------------------                  --------------------//

//Core classes to be used in all feature selection algorithms
class Gene
{
public:
	Gene();
	Gene(UNSIGNED_4B_TYPE);
	Gene(REALNUM_TYPE, REALNUM_TYPE);
	Gene(const Gene &G);

	virtual ~Gene();

	void MutateGene();
	

	//gene's body
	UNSIGNED_4B_TYPE	GD;		//body for the gene in case of the limited pool of genes
	UNSIGNED_4B_TYPE	N;		//size of the pool of genes if it's limited
	
	REALNUM_TYPE		GR;		//body for the gene in case of unlimited pool of genes
	REALNUM_TYPE		A, B;	//range of variation of the genes in case of unlimited pool

	//switcher
	bool gen_pool_limited;		//defines whether the number of gens is limited
};

//comparison operators for class Gene
bool operator == (const Gene &, const Gene &);
bool operator != (const Gene &, const Gene &);
bool operator <  (const Gene &, const Gene &);
bool operator >  (const Gene &, const Gene &);
bool operator >= (const Gene &, const Gene &);
bool operator <= (const Gene &, const Gene &);

class Individual
{
public:
	Individual();
	Individual(UNSIGNED_4B_TYPE);
	Individual(UNSIGNED_4B_TYPE, UNSIGNED_4B_TYPE, bool);
	Individual(UNSIGNED_4B_TYPE, REALNUM_TYPE, REALNUM_TYPE, bool);
	Individual(const Individual &I);

	virtual ~Individual();	

	REALNUM_TYPE Compare(Individual &);						//calculates disimilarity of two individual solutions
	bool CompareExact(Individual &);						//for integer genes only, checks if soluations are the same
	void	MutateIndividual();								//member function that performs mutation
	void	EvaluateIndividual(GENERIC_POINTER DFunct);		//calls  D-Function to evaluate individual
	void	SetGenesTo(UNSIGNED_4B_TYPE = 0);
	void	SetRGenesTo(REALNUM_TYPE);
	
	apvector<class Gene> body;

	REALNUM_TYPE quality;									//fittness of the individual, 															
															//DFunction is used to calculate it. 
															//Greater value implies better fit individual.
															// Should be always positive. -1 value is used as a flag of undefined state!

	bool SameGeneShouldNotOccur;							//restricts an individual's genom to have only 1 gene of each type
};


/*Genetic Algorithm (GA)
//	Configure()				use it to set up all parameters in a way you need
//							and do not forget to provie Dfunction!
//	InsertIndividual();		use it to insert some favourite
//							individuals into population before 
//							optimization starts.
//	PerformOptimization();	use it to perform whole optimization.
//							both parameters can be NULL, but if not, then
//							first parameter-function will be called once per 
//							cycle in an optimization loop. Its return code values 
///							can abort/resume or stop optimization. 
//							Second parameter-function is used for output
//							Look at the implementation in cpp-file for 
//							specific details and format of these functions
//
//  NOTE:
//							Sharing, Crowding, Parallelism and Hybridisation
//							modes extend the basic scheme to fight 
//							premature convergence. ClearDuplIndividuals() does Crowding.
//							
//							SHARING implies that if some individual is very cose
//							to another one, then it will be either forced to undergo 
//							a mutation or its fitness value will be lowered.
//
//							CROWDING mode works in a way that if newly generated
//							individual is too close to one of existing in a pool, 
//							it will replace it.
//						
//							PARALLELISM employes several small sub-populations instead 
//							of one big one.  Those small subpopulation can exchange genetic 
//							material at some timepoints. Such a mode limits the influence 
//							of a very good individual to the borders of one subpopulation.
//
//							HYBRIDISATION consists of quick-and-dirty local optimization
//							for each trial individual in the evaluation phase instead of
//							evaluating its fitness factor. Then optimized intividual will
//							replace the initial one. In fact it can be done by proper
//							coding of the D-Function.
*/

class GenAlgorithm  
{
public:
	GenAlgorithm();
	virtual ~GenAlgorithm();

	bool InsertIndividual(Individual);
	void Configure(
		//first parameters that describe individual go
							 bool 				= true,
							 bool 				= GENE_INTEGER,
							 bool 				= DEFAULT_VARIATIONOFSIZE_MODE,
							 UNSIGNED_4B_TYPE 	= DEFAULT_POOL,
							 REALNUM_TYPE 		= 0,
							 REALNUM_TYPE 		= 1,
							 UNSIGNED_4B_TYPE 	= DEFAULT_IND_SIZE,
							 UNSIGNED_4B_TYPE 	= DEFAULT_MAX_IND_SIZE,
							 UNSIGNED_4B_TYPE 	= DEFAULT_MIN_IND_SIZE,

							 //parameters of the algorithm
							 UNSIGNED_4B_TYPE 	= DEFAULT_SIZE, 
							 UNSIGNED_4B_TYPE 	= DEFAULT_MAXGEN,							 
							 REALNUM_TYPE 		= DEFAULT_MUTATION_RATE,
							 REALNUM_TYPE 		= DEFAULT_CROSSOVER_CHANCE,							 
							 REALNUM_TYPE 		= DEFAULT_MIN_DIVERSION,
							 UNSIGNED_4B_TYPE	= DEFAULT_STABILITY_RANGE,
							 REALNUM_TYPE 		= DEFAULT_IDEAL_FITNESS,
							 GENERIC_POINTER DFUNCT		= NULL,
							 //modes
							 UNSIGNED_4B_TYPE 	= DEFAULT_CROSSOVER_MODE,
							 UNSIGNED_4B_TYPE 	= DEFAULT_PARENTSELECTION_MODE, 							 
							 bool 				= DEFAULT_ELITISM_MODE,
							 //some additional parameters needed for some modes
							 UNSIGNED_4B_TYPE 	= DEFAULT_GROUP_SIZE,
							 REALNUM_TYPE 		= DEFAULT_ELITEPART_SIZE);

	void PerformOptimization(GENERIC_POINTER, GENERIC_POINTER);		//performs whole optimization, calls some user-defined function as well
	bool PrepareOptimization();						//prepares everything for the start
	bool PerformOneCycle(GENERIC_POINTER);					//performs one cycle of the optimization (1 generation)
		
	UNSIGNED_4B_TYPE GetCurrGener();				//Encapsulation of currGener field
	void SetCurrGener(UNSIGNED_4B_TYPE);

private:
	bool CheckSettings();
	void PickParent(UNSIGNED_4B_TYPE &);
	void CrossOver(UNSIGNED_4B_TYPE, UNSIGNED_4B_TYPE, UNSIGNED_4B_TYPE &);		//peforms crossover of two parents
	void ProduceNewGeneration();
	void EvaluateIndividuals(apvector<Individual> &);

public:
	apvector<class Individual> population;
	apvector<class Individual> next_generation;

private:
	UNSIGNED_4B_TYPE currGener;					//current generation ID
	GENERIC_POINTER	DFUNCTION;					//pointer to D-Function
	UNSIGNED_4B_TYPE p_size;					//Size of the population

	REALNUM_TYPE mutation_rate;					//probability for the mutation in a child,			should be <=1, usual rate is 0.001 - 0.05
	REALNUM_TYPE crossover_chance;				//probability for the crossover	between parents,	should be <=1, if crossover did not take place, then children will be identical to parents. usual range is 0.5-0.9
	
	REALNUM_TYPE Diversion;						//difference in "fitness" parameter between current generation's best individual and that one from the previous generation, when becomes negligible, stable conditions are reached.	
	UNSIGNED_4B_TYPE Range;						//work variable for StableRange and MinDiversion mode

	//the termination criteria
	UNSIGNED_4B_TYPE MaxNGener;					//maximal number of generations
	REALNUM_TYPE MinDiversion;					//if difference of"fitnesses" of the previous and the current generations is less than this value, then optimization is finished.
	UNSIGNED_4B_TYPE StabilityRange;			//span of generations through which  MinDiversion should be maintained for the STOP of optimization
	REALNUM_TYPE IdealFitness;					//if "fitness" factor of a certain individual in the population is greater or equal than this, optimisation is finished

	//various modes of work		
	bool IndSizeConstant;						//specifies whether the number of genes in each individual is constant or not.
	bool Elitism;								//best individuals of the current generation are atomatically go into next generation (ElitePart is used for that)

	//Parents Selection modes
	UNSIGNED_4B_TYPE	PSelMode;				//RANK_BASED, ROULETTE, TOURNAMENT
	UNSIGNED_4B_TYPE	GroupSize;				//size of the randomly picked group of individuals amongst which the best one is picked as a parent. Used if PSelMode = TOURNAMENT_BASED 
	UNSIGNED_4B_TYPE	MaxIndSize;				//maximal size of the individual (used only when IndSizeConstant is ON)
	UNSIGNED_4B_TYPE	MinIndSize;				//minimal size of the individual (used only when IndSizeConstant is ON)
	REALNUM_TYPE		ElitePart;				//this parameter is used for picking the best individuals from the population percentage of "fittest" individual	should be <=1. Used if PSelMode = RANK_BASED

	UNSIGNED_4B_TYPE	CrossOverMode;			//ONE_POINT, TWO_POINT, UNIFORM, TWO_OR_ONE_POINT
};

//Simulated Annealing (SA)
class SimAnneal
{
public:
	SimAnneal();
	virtual ~SimAnneal();

	void Configure(
	//parameters that describe individual solution
	bool ifCanSameG = false,
	UNSIGNED_4B_TYPE nGPool = DEFAULT_POOL,
	REALNUM_TYPE rtGs = 0, //real value, start margin
	REALNUM_TYPE rtGe = 1, //real value, end margin
	UNSIGNED_2B_TYPE nSize = DEFAULT_IND_SIZE,

	//optimization parameters
	bool ifCM = false,
	UNSIGNED_2B_TYPE nB = 1,		//number of best models to store
	GENERIC_POINTER gpFunc = NULL,	//score function
	REALNUM_TYPE rtIS = 1.0,		//ideal value of score-function

	REALNUM_TYPE rtPM = 0.1,		//probability of mutation per gene of solution
	UNSIGNED_4B_TYPE nI = 100,		//number of mutation trials to do before changing temperature
	REALNUM_TYPE rtT_0 = 3, REALNUM_TYPE rtT_END = -6, 
	REALNUM_TYPE rtT_CONV = -6, REALNUM_TYPE rtK =  0.9);

	bool Insert(Individual &I);		//seed a starting solution, should be called after Configure()
	void Run(GENERIC_POINTER);		//performs optimization, calls a user-defined function each cycle

private:
	GENERIC_POINTER	score_func;		//pointer to a score-function
	REALNUM_TYPE rtPMutate;			//Mutation probability of one gene
	UNSIGNED_2B_TYPE	nModel;		//size of the model

	//the optimization criteria
	bool ifContMode;				//flag to continue at the same T after a good mutation
	UNSIGNED_4B_TYPE nIter;			//number of mutation trials before cooling
	REALNUM_TYPE T_0, T_END, K, T_CONV;	//T(i+1) = K*Ti; defs: 10^3, 10^-6 and 0.9 resp; T_CONV - convergence span

	UNSIGNED_2B_TYPE nBest;			//to store n best solutions, default is 1
	REALNUM_TYPE ideal_score;		//should be from 0 to 1

public:
	apvector<class Individual> models;
};


//Artificial Ant Colony Optimization (ACO)

class ACO
{
public:
	ACO();
	virtual ~ACO();

	void Configure(
	//configuration parameters 
		UNSIGNED_2B_TYPE 	= DEFAULT_POOL,				//nFeatures
		UNSIGNED_2B_TYPE 	= DEFAULT_SIZE, 			//nAnts
		UNSIGNED_2B_TYPE 	= DEFAULT_MIN_IND_SIZE,
		UNSIGNED_2B_TYPE 	= DEFAULT_MAX_IND_SIZE,		//Min and Max Solution size 
		UNSIGNED_2B_TYPE 	= DEFAULT_MAXGEN,			//nGenerations
		REALNUM_TYPE 		= DEFAULT_MIN_DIVERSION,	//Convergence criteria
		UNSIGNED_2B_TYPE	= DEFAULT_STABILITY_RANGE,	//nGenerations without change
		REALNUM_TYPE 		= DEFAULT_IDEAL_FITNESS,	//ideal value of the score-function
		GENERIC_POINTER 	= NULL,		//score function
		REALNUM_TYPE 		= 0.25,		//pheromone volatility rate (evaporating fraction)
		REALNUM_TYPE 		= 1.0,		//pheromone power-factor (importance of pheromone marks)
		REALNUM_TYPE 		= 1.0,		//feature importance power-factor (affects pheromone increament)		
		REALNUM_TYPE 		= 0.0, 		//min. pheromone level
		REALNUM_TYPE 		= 1.0, 		//max. pheromone level
		UNSIGNED_1B_TYPE	= 0,		//Workmode for pheromone schelduling(bitwise): 0 - regular update; 1 - only by iteration best, 2 - only by global best, 
		//4 - additional update by global best every iter, 8 - ASRank (pheromone increments are by ranking scheme), 16 - final run with Tmax after convergence
		UNSIGNED_2B_TYPE 	= 1);		//number of best models to remember

	bool Insert(Individual &);		//seed a starting solution, should be called after Configure()
	void Run(GENERIC_POINTER);		//performs optimization, calls a user-defined function each cycle

private:

	void update(Individual &);		//updates pheromon level from Individual solution I
	void resample(Individual &);	//resample descriptors based on pheromarks

	GENERIC_POINTER		sfunc;		//pointer to a score-function
	UNSIGNED_2B_TYPE	mnN, mxN;	//solution size boundaries
	UNSIGNED_2B_TYPE	nCyc;		// max.iterations

	//the optimization criteria
	REALNUM_TYPE p, tmin, tmax;		//pheromone volatiliy and min. and max. levels
	REALNUM_TYPE alpha, beta;		//power factors controlling pheromone importance and increments
	REALNUM_TYPE ph_sum;			//cumulative level of pheromon, used in resample()
	UNSIGNED_1B_TYPE phUpdate;		//pheromone schedule mode: 0 - regular update; 1 - only by iteration best, 2 - only by global best, 
									//4 - additional update by global best every iter, 8 - ASRank (pheromone increments are by ranking scheme), 16 - final run with Tmax after convergence	

	REALNUM_TYPE ideal_score;		//should be from 0 to 1; not used if INVALID
	REALNUM_TYPE converg_thr;		//minimum difference of scoring f() as a sign of convergence
	UNSIGNED_2B_TYPE stbl_range;	//max.number of iterations without improvement;
	
	apvector<class Individual> anthill;		//solutions explored by ants
	apvector<SIGNED_4B_TYPE> ipheromarks;	//indices of descriptors sorted by their pheromon levels
public:
	apvector<REALNUM_TYPE> pheromarks;
	apvector<class Individual> models; //to store n best solutions, default is 1
};



//Particle Swarm Optimization (PSO), includes binary mode as well
//new v = wv + c1*(x0 - x) + x2(localbest - x)
class particle
{
public:
	particle();
	virtual ~particle();
	void toss(UNSIGNED_2B_TYPE);
	void navigate(apvector<REALNUM_TYPE> &, REALNUM_TYPE, REALNUM_TYPE, REALNUM_TYPE, REALNUM_TYPE);
	void fly();

	REALNUM_TYPE q0, q;		//q - current and q0 - best-so-far quality values
	apvector<REALNUM_TYPE> x0, x, v;	//x - current position, x0 - best-so-far, v - velocity
};

/*
Configuration defaults for w, c1, c2 are from Pedersen, M.E.H. 2010 Thesis. 
However, other settings are possible, i.e. w=-0.6 and negative c1 to stimulate exploration (Pedersen, 2008)

TODO: 
Probability based binary PSO (BPSO), using sigmoid curve and two sets of positions and velocities, for each bit (0 and 1)
*/

class PSO
{
public:
	PSO();
	virtual ~PSO();

	void Configure(		
		UNSIGNED_2B_TYPE 	= DEFAULT_POOL,				//#Features
		UNSIGNED_2B_TYPE 	= DEFAULT_SIZE, 			//#particles
		UNSIGNED_2B_TYPE 	= DEFAULT_MAXGEN,			//#cycles
		REALNUM_TYPE 		= DEFAULT_MIN_DIVERSION,	//Convergence threshold
		UNSIGNED_2B_TYPE	= DEFAULT_STABILITY_RANGE,	//max #cycles without change
		GENERIC_POINTER 	= NULL,						//score function
		REALNUM_TYPE 		= 0.73,						//w, inertia
		REALNUM_TYPE 		= 1.5,						//c1, weight for  cognitive term
		REALNUM_TYPE 		= 1.5,						//c2, weight for social term
		REALNUM_TYPE 		= 3.0, 						//Vmax. speed limit
		UNSIGNED_2B_TYPE 	= 0);						//k, neighborhood for socializing, 0 means global

	//Insert() seeds a starting position, should be called after Configure()
	bool Insert(Individual &);		
	bool Insert(apvector<REALNUM_TYPE> &);

	//performs particle swarm optimization, calls a user-defined function each cycle
	void Run(GENERIC_POINTER);		

private:
	UNSIGNED_2B_TYPE ndims;
	REALNUM_TYPE vmax;
	REALNUM_TYPE w, c1, c2;		//weights for inertia, and for cognitive and socializing terms
	SIGNED_4B_TYPE allbest;		//index of the global best particle
	UNSIGNED_1B_TYPE k;			//number of neighbours to use, if k==0 then allbest is used instead
	 
	GENERIC_POINTER		sfunc;		//pointer to a score-function	

	UNSIGNED_2B_TYPE	nCyc;		// max.iterations
	REALNUM_TYPE converg_thr;		//minimum difference of scoring f() as a sign of convergence
	UNSIGNED_2B_TYPE stbl_range;	//max.number of iterations without improvement;
public:
	apvector<particle> swarm;
};

//service function to export solutions
void ClearDuplIndividuals(apvector<Individual> &, REALNUM_TYPE, REALNUM_TYPE, bool = true, bool = false);
void SortIndividuals(apvector<Individual> &);
void CopyGenes2Set(apvector<class Gene> &, set &);

#endif // !defined(FSEL_ALG_CLASS)

