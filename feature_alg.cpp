//feature_alg.cpp
//Implementation of feature-selection algorithms
/*
//------------------------------------------
//  Description of the GA algorithm
//	1)  Form the initial population of individuals using pool of genes, 
//		this process is called initialization, making up of individuals is random
//	2)	then using various algorithms two parents are selected and children are produced.
//		next generation is created by the process of randomly picking pairs of "parents" 
//		among these elite individuals and producing children individuals by the process 
//		of crossover and mutation. Only children form a new generation, but if crossover 
//		rate is not 100% then some children can be identical to their parents
//		
//	3)  Fitness function MUST BE coded in a way that higher values reflect better fitted individuals
//	
//	For more details about how the GA algorithm works please look fallowing reference:
//		Reference:	Ron Wehrens, Lutgarde M.C.Buydens, Evolutionary optimisation: a tutorial,
//				Trends in Analytical Chemistry, volume 17, N4, 1998, p193-203

//------------------------------------------
//	Simulated Annealing (SA) Algorithm 2008-2009
//  NB: 
//	1) Current solution is compared with best or current 
//	2) solution of equal quality with the best is always selected due to probability function: exp(-dq/Tc)
//
//------------------------------------------
//	Ant Colony Optimization (ACO) 2010
//	Population of N ants builds models simultaneously. Selected descriptors then marked by pheromone level to relfect resulting quality of the models
//	
//////////////////////////////////////////////////////////////////////
//	TODO:
//	1) Maybe for GA specify mutation rate on per gene basis also? 
//	2) Parallelism for GA?
*/

#include "feature_alg.h"

//-------   memory leaks catcher for the current source-file  --------
#ifdef ADV_LEAK_CATCHER
#ifdef _DEBUG 
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
#endif
//--------------------------------------------------------------------


UNSIGNED_4B_TYPE RandomN(UNSIGNED_4B_TYPE N)
{//local random integer number generator in the range [0.. N-1] N should be less than RAND_MAX
	return GetRandomNumber(N);
};

REALNUM_TYPE RandomRN(REALNUM_TYPE S, REALNUM_TYPE E)
{
	REALNUM_TYPE R  = rand();

	R /= RAND_MAX;
	R *= E - S;
	R += S;

	return R;
};

//////////////////////////////////////////////////////////////////////
//Gene Class
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Gene::Gene()
{
	GD = N  = 0;
	A = B = GR = 0;
	gen_pool_limited = GENE_INTEGER;
}

Gene::Gene(UNSIGNED_4B_TYPE NumGen)
{	
	N = GD = 0;
	A = B = GR = 0;
	gen_pool_limited = GENE_INTEGER;

	if (NumGen)
	{
		N  = NumGen;
		GD = RandomN(N);				
	};
}

Gene::Gene(REALNUM_TYPE START, REALNUM_TYPE END)
{
	N = GD = 0;
	A = B = GR = 0;
	gen_pool_limited = GENE_INTEGER;
		
	if (START < END)
	{
		A = START;
		B = END;

		GR = RandomRN(A, B);
		gen_pool_limited = GENE_REAL;
	};	
}

Gene::Gene(const Gene &G)
{
	GD = G.GD;
	N = G.N;
	
	GR = G.GR;
	A = G.A;
	B = G.B;	

	gen_pool_limited = G.gen_pool_limited;
}

Gene::~Gene()
{
	N = GD = 0;
	A = B = GR = 0;
	gen_pool_limited = GENE_INTEGER;
}

//////////////////////////////////////////////////////////////////////
// End of Construction/Destruction
//////////////////////////////////////////////////////////////////////

//comparison operators for Gene class
bool operator == (const Gene &A, const Gene &B)
{
	if (A.gen_pool_limited != B.gen_pool_limited)
		return false;

	if (A.gen_pool_limited)
	{
		if (A.N != B.N)
			return false;

		return (A.GD == B.GD);
	}
	else
	{
		if ((fabs(A.A - B.A) > SMALL_NUMBER) || (fabs(A.B - B.B) > SMALL_NUMBER))
			return false;

		return (fabs(A.GR -  B.GR) < SMALL_NUMBER);	
	};

	return false;
}

bool operator != (const Gene &A, const Gene &B)
{
		if (A.gen_pool_limited != B.gen_pool_limited)
		return true;

	if (A.gen_pool_limited)
	{
		if (A.N != B.N)
			return true;

		return (A.GD != B.GD);
	}
	else
	{
		if ((fabs(A.A - B.A) > SMALL_NUMBER) || (fabs(A.B - B.B) > SMALL_NUMBER))
			return false;

		return (fabs(A.GR -  B.GR) > SMALL_NUMBER);
	};

	return true;
}

bool operator < (const Gene &A, const Gene &B)
{
		if (A.gen_pool_limited != B.gen_pool_limited)
		return false;

	if (A.gen_pool_limited)
	{
		if (A.N != B.N)
			return false;

		return (A.GD < B.GD);
	}
	else
	{
		if ((fabs(A.A - B.A) > SMALL_NUMBER) || (fabs(A.B - B.B) > SMALL_NUMBER))
			return false;

		return (A.GR < B.GR);	
	};

	return false;
}

bool operator > (const Gene &A, const Gene &B)
{
		if (A.gen_pool_limited != B.gen_pool_limited)
		return false;

	if (A.gen_pool_limited)
	{
		if (A.N != B.N)
			return false;

		return (A.GD > B.GD);
	}
	else
	{
		if ((fabs(A.A - B.A) > SMALL_NUMBER) || (fabs(A.B - B.B) > SMALL_NUMBER))
			return false;

		return (A.GR > B.GR);	
	};

	return false;
}

bool operator >= (const Gene &A, const Gene &B)
{
		if (A.gen_pool_limited != B.gen_pool_limited)
		return false;

	if (A.gen_pool_limited)
	{
		if (A.N != B.N)
			return false;

		return (A.GD >= B.GD);
	}
	else
	{
		if ((A.A != B.A) || (A.B != B.B))
			return false;

		return (A.GR >= B.GR);	
	};

	return false;
}

bool operator <= (const Gene &A, const Gene &B)
{
		if (A.gen_pool_limited != B.gen_pool_limited)
		return false;

	if (A.gen_pool_limited)
	{
		if (A.N != B.N)
			return false;

		return (A.GD <= B.GD);
	}
	else
	{
		if ((A.A != B.A) || (A.B != B.B))
			return false;

		return (A.GR <= B.GR);	
	};

	return false;
}

void Gene::MutateGene()
//performs mutation on the gene. 
//Note! May be it makes sense to introduce minimal 
//		distance between mutated and previous genes
//		in case of unlimited gene values (REAL) ? 
//		Right now such a mutation can lead to almost the same invidiual
{
	UNSIGNED_4B_TYPE OLD_GD = GD;
	REALNUM_TYPE OLD_GR = GR;
	if (gen_pool_limited)
	{
		do {
			GD = RandomN(N);
		} while (GD == OLD_GD);

	}
	else
	{
		do {
			GR =  RandomRN(A, B);
		} while (GR == OLD_GR);
	}
}

//////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////
//Individual Class
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Individual::Individual()
{
	quality = INVALID;
	SameGeneShouldNotOccur	= false;

	body.resize(0);
}

Individual::Individual(UNSIGNED_4B_TYPE N)
{
	quality = INVALID;
	SameGeneShouldNotOccur	= false;
	body.resize(N);
}

Individual::Individual(UNSIGNED_4B_TYPE N, UNSIGNED_4B_TYPE GENEPOOL, bool SameGeneF = false)
//randomizes genes automatically!
{
	set Genes;	
	
	quality = INVALID;
	SameGeneShouldNotOccur	= SameGeneF;
	body.resize(N);

	if (N > GENEPOOL)		//if gene sequence is longer than size of gene pool then put the flag is OFF
		SameGeneShouldNotOccur = false;

	for (UNSIGNED_4B_TYPE i=0; i<N; i++)
	{
		Gene Temp(GENEPOOL);
		body[i] = Temp;

		if (SameGeneShouldNotOccur)
		{
			while (Genes.IsInSet(body[i].GD))
				body[i].MutateGene();

			Genes.PutInSet(body[i].GD);
		};
	};	
}

Individual::Individual(UNSIGNED_4B_TYPE N, REALNUM_TYPE A, REALNUM_TYPE B, bool SameGeneF = false)
{
	UNSIGNED_4B_TYPE i, j;
	
	quality = INVALID;
	SameGeneShouldNotOccur	= SameGeneF;
	body.resize(N);
	
	for (i=0; i<N; i++)
	{
		Gene Temp(A, B);
		body[i] = Temp;

		if (SameGeneShouldNotOccur)		
			do 
			{
				for (j=0; j<i; j++)
					if (body[j].GR == body[i].GR)					
						break;
			
				if (j<i)
					body[i].MutateGene();
			} while (j<i);		
	};	
}

Individual::Individual(const Individual &I)
:body(I.body)
{
	quality = I.quality;
	SameGeneShouldNotOccur	= I.SameGeneShouldNotOccur;
}

Individual::~Individual()
{
	body.resize(0);
	quality = INVALID;
	SameGeneShouldNotOccur	= false;
}

//////////////////////////////////////////////////////////////////////
// End of Construction/Destruction
//////////////////////////////////////////////////////////////////////

void	Individual::MutateIndividual()
//member function that performs mutation
{
	UNSIGNED_4B_TYPE bl = (UNSIGNED_4B_TYPE)body.length();
	UNSIGNED_4B_TYPE s, uG = RandomN(bl);

	body[uG].MutateGene();
	quality = INVALID;
	
	if (!SameGeneShouldNotOccur) return;

	do
	{//check for duplicate gene-values, but using "set" class is not reasonable!
		for (s = 0; s < bl; s++)
		{
			if (s == uG) continue;
			
			if ( (body[uG].gen_pool_limited && (body[s].GD == body[uG].GD)) ||
				((!body[uG].gen_pool_limited) && (body[s].GR == body[uG].GR)) )
			{
				body[uG].MutateGene();
				break;
			}
		}
	} while (s < bl);
}

void	Individual::EvaluateIndividual(GENERIC_POINTER DFunct)
//calls  D-Function to evaluate individual
{
	REALNUM_TYPE (*DFunction)(apvector<class Gene> &);
	
	if (DFunct == NULL)
		return;

	DFunction = (REALNUM_TYPE (*)(apvector<class Gene> &)) DFunct;
	quality = DFunction(body);
}

void	Individual::SetGenesTo(UNSIGNED_4B_TYPE DefD)
{
	if (SameGeneShouldNotOccur) return;
	for (SIGNED_4B_TYPE s = 0; s < body.length(); s++)	body[s].GD = DefD;
}

void	Individual::SetRGenesTo(REALNUM_TYPE DefR)
{
	if (SameGeneShouldNotOccur) return;
	for (SIGNED_4B_TYPE s = 0; s < body.length(); s++)	
	{
		if (body[s].A > DefR) continue;
		if (body[s].B < DefR) continue;
		body[s].GR = DefR;
	}
}
bool Individual::CompareExact(Individual &In)
//for integer genes only, checks if soluations are the same
{
	if (body.length() != In.body.length()) return false;
	set xx1, xx2;
	CopyGenes2Set(In.body, xx1);
	CopyGenes2Set(body, xx2);
	xx1 -= xx2;
	return xx1.IsEmpty();
}

REALNUM_TYPE Individual::Compare(Individual &In)
{//calculates similarity to In (in terms of probablity of being unrelated)
	REALNUM_TYPE p = 1.0, fx, divz;
	UNSIGNED_4B_TYPE k, L1 = body.length(), L2 = In.body.length(), M;
	bool md = (body[ZERO].gen_pool_limited == GENE_REAL);

	M = max(L1, L2);
	if (md)
		divz = body[ZERO].B - body[ZERO].A;
	else
		divz = body[ZERO].N;

	for (k = ZERO; k< M; k++)
	{
		if (md)
		{
			fx = ZERO;
			if (k < L1) fx = body[k].GR;
			if (k < L2) fx -= In.body[k].GR;
			fx = fabs(fx);
		}
		else
		{
			fx = 1;
			if ((k < L1) && (k < L2)) 
				if (body[k].GD == In.body[k].GD) 
					fx = 0;
		}

		fx /= divz;		
		p *= 1 - fx;
	}
	return (1 - p);
}


//***********************************************************************************************
//service functions for operations with class Individual

void SortIndividuals(apvector<Individual> &P)
//precondition:		fitness factor (quality-field had to be evaluated for each individual)
//descriptioni:		sorts individuals by the fitness factor
//postcondition:	individuals are sorted, better fitted in the beginning of the array.

//Quicksort
//description:		sorts an array P into descending numerical order of its quality-field
//					using the Quicksort algorithm. 
//					P is replaced on output by its sorted rearrangement.
{
	SIGNED_4B_TYPE N = P.length();
	SIGNED_4B_TYPE *A	= GRAB_MEM_BLOCKS(SIGNED_4B_TYPE, N), ic = ZERO;
	REALNUM_TYPE *F		= GRAB_MEM_BLOCKS(REALNUM_TYPE, N);
	
	QSortScore = F;
	for (; ic < N; ic++)
	{
		A[ic] = ic;	
		F[ic] = P[ic].quality; 
	}
	qsort(A, (size_t)N, sizeof(SIGNED_4B_TYPE), QSortCompareLess);
	
	DROP_MEM_BLOCKS(F);
	QSortScore = NULL;

	apvector<Individual> NewP(N);
	for (ic = ZERO; ic < N; ic++)	NewP[ic] = P[A[ic]];
	DROP_MEM_BLOCKS(A);
	P = NewP;
}

void CopyGenes2Set(apvector<class Gene> &I, set &GeneSet)
{//works with descrete genes only!
	UNSIGNED_4B_TYPE i, Nm = I.length();
	GeneSet.Dump();
	for (i = 0; i < Nm; i++)	GeneSet.PutInSet(I[i].GD);	
}

void ClearDuplIndividuals(apvector<Individual> &P, REALNUM_TYPE MinQualDiff, REALNUM_TYPE MinPDiff, bool ifResize, bool ExactMatch)
//precondition:		individuals have to be sorted by their fitness.
//description:		deletes individuals that have identical gene sequences (order matters too)
{
	bool Duplicates = false, Exo = (ExactMatch || P[0].SameGeneShouldNotOccur);	//Dec 24 2010, to force exact match if genes are not on/off but coded by id (so different orders could make false mismatches)
	UNSIGNED_4B_TYPE i, j, L = P.length();

	i = 0;
	while (i < L-1)	
	{
		j = i+1;
		while (j < L)
		{
			if ( (P[i].quality - P[j].quality) > MinQualDiff ) break;

			bool isMatch = false;
			if (Exo)
				isMatch = P[i].CompareExact(P[j]);			
			else
				isMatch = (P[i].Compare(P[j]) <= MinPDiff);

			if (isMatch)
			//sequences are identical. Now let's mark individual j as deletable
			{
				P[j].quality = INVALID;
				Duplicates = true;
			};
			j++;
		};//while j

		i = j;
	};//while i

	if (!Duplicates) return;

	//now purge those individuals that have INVALID in their quality fields.
	i = j = 0;
	while (i<L)
	{
		if (P[i].quality != INVALID)
		{
			if (i!= j)	P[j] = P[i];
			j++;
		};
		i++;
	};

	if (ifResize) //Crowding mode
		P.resize(j);
	else
		for (; j<L; j++)
		{//Sharing mode
			i = P[j].body.length() >> 2; //mutate 25%
			while (i--) P[j].MutateIndividual();
		}
}

//***********************************************************************************************
//***********************************************************************************************
//***********************************************************************************************

//GenAlgorithm Class

GenAlgorithm::GenAlgorithm()
{
	Configure();
}

GenAlgorithm::~GenAlgorithm()
{
	population.resize(0);
	next_generation.resize(0);
}

void GenAlgorithm::PickParent(UNSIGNED_4B_TYPE &ID)
//precondition: quality-field of individuals had to be evaluated and contain "fittness" parameter
//				for RANK_BASED mode individuals also had to be sorted in the population by the fitness value.
//				ID contains either another parent which was picked already, or should be set to value bigger than population.length()
//description:  function overwrites ID with ID of picked individual
//
{
	UNSIGNED_4B_TYPE j, IDNEW, Elite, L = population.length();
	REALNUM_TYPE El = ElitePart*L;
	
	switch (PSelMode)
	{
		case RANK_BASED:
			Elite = (UNSIGNED_4B_TYPE)El;
			do  IDNEW = RandomN(Elite); while (IDNEW == ID);
			ID = IDNEW;
			break;

		case ROULETTE:
			do 			
				do 	Elite = RandomN(L); while (Elite == ID);
			while (population[Elite].quality < RandomRN(population[L-1].quality, population[0].quality));
			ID = Elite;
			break;

		case TOURNAMENT:
			do 	Elite = RandomN(L); while (Elite == ID);

			for (j=1; j<GroupSize; j++)		
			{
				do { IDNEW = RandomN(L);} while (IDNEW == ID);
				
				if (population[IDNEW].quality > population[Elite].quality)
					Elite = IDNEW;
			};

			ID = Elite;			
			break;

		default:
			//set the error
			ID = population.length()+1;
			break;
	};//end of switch (PSelMode)
}

void GenAlgorithm::CrossOver(UNSIGNED_4B_TYPE A, UNSIGNED_4B_TYPE B, UNSIGNED_4B_TYPE &POS)
//description:		peforms crossover of two parents
//					automatically adds 2 generated children to next_generation
//					also performs mutation of those children

//NB:				crossover can introduce duplicate genes (if 'SameGeneShouldNotOccur' is true,
//					'SameGeneShouldNotOccur' is enforced only for the mutated genes)
//					these duplicates must be handled in the fitness f()
{
	UNSIGNED_4B_TYPE i, j, k, la = population[A].body.length(), lb = population[B].body.length();

	if (la > lb) { CrossOver(B, A, POS); return; };
	//now we can assume that la < lb
	Individual N1(population[A]), N2(population[B]);

	UNSIGNED_4B_TYPE crs = CrossOverMode, cpz = RandomN(lb - la + 1);
	if (crs & CRS_ALIGN_HEADS)
	{//do not randomize alignment
		cpz = 0;
		crs -= CRS_ALIGN_HEADS;
	};
	//if (MinIndSize < MaxIndSize) //size-change during crossover is probably not needed

	if (crossover_chance > RandomRN(0, 1))
	{
			N2.quality = N1.quality = INVALID;
			switch (crs)
			{
				case ONE_POINT:
					i = RandomN(la - 1);					
					for (j=0; j<=i; j++)
					{
						N1.body[j] = population[B].body[cpz + j];
						N2.body[cpz + j] = population[A].body[j];
					};			
				break;

				case TWO_POINT: //always two-points crossover
					i = RandomN(la - 2) + 1;
					k = RandomN(la - i - 1) + i; //if k=i, then middle area = 1
					for (j=i; j<=k; j++)
					{
						N1.body[j] = population[B].body[cpz + j];
						N2.body[cpz + j] = population[A].body[j];
					};
				break;

				case TWO_OR_ONE_POINT:
					i = RandomN(la -1);
					k = RandomN(la -1);

					if (k>i)
					//two-point crossover
						for (j=i+1; j<=k; j++)
						{
							N1.body[j] = population[B].body[cpz + j];
							N2.body[cpz + j] = population[A].body[j];
						}
					else
					//one-point crossover only
						for (j=0; j<=i; j++)
						{
							N1.body[j] = population[B].body[cpz + j];
							N2.body[cpz + j] = population[A].body[j];
						};			
				break;

				case UNIFORM:
					for(j=0; j<la; j++)
					if (RandomN(2))
					{//50% chance to switch two genes
						N1.body[j] = population[B].body[cpz + j];
						N2.body[cpz + j] = population[A].body[j];
					};
				break;
			}//switch
	}	//if (crossover_chance > RandomRN(0, 1))

	//perform mutations
	if (mutation_rate > RandomRN(0, 1))
		N1.MutateIndividual();

	if (mutation_rate > RandomRN(0, 1))
		N2.MutateIndividual();
	
	//saving children to the next_generation array	
	next_generation[POS++] = N1;
	next_generation[POS++] = N2;
}

void GenAlgorithm::ProduceNewGeneration()
//precondition:		individual in population had to be sorted by their fitness
//postcondition:	next_generation contain new generation of individuals
//					next_generation fitness is evaluated, and individuals are
//					sorted by this fitting factor.
//					NOTE: next generation can be less the p_size if some duplicate 
//					individuals were detected and purged
{
	UNSIGNED_4B_TYPE i, j, l, n, k;
	
	l = population.length();
	k = n = 0;
	
	next_generation.resize(0);
	next_generation.resize(p_size+1);

	if (Elitism)
	{//copy best parents to the next generation automatically

		k = n = (UNSIGNED_4B_TYPE)(ElitePart*l);
		
		for (i = 0; i< n; i++)		
			next_generation[i] = population[i];		
	};	

	do 
	{
		if ( next_generation.length() < SIGNED_4B_TYPE(p_size+1) )
			next_generation.resize(p_size+1);

		while (n < p_size)
		{
			j = i = l;
			PickParent(j);
			PickParent(i);
			CrossOver(i, j, n);
		};

		if (n == p_size)
		//copy a best individual from the previous population if 1 place is still available
			next_generation[p_size] = population[k];
		
		EvaluateIndividuals(next_generation);
		SortIndividuals(next_generation);
		ClearDuplIndividuals(next_generation, MinDiversion, MinDiversion);
		n = next_generation.length();

	}while ((p_size +1 - n) > (p_size >> 2));		//if new generation is less then 75% of the previous one then continue multiplying	
}

void GenAlgorithm::EvaluateIndividuals(apvector<Individual> &P)
{
	if (DFUNCTION == NULL)
		return;

	for (SIGNED_4B_TYPE j = 0; j<P.length(); j++)
		if (P[j].quality == INVALID)
			P[j].EvaluateIndividual(DFUNCTION);
};

bool GenAlgorithm::InsertIndividual(Individual I)
//precondition:		should be used after call of Configure()
//description		inserts an individual (if it's valid) into the  
//					population instead of one of existing individuals
//					which are randomly generated by Configure()
{
	UNSIGNED_4B_TYPE bL = I.body.length();
	if ( (MinIndSize > bL) || (MaxIndSize < bL) )
		return false;

	//if given valid is valid then let's insert it in a random place
	population[RandomN(population.length())] = I;

	return true;
}

void GenAlgorithm::Configure
(
 //first parameters that describe individual go
 bool SAMEGENEALLOWED,
 bool GENEPOOLLIMITED,
 bool VARIATIONOFSIZE_MODE,
 UNSIGNED_4B_TYPE GENES_POOL,
 REALNUM_TYPE GENES_RANGE_S,
 REALNUM_TYPE GENES_RANGE_E,
 UNSIGNED_4B_TYPE INDIVIDUAL_SIZE,
 UNSIGNED_4B_TYPE MIN_INDIVIDUAL_SIZE,
 UNSIGNED_4B_TYPE MAX_INDIVIDUAL_SIZE,

 //parameters of the algorithm
 UNSIGNED_4B_TYPE POPULATION_SIZE, 
 UNSIGNED_4B_TYPE MAX_N_GENERATIONS, 
 REALNUM_TYPE MUTATION_RATE,
 REALNUM_TYPE CROSSOVER_CHANCE,
 
 REALNUM_TYPE MIN_DIVERSION,
 UNSIGNED_4B_TYPE STABILITY_RANGE,
 REALNUM_TYPE IDEAL_FIT,
 GENERIC_POINTER DFUNCT,

 //modes
 UNSIGNED_4B_TYPE CROSSOVER_MODE,
 UNSIGNED_4B_TYPE PARENTSELECTION_MODE,  
 bool ELITISM_MODE,

 //some additional parameters needed for some modes
 UNSIGNED_4B_TYPE GROUP_SIZE,
 REALNUM_TYPE ELITEPART_SIZE)
//
//description :		1) Sets work-parameters
//					2) Seeds initial population of individuals
//					3) Randomizes starting point of the pseudo-
//					   random numbers generator
//
//postcondition:	parameters of the genetic algorithm are all 
//					configured and population of individuals is ready to use
//
{
	UNSIGNED_4B_TYPE i, ind_size;

	if (VARIATIONOFSIZE_MODE)
	{
		MaxIndSize = MAX_INDIVIDUAL_SIZE;
		MinIndSize = MIN_INDIVIDUAL_SIZE;	
	}
	else
		MaxIndSize = MinIndSize = INDIVIDUAL_SIZE;

	p_size = POPULATION_SIZE;	

	population.resize(p_size);

	//initialization of the population
	for (i = 0; i<p_size; i++)
	{
		if (VARIATIONOFSIZE_MODE)
			//if size can be changed then set random length
			ind_size = RandomN(MaxIndSize - MinIndSize + 1) + MinIndSize;
		else
			ind_size = INDIVIDUAL_SIZE;

		if (GENEPOOLLIMITED)
			population[i] = Individual(ind_size, GENES_POOL, !SAMEGENEALLOWED);
		else
			population[i] = Individual(ind_size, GENES_RANGE_S, GENES_RANGE_E, !SAMEGENEALLOWED);
	};

	
	MaxNGener			= MAX_N_GENERATIONS;	

	mutation_rate		= MUTATION_RATE;
	crossover_chance	= CROSSOVER_CHANCE;
	
	DFUNCTION			= DFUNCT;	
	
	MinDiversion		= MIN_DIVERSION;
	StabilityRange		= STABILITY_RANGE;
	IdealFitness		= IDEAL_FIT;
		

	CrossOverMode		= CROSSOVER_MODE;
	PSelMode			= PARENTSELECTION_MODE;

	IndSizeConstant		= VARIATIONOFSIZE_MODE;	
	Elitism				= ELITISM_MODE;	
	
	GroupSize			= GROUP_SIZE;
	ElitePart			= ELITEPART_SIZE;
}

bool GenAlgorithm::CheckSettings()
{
	if ( (p_size == 0) || (SIGNED_4B_TYPE(p_size) != population.length()) )	return false;
	if ( (MaxIndSize < MinIndSize) || (MinIndSize < 1) ) return false;
	if (DFUNCTION	== NULL)	return false;
	if (GroupSize < 2)	return false;
	if ( (mutation_rate<= 0) || (crossover_chance<= 0) ) return false;
	if (Elitism && (ElitePart <= 0)) return false;

	for (UNSIGNED_4B_TYPE sc = 0; sc < p_size; sc++)
	for (SIGNED_4B_TYPE sb = 0; sb < population[sc].body.length(); sb++)
	if (population[sc].SameGeneShouldNotOccur)
	if (population[sc].body[sb].gen_pool_limited)
	if (population[sc].body.length() > SIGNED_4B_TYPE(population[sc].body[sb].N) )
		return false; //#genes should be <= gene pool size!
	return true;
}

bool GenAlgorithm::PrepareOptimization()
//description:		1) checks settings and returns FALSE if some 
//					error was detected
//					2) evaluates fitness factor for each individual 
//					of the initial population 
//					3) Purges duplicates if any
//postcondition:	counter of generations currGener is set to 0
//					Range varialbe is set to Stability Range
{
	currGener = ZERO;
	Range	  = StabilityRange;

	if (!CheckSettings())
		return false;
	
	EvaluateIndividuals(population);
	SortIndividuals(population);
	ClearDuplIndividuals(population, MinDiversion, MinDiversion);

	return true;
}

bool GenAlgorithm::PerformOneCycle(GENERIC_POINTER OUTPUT_FUNCT = NULL)
//description:	produces one geenration and evaluate it
//				returns true if final conditions are met and optimization is done
//				PrepareOptimization() had to be run at least 
//				once before any call of PerformOneCycle
{
	bool STOP = false, Output = false;
	char buffer[200];
	STRING_TYPE S;

	void (*OutputF)(STRING_TYPE) = (void (*)(STRING_TYPE)) OUTPUT_FUNCT;
	if (OUTPUT_FUNCT != NULL) Output = true;
		
	ProduceNewGeneration();

	//fabs() is needed in case if next generations has worse quality!
	Diversion = fabs(next_generation[0].quality - population[0].quality);
	population = next_generation;
	currGener++;	
	
	if (Output)
	{
		sprintf(buffer, "Generation N%d produced.", currGener);
		OutputF(buffer);	
		sprintf(buffer, "Best Fitness Factor: %12.9f", population[0].quality);
		OutputF(buffer);
		sprintf(buffer, "Diversion b\\w 2 populations: %12.9f", Diversion);
		OutputF(buffer);
		OutputF("");
	};
	
	//--- check conditions for the end of optimization ---//
	if (Diversion <= MinDiversion)	Range--;	else	Range = StabilityRange;

	if (Range == ZERO)
	{
		STOP = true;
		if (Output)
		{
			sprintf(buffer, "***   During last %d generations diversion between them was less than minimal diversion (%12.9f)", StabilityRange, MinDiversion);
			OutputF(buffer);
		}
	};

	if ((population[0].quality >= IdealFitness) && (IdealFitness != INVALID))
	{
		STOP = true;
		if (Output)
			OutputF("***   Ideal fitness factor is reached.");
	};

	if (currGener == MaxNGener)
	{
		if (Output)
			OutputF("***   Maximum number of generations is reached.");
		STOP = true;
	};
	
	if (STOP && Output)						// ----- output	
		OutputF("Final conditions met. Optimization is finished.");

	return (STOP);
}

void GenAlgorithm::PerformOptimization(GENERIC_POINTER CF = NULL, GENERIC_POINTER OF = NULL)
//precondition:		Configure() had to be called before
//					any calls of InsertIndividual() also should be 
//					done before	the call of this function if user 
//					wants to have some certain individuals to be 
//					present in the initial population.
//
//parameters-functions:
//					OF - output function :   void			OF (STRING_TYPE); 
//					CF - control function:   SIGNED_2B_TYPE CF (void);
//
//					OF is used to handle program messages and it is called occasionally
//					throughout a whole optimization procedure.
//
//					if CF is not NULL it will be called once per 
//					optimization cycle.
//					if return-value of CF has bit 0 set :			OPTIMIZATION is ABORTED
//					if return-value of CF has bit 1 set :			PAUSE flag ON
//					if return-value of CF has bit 1 is not set:		PAUSE flag OFF
//					Thus CF can be used to govern the optimization 
//					if it was opened in another thread and\or to 
//					output some statistics (plot a diagram and so on)
//		
//description:		Being a main tool for GenAlgorithm class this
//					function governs whole optimization process
{
	bool STOP = false, PAUSE = false;
	SIGNED_2B_TYPE CODE = ZERO;
		
	SIGNED_2B_TYPE (*FlagFunction)(void);
	
	if (CF != NULL)
		FlagFunction = (SIGNED_2B_TYPE (*)(void)) CF;
	else
		FlagFunction = NULL;

	if (!PrepareOptimization())
		return;

	do 
	{//main cycle

		if (CF != NULL)	CODE = FlagFunction();
		if (CODE & MODE_PAUSE)	PAUSE = true;	else	PAUSE = false;		
		if (CODE & MODE_ABORT)	STOP = true;		
		if (!PAUSE)		STOP = STOP || PerformOneCycle(OF);

	}while (!STOP);	
}

//Encapsulation of currGener field
UNSIGNED_4B_TYPE GenAlgorithm::GetCurrGener()
{
	return (currGener);
}

void GenAlgorithm::SetCurrGener(UNSIGNED_4B_TYPE C)
{
	currGener = C;
}



//***********************************************************************************************
//***********************************************************************************************
//***********************************************************************************************
//class for the Simulated Annealing algorithm

SimAnneal::SimAnneal()
{
	Configure();
}

SimAnneal::~SimAnneal()
{
	models.resize(0);
}

void SimAnneal::Configure(
	//parameters that describe individual solution
	 bool ifCanSameG, UNSIGNED_4B_TYPE nGPool,
	 REALNUM_TYPE rtGs, REALNUM_TYPE rtGe,
	 UNSIGNED_2B_TYPE nSize,
	//optimization parameters
	bool ifCM,
	UNSIGNED_2B_TYPE nB,
	GENERIC_POINTER gpFunc, REALNUM_TYPE rtIS,
	REALNUM_TYPE rtPM, UNSIGNED_4B_TYPE nI,
	REALNUM_TYPE rtT_0, REALNUM_TYPE rtT_END, 
	REALNUM_TYPE rtT_CONV, REALNUM_TYPE rtK)
{
	score_func = gpFunc;

	if ( (rtPM > 0) && (rtPM < 1.0) )
		rtPMutate = rtPM;

	nIter = max(1, nI);

	if (rtT_0 > rtT_END)
	{
		T_0 = pow(10, rtT_0);
		T_END = pow(10, rtT_END);
	}

	if ((rtK > 0.0) && (rtK < 1.0))	K = rtK;
	if (rtT_CONV < 0.0)	T_CONV	 = pow(10, rtT_CONV);

	nBest = nB;	//to store n best solutions, default is 1
	models.resize(nBest + 1);
	ideal_score = rtIS;
	ifContMode =  ifCM;

	for (UNSIGNED_4B_TYPE i = 0; i <= nBest; i++)
	{		
		
		Individual I_D(nSize, nGPool, !ifCanSameG),
		I_R(nSize, rtGs, rtGe, !ifCanSameG);

		if (nGPool)
			models[i] = I_D;
		else
			models[i] = I_R;
	}//for i
}

void SimAnneal::Run(GENERIC_POINTER OutputF)
{
	char bfr[256];
	bool ifOutput = false, ifDo = true;
	void (*outF)(STRING_TYPE) = (void (*)(STRING_TYPE)) OutputF;	

	if (OutputF != NULL) ifOutput = true;
	
	REALNUM_TYPE dQ, T_C = T_0, T_B = T_0;
	UNSIGNED_4B_TYPE nT, nMut, l = 0, kl;
	Individual indTmp;
	nT = UNSIGNED_4B_TYPE(rtPMutate * REALNUM_TYPE(models[nBest].body.length()) );

	for (; l <= nBest; l++)	models[l].EvaluateIndividual(score_func);	

	while (ifDo)
	{
		nMut = nIter;
		while (nMut)
		{
			indTmp = models[nBest];
			for (l = 0; l < nT; l++)	models[nBest].MutateIndividual();
			models[nBest].EvaluateIndividual(score_func);

			//1) we compare with the best of the stored 2) new candidate of equal quality is always selected over the old one!
			dQ = models[nBest].quality - models[0].quality;
			if ( RandomRN(0, 1) < exp(dQ/T_C) )
			{
				for (l = 0; l < nBest; l++)
				{
					if (models[l].quality > models[nBest].quality) continue;
					if (models[l].CompareExact(models[nBest])) break;
				
					if (models[l].quality < models[nBest].quality)
					{//if there are multiple equal-quality solutions, all will be checked above first, Oct.2010
						kl = nBest;
						while (--kl > l)	models[kl] = models[kl-1];
						models[l] = models[nBest];
						T_B = T_C;
						break;
					}
				}
				//finish the cycle if dQ > 0 and the workmode is permitting
				if (ifContMode || (dQ < 0))	nMut--;	else nMut = 0;
			}
			else
			{//put it back
				models[nBest] = indTmp;
				nMut--;
			}
		}//while (nMut)

		if (ifOutput)
		{
			sprintf(bfr, "current T: %f; best score: %f", T_C, models[0].quality);
			outF(bfr);
		}

		if (models[nBest].quality < models[0].quality) //seed best solution for next cycle
			models[nBest] = models[0];

		T_C *= K;
		if ( (ideal_score == models[0].quality) || (T_C < T_END) || (T_C/T_B < T_CONV) )
			ifDo = false;
	}// while (ifDo)
}

bool SimAnneal::Insert(Individual &I)
//precondition:		should be used after call of Configure()
//description		inserts an individual (if it is a valid one) to be used
//					as a starting point of optimization
{
	if (models[nBest].body.length() == I.body.length() )
	{
		models[nBest] = I;
		return true;
	}
	return false;
}

//***********************************************************************************************
//***********************************************************************************************
//***********************************************************************************************
ACO::ACO()
{
	Configure();
}

ACO::~ACO()
{
	models.resize(0);
	anthill.resize(0);
	pheromarks.resize(0);
	ipheromarks.resize(0);
}

void ACO::Configure(
	//configuration parameters 
		UNSIGNED_2B_TYPE nFeatures,
		UNSIGNED_2B_TYPE nAnts,
		UNSIGNED_2B_TYPE 	MinSize,
		UNSIGNED_2B_TYPE 	MaxSize,		//Min and Max Solution size 
		UNSIGNED_2B_TYPE 	nIterations,	//nGenerations
		REALNUM_TYPE 		minDiff,		//Convergence criteria
		UNSIGNED_2B_TYPE	stabilityRange,	//nGenerations without change
		REALNUM_TYPE 		idealFit,		//ideal value of the score-function
		GENERIC_POINTER 	dFunct,			//score function
		REALNUM_TYPE 		volatility,		//pheromone volatility rate (evaporating fraction)
		REALNUM_TYPE 		a,				//pheromone power-factor (importance of pheromone marks)
		REALNUM_TYPE 		b,				//feature importance power-factor (affects pheromone increament)		
		REALNUM_TYPE 		ph_min, 		//min. pheromone level
		REALNUM_TYPE 		ph_max, 		//max. pheromone level
		UNSIGNED_1B_TYPE	ph_mode,		//Workmode for pheromone schelduling(bitwise): 0 - regular update; 1 - only by iteration best, 2 - only by global best, 
		//4 - additional update by global best every iter, 8 - ASRank (pheromone increments are by ranking scheme), 16 - final run with Tmax after convergence
		UNSIGNED_2B_TYPE 	nBest)		//number of best models to remember
{
	alpha = a; 
	beta = b;
	
	sfunc	= dFunct;	
	
	converg_thr		= minDiff;
	stbl_range		= stabilityRange;
	ideal_score		= idealFit;
	
	anthill.resize(nAnts);
	pheromarks.resize(nFeatures);
	ipheromarks.resize(nFeatures);
	models.resize(nBest);
	
	mnN = MinSize;
	mxN = MaxSize;
	nCyc = nIterations;
	p = volatility; 
	
	tmin = ph_min;
	tmax = ph_max;
	
	//normalize by expected sampling power
	tmax *=	nAnts;
	tmax /= nFeatures;
	//tmax *= MinSize + MaxSize;	//no need to accout for solution size since pheromone contributions are size-scaled
	//tmax /= 2;
	
	
	UNSIGNED_2B_TYPE i = 0, N;
	for (i = 0; i < nFeatures; i++)	pheromarks[i] = tmax;
	for (i = 0; i < nBest; i++) models[i].quality = INVALID;

	for (i = 0; i < nAnts; i++)
	{
		N = (mnN == mxN) ? mnN : (mnN + RandomN(mxN - mnN + 1));
		anthill[i] = Individual(N, nFeatures, true);
	}

	for (i = 0; i < nBest; i++) models[i] = anthill[i];	//needed to ensure proper format of the invidial solution

	//check the settings of pheromone depositing mode:
	ph_sum = 0;
	phUpdate = ph_mode;	
	//validate:
	if ((phUpdate & 3) == 3) phUpdate -= 2;
	if ((phUpdate & 6) == 6) phUpdate -= 4;
	if ((phUpdate & 11) > 8) phUpdate -= 8;
}

bool ACO::Insert(Individual &I)
//seed a starting solution, should be called after Configure()
{
	if ( (mnN > I.body.length()) || (mxN < I.body.length()) ) return false;

	anthill[ RandomN(anthill.length()) ] = I;
	return true;
}

void ACO::update(Individual &I)
//volatility rate should have been applied before update() calls!
{
	REALNUM_TYPE c = I.quality/pow(I.body.length(), beta);

	for (UNSIGNED_2B_TYPE g = 0; g < I.body.length(); g++)
		pheromarks[ I.body[g].GD ] += c;
}

void ACO::resample(Individual &I)
{//TODO: sort by probability when sampling????
	UNSIGNED_2B_TYPE i, j, k, N = I.body.length(), F = pheromarks.length();
	
	if (ph_sum == 0)
	{//calculate cumulative pheromon level
		for (i=0; i < F; i++)	ph_sum += pow(pheromarks[i], alpha);
		if (ph_sum == 0) return;
		
		//sort descriptors
		SIGNED_4B_TYPE *AA = GRAB_MEM_BLOCKS(SIGNED_4B_TYPE, F);
		REALNUM_TYPE *FF	= GRAB_MEM_BLOCKS(REALNUM_TYPE, F);
		QSortScore = FF;
		for (i=0; i < F; i++)
		{
			AA[i] = i;	
			FF[i] = pheromarks[i]; 
		}
		qsort(AA, (size_t)F, sizeof(SIGNED_4B_TYPE), QSortCompareGreater);		
		
		//store sorted descriptors in ipheromarks[]:
		for (i=0; i < F; i++) ipheromarks[i] = AA[i];

		DROP_MEM_BLOCKS(AA);
		DROP_MEM_BLOCKS(FF);		
		QSortScore = NULL;
	}	

	I.quality = INVALID;	//invalidate I's fitness value
	set Picked;
	REALNUM_TYPE pr, prlim, RANDMX = RAND_MAX;
	RANDMX += 1;
	for (j=i=0; i < N; i++)
	{
		do
		{
			prlim  =	rand();
			prlim /=	RANDMX;
			prlim *=	ph_sum;
			pr =	0;

			for (k=0; k < F; k++)
			{	
				j = ipheromarks[k];
				pr += pow(pheromarks[j], alpha);
				if (pr > prlim)
					if (!Picked.IsInSet(j)) break;
			}
		}while (k == F); //repeat if descriptor was not found (e.g. if high rand() value and all top ones were already picked)

		I.body[i].GD = j;
		Picked.PutInSet(j);
	}
}

void ACO::Run(GENERIC_POINTER OutputF)
//performs optimization, calls a user-defined function each cycle
{
	char bfr[256];
	bool ifOutput = false;
	void (*outF)(STRING_TYPE) = (void (*)(STRING_TYPE)) OutputF;	

	if (OutputF != NULL) ifOutput = true;
	
	UNSIGNED_2B_TYPE l, m, n, iter, stb_iter;
	REALNUM_TYPE f, cbest = INVALID;

RERUN:
	iter = stb_iter = 0;

	for(;;)
	{
		iter++;
		
		for (l=0; l < anthill.length(); l++)	anthill[l].EvaluateIndividual(sfunc);	
		SortIndividuals(anthill);
		
		//--------------------------------------------
		//update global best models
		for (m = 0; m < anthill.length(); m++)
		for (l = 0; l < models.length(); l++)		
		{
			if (models[l].quality > anthill[m].quality) continue;
			if (models[l].CompareExact(anthill[m])) break;
			if (models[l].quality < anthill[m].quality) 
			{//if there are multiple equal-quality solutions, all will be checked above first, Oct.2010
				for (n=models.length() -1; n > l; n--) models[n] = models[n-1];
				models[l] = anthill[m];
				break;
			}
		}
		//--------------------------------------------
		
		//compare and update termination parameters
		if (models[0].quality > cbest + converg_thr) 
		{
			cbest = models[0].quality; 
			stb_iter = 0;
		}
		else stb_iter++;
		
		//update pheromon levels
		for (l = 0; l < pheromarks.length(); l++)
			pheromarks[l] *= (1 - p);	//volatility update

		switch (phUpdate & 15)
		{
			case 1:	update(anthill[0]);	break;	//only by iteration best
			case 2: update(models[0]); break;	//only by global best

			case 8:	//ASRank
			{
				n = 1;		//rank
				f = cbest;	//best score
				for (l=0; l < anthill.length(); l++) 
				{//replace quality scores with generated ranks
					if ( (anthill[l].quality + converg_thr) < f)
					{
						n++; 
						f = anthill[l].quality;
					}
					anthill[l].quality = (1.0/n);	//inversely proportional to rank
				}
			} //default: update below is next

			default:
				for (l=0; l < anthill.length(); l++) update(anthill[l]);
			break;
		}

		if (phUpdate & 4) //additional update from the global best
		{
			if (phUpdate & 8) models[0].quality = 1.0;	//ASRank
			update(models[0]);	
			models[0].quality = cbest;
		}

		//reset pheromarks that are too low
		for (l=0; l < pheromarks.length(); l++)
			if (pheromarks[l] < tmin) 
				pheromarks[l] = tmin;

		if (ifOutput)
		{
			sprintf(bfr, "Current cycle: %d; best score: %6.4f", iter, cbest);
			outF(bfr);
		}

		//check termination criteria
		if (stb_iter == stbl_range) break;
		if (iter == nCyc) break;
		if ( ideal_score > INVALID) 
			if (ideal_score == models[0].quality) break;

		//sample new set of solutions
		ph_sum = 0;
		for (l=0; l < anthill.length(); l++) resample(anthill[l]);
	}//for(;;)
	
	if (phUpdate & 16)
	{//after convergence, additional post-optimization with tmax seeding 
		phUpdate -= 16;
		for (l = 0; l < anthill.length(); l++)
		{
			n = (mnN == mxN) ? mnN : (mnN + RandomN(mxN - mnN + 1));
			anthill[l] = Individual(n, pheromarks.length());
		}

		for (l = 0; l < pheromarks.length(); l++) 
			pheromarks[l] = tmax;

		goto RERUN;
	}
}

//***********************************************************************************************
//***********************************************************************************************
//***********************************************************************************************
particle::particle()
{
	q0 = q = INVALID;
}

particle::~particle()
{
	x0.resize(0);
	x.resize(0);
	v.resize(0);
}

void particle::toss(UNSIGNED_2B_TYPE N)
{
	x.resize(N);
	v.resize(N);
	for (UNSIGNED_2B_TYPE i = 0; i < N; i++)
	{
		x[i]  = rand();
		x[i] /= RAND_MAX - 1;

		v[i]  = rand();
		v[i] /= RAND_MAX - 1;
	}

	x0 = x;
	q0 = q = INVALID;
}

void particle::navigate(apvector<REALNUM_TYPE> &Best, REALNUM_TYPE w, REALNUM_TYPE c1, REALNUM_TYPE c2, REALNUM_TYPE Vmax)
//calculates velocity vector
{
	UNSIGNED_2B_TYPE i = 0, N = v.length();
	REALNUM_TYPE f1 = c1*rand(), f2 = c2*rand();
	f1	/= RAND_MAX - 1;
	f2	/= RAND_MAX - 1;

	for (; i < N; i++)
	{
		v[i]	*= w;
		v[i]	+= f1*(x0[i] - x[i]);
		v[i]	+= f2*(Best[i] - x[i]);

		if (v[i] > Vmax) v[i] = Vmax;
		if (v[i] < -Vmax) v[i] = -Vmax;
	}
}

void particle::fly()
//moves the particle
{
	UNSIGNED_2B_TYPE i = 0, N = v.length();
	for (; i < N; i++)	x[i]  += v[i];
	q = INVALID;
}

PSO::PSO()
{
	Configure();
}

PSO::~PSO()
{
	swarm.resize(0);
}

void PSO::Configure(		
		UNSIGNED_2B_TYPE nDims,				//#Features
		UNSIGNED_2B_TYPE nParticles, 		//#particles
		UNSIGNED_2B_TYPE nCycles,			//#cycles
		REALNUM_TYPE 	rtMinDiff,			//Convergence threshold
		UNSIGNED_2B_TYPE	nStabilityRange,//max #cycles without change
		GENERIC_POINTER 	dFunct,			//score function
		REALNUM_TYPE 		rtW,			//w, inertia
		REALNUM_TYPE 		rtC1,			//c1, weight for  cognitive term
		REALNUM_TYPE 		rtC2,			//c2, weight for social term
		REALNUM_TYPE 		rtVmax, 		//Vmax. speed limit
		UNSIGNED_2B_TYPE 	nK)				//k, neighborhood for socializing, 0 means global
{

	ndims = nDims;
	vmax = rtVmax;
	k = nK;
	//weights for inertia, and for cognitive and socializing terms
	w = rtW;
	c1 = rtC1;
	c2 = rtC2;		
	 
	sfunc = dFunct;		//pointer to a score-function	

	nCyc = nCycles;
	converg_thr = rtMinDiff;
	stbl_range = nStabilityRange;	

	swarm.resize(nParticles);
	for (UNSIGNED_2B_TYPE i = 0; i< nParticles; i++)
		swarm[i].toss(ndims);
	
	allbest = INVALID;
}

bool PSO::Insert(Individual &I)
//seed a starting solution, should be called after Configure(). Individual I should have binary genes.
{	
	if ( !swarm.length() ) return false;
	SIGNED_4B_TYPE i, ri = RandomN(swarm.length());

	for (i = 0; i < min(ndims, I.body.length()); i++)
	{
		swarm[ ri ].x[i] = I.body[i].GD; //copy 1 or 0
		swarm[ ri ].x0[i] = swarm[ ri ].x[i];
	}
	
	swarm[ ri ].q = swarm[ ri ].q0 = INVALID;
	return true;
}

bool PSO::Insert(apvector<REALNUM_TYPE> &R)
//seed a starting solution, should be called after Configure()!!
{	
	if ( !swarm.length() ) return false;
	SIGNED_4B_TYPE i, ri = RandomN(swarm.length());

	for (i = 0; i < min(ndims, R.length()); i++)
	{
		swarm[ ri ].x[i] = R[i];
		swarm[ ri ].x0[i] = swarm[ ri ].x[i];
	}
	
	swarm[ ri ].q = swarm[ ri ].q0 = INVALID;
	return true;
}

void PSO::Run(GENERIC_POINTER OutputF)
//performs optimization, calls a user-defined function each cycle
{
	char bfr[256];
	bool ifOutput = true;
	void (*outF)(STRING_TYPE) = (void (*)(STRING_TYPE)) OutputF;	
	REALNUM_TYPE (*DFunction)(apvector<REALNUM_TYPE> &);
	
	if (sfunc == NULL)	return;
	if (OutputF == NULL) ifOutput = false;

	DFunction = (REALNUM_TYPE (*)(apvector<REALNUM_TYPE> &)) sfunc;
	
	REALNUM_TYPE f;	//to calculate convergence
	UNSIGNED_2B_TYPE l, m, n, iter = 0, stb_iter = 0, sl=swarm.length();
	matrix<REALNUM_TYPE> dst;
	if (k) dst.SetSize(sl, sl);
	
	allbest = 0; //set by default to the first particle
	for(;;)
	{
		iter++;
		f = 0;
		for (l=0; l < sl; l++)
		{
			swarm[l].q = DFunction(swarm[l].x);
			f += fabs(swarm[l].q - swarm[l].q0);
			if (swarm[l].q > swarm[l].q0)
			{				
				swarm[l].q0 = swarm[l].q;
				swarm[l].x0 = swarm[l].x;

				if (swarm[l].q > swarm[allbest].q0) allbest = l; //update global best
			}
		}
		f /= sl; //mean absolute difference

		if (k)
		{//calculate distance matrix to find local neighbors
			for (l=0; l < sl; l++)
			{
				dst(l,l) = ndims;	//quite a good approximation for maximum distance
				for (m=l+1; m < sl; m++)
					dst(l,m) = dst(m,l) = getMetricDistance(swarm[l].x, swarm[m].x);
			}
		}

		for (l=0; l < sl; l++) 
		{//calculate velocities
			if (k)
			{
				apvector<SIGNED_4B_TYPE> kn(k+1,l);

				for (m=0; m < sl; m++) 
				{
					n = k;
					while (dst(kn[--n], l) > dst(m, l))
					{
						kn[n+1] = kn[n];
						kn[n] = m;
						if (n == 0) break;
					}
				}

				m = l;
				for (n=0; n < k; n++) 
					if (swarm[m].q0 < swarm[kn[n]].q0) 
						m = kn[n];
			}
			else 
				m = allbest;
			swarm[l].navigate(swarm[m].x0, w, c1, c2, vmax);
		}
		
		//compare and update termination parameters
		if (f >  converg_thr)
			stb_iter = 0;
		else 
			stb_iter++;
	
		if (ifOutput)
		{
			sprintf(bfr, "Current cycle: %d; best score: %6.4f; diff: %8.6f", iter, swarm[allbest].q0, f);
			outF(bfr);
		}

		//check termination criteria
		if (stb_iter == stbl_range) break;
		if (iter == nCyc) break;
		
		for (l=0; l < sl; l++)	swarm[l].fly();			//update positions
	}//for(;;)
	
}
