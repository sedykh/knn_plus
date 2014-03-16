/* knn+.cpp 2008-2013

Console application (Linux, Windows). 
Based on dataset and knn classes

//---------------------------------------------------------------------
Class kNN:
pept1_moe_mdl0 -M=CLS -PNL=0.5 -SEP=0.2 -D=5@40@6 -EVAL=EV0.7@0.7 -LOGALL -OUT=pept1_class
pept1_moe_mdl0 -M=CLS -PNL=0.5 -SEP=0.2 -D=5@40@6 -EVAL=EV0.7@0.7 -LOGALL -OUT=pept1_class
pept1_moe_mdl0 -M=CLS -PNL=0.5 -D=5@40@6 -EVAL=EV0.7@0.7 -LOGALL -OUT=pept1_class

Tried weight schemes:
-WT=EBS.5 -WT=EBQ2.5 -WT=ERS0.5 -WT=KHRS2 -WT=QMB2.2 -WT=NK

PSO:
pept1_moe_mdl0 -M=CTG -PNL=0.5 -D=5@40@5 -O=PS -PS@B=3@N=30@D=100@S=10@V=-2 -SRND=8720 -EVL=0.8@0.8 -LOGALL -OUT=pept1_PS

AC:
pept1_moe_mdl0 -M=CTG -PNL=0.5 -D=5@40@5 -O=AC -WT=NK -LOGALL -OUT=pept1_ACaaa
pept1_moe_mdl0 -M=CTG -PNL=0.5 -D=5@40@5 -O=AC -WT=NK -OUT=pept1_ACaaa -DOMODS -DONEIBS -LOGALL

GA:
pkb20_mdl.t2t -SEED=MID0001_pkb20_mdl.mod@7 -M=CTG -D=5@50 -EVL=0.5@0.4_0.3@0.3 -LOGALL -2OLD -O=GA1 -GA@N=20@D=200@S=10@V=0.01@G=7@E=1@Z=0.05
pkb20_mdl.t2t -SEED=MID0001_pkb20_mdl.mod@7 -M=CTG -D=5@50 -EVL=0.5@0.4_0.3@0.3 -LOGALL -2OLD -O=GA -GA@N=20@D=200@S=10@V=0.01@G=7@E=1

category:
pkb20_mdl.t2t -M=CTG -D=5@7 -EVL=0.5@0.4_0.3@0.3 -SA@N=1@B=3@K=0.6@DT=-1@TE=-2 -LOGALL -2OLD
pkratbin20 -M=CTG -EVL=0.5@0.4_0.3@0.3 -SA@N=2@B=3@K=0.6@DT=-1@TE=-2 -LGO=0.1 -LOGALL

continuous:
KNN+ aa_48 -EVL=0.3@0.1 -SA@N=2@B=3@K=0.6@DT=-1@TE=-2 -LGO=0.1 -LOGALL

//---------------------------------------------------------------------
TODO:


1) implement optimization of RNN (alternative of kNN)
2) add class kNN with non-exclusive memberships? Will be essentially a category MTL (multi-task learning) with adjustable class (i.e., task) weights

NB: for a 2-category model (A and B) of unequal size
CCR > Accuracy if minority class is better predicted than major one
CCR < Accuracy if majority class is better predicted than minor one

//---------------------------------------------------------------------
Versions History:
2.9
June 8 2013				post-optimization mode was expanded to include SA method

2.83
May 27 2013
						minor fix for the running with everything default (would print help instead)
2.82
January 7 2012			-O=XX command key for combinatorial optimization. 
						(old: -O=** gives problem with regular expressions in some scripts)

2.81
November 3 2011			knn::calc_mean_class_separation() added into training evaluation
						class_sep_cf evaluation term added for class-kNN
						its input is enabled through '-SEP=' command line option, help updated.

2.8
October 30 2011			class-mode added (i.e., datapoints have to be with mutually exclusive class memberships)
						added normcdf() to QSAR class to calculate p-values

2.7
October 5-6 2011		q2F13() added to qsar class
						help restructured, 
						new evaluation metric for test sets added: '-EVL=F...'
						MIN_GROUP_N is now also controlling min.size of test set for continuous mode
						in case of continuous models, post-evaluation check is now done only in KNN_EVAL_EXT_R2 mode
2.62
April 11 2011			added logging of the randomizer seed even if not user-specified

2.61
February 20-24	2011	Default parameters of Feature selection algorithms were adjusted
						(based on benchmarking results).
						SA_K =  0.75
						MIN_DIVERSION = 0.001
						POPULATION_SIZE	= 200;
						SWARM_SIZE = 100;
						AC_SIZE = 100;

						minor fix: RAND_MAX+1 possible overflow in ACO_resample()

2.6
February 17 2011		-LOGALL -DOMODS -DONEIBS now properly control output 
						(i.e.logging and secondary files).
						PrintKNNHelp() was expanded

2.5 - 2.52
February 2 2011			Ant Colony solutions size sampling adjusted, 
						to avoid disproportionately small ranges:
						if ( (1 + ((NN_DIMS_MAX - u) << 1)) < NN_DIMS_INC ) 
							u = NN_DIMS_MAX;

January 1 2011
						Neighborhhood weighing functions expanded\overhauled.
						NN_WT_K = 1.0;
						NN_WT_MODE = KNN_WT_K_VOTE | KNN_WT_F_MNK; 
						(these changes correspond to the same default behavior as was set before)

2.42
December 24 2010
						minor fix in ClearDuplIndividuals() to use exact match if GA individual's genes are id-coded
						explicit Gene::Gene(const Gene &) copy-constructor added to feature_alg.h/.cpp
2.41
December 17 2010		minor fix in the solution size cycle for Ant Colony Optimization
2.4
December 7 2010			Particle Swarm Optimization (PSO) added.

2.34
December 2 2010			Sampling improvement for Ant Colony
2.3 - 2.33
October 2010			Ant Colony optimization added, genalgorithm.cpp renamed into feature_alg.cpp
						random number generator overhauled: GetRandomNumber() in core class
						error-based calculation for classes fixed in get_ccr() in qsar class

2.2 - 2.21
July 7 2010				minor fix for -KR= option to allow NN_K_MIN = NN_K_MAX,
						and likewise for NN_DIMS_INC, NN_DIMS_MIN = NN_DIMS_MAX,

June 28 2010			fix in qsar::get_ccr() for the cost-functions based on av.err and max.err.
						(debugged for 3-bins with -EVL=E0.6@0.65 and -EVL=EV0.6@0.65 options)

2.1-2.11
June 2 2010
						minor fix for the warning message:
						sprintf(sbuffr, "'Act.group = %d' has < %d datapoints in test set.", ni, MIN_GROUP_N);

April 29 2010			
						minor fix in PrintOutPredictionsHeader()
						negative q2-R2 evaluation parameters are now allowed by using fabs(rtX) comparison

2.0
April 20 2010			Major restructuring in the workflow and combinatorial functions.
						Post-optimization of created models is enabled (.tbl file as input)

1.43
April 10 2010			MIN_GROUP_N introduced to control coverage of test set categories

1.42
April 8					critical fix for reporting neibors
April 3 2010			dataset-class update

1.37-1.40
March 24-30 2010		Prediction-output file-name adjusted. Case sensitive input filename allowed.

1.36
Feb 28 2010				dataset-class updated (normalization functions, etc.)

1.31 - 1.35

Jan 21 2010				fix of the ccr with penalty, was
						ifPenlt = (PN != 0) //before was: ifPenlt = (PN == 0)

Jan 15-16 2010				
						log-filename modified to be output-specific
						(this avoids multiple runs writing to the same log)
						in SimAnneal.Run() a check was added to make
						#mutated genes at least 1 regardless of mutation rate

						fixed infinite loop possibility in SA knn:
						for z ... ; z = max(1, min(NN_DIMS_INC, NN_DIMS_MAX - z))
						before it was z = min(NN_DIMS_INC, NN_DIMS_MAX - z)

						fix in the R2 filter (in the last check it used q2 instead)

						help and report info was updated

1.30
Jan 10 2010				- overhauling source file system to make separate
						  class files reusable between multiple projects:
						  knn+, datasplit, subgrafs.
						  Stack class and bonds.h were therefore included.
1.20-1.21
Dec 18-21 2009			- LoadDataset() altered to allow loading of x-files only,
						  load() in dataset is altered to skip empty lines.
						  
						  Fixed truncating of output filename extension.
						  
						  New command key -LAXDIMS lets descriptors of a predicted file 
						  be in different order, but labels must match

						  Neighbors are reported for each individual prediction in a
						  special file: *.neibs

1.10-1.14
Dec 8 2009				- extra check for dims mismatch b/w model and test
						  fix in get_dscr_pos() of dataset-class 
Dec 3  2009				- fix of descriptor check in load() in dataset class
Dec 1  2009				- minor fix of x-file loading when a-file is missing
Nov 18 2009				- overhaul of applicability domain options
Nov 10 2009				- update of dataset-class (svm support, etc.)

1.09
Oct 25 2009
						- fixed reading of xa-based result tables
						- more output info added to the log file (version, prediction mode)
1.08
Sept 11 2009 
						- user controlled seed for randomizer added -SRND=

1.06-1.07
July 23-26 2009
						- fixing max.pred.act. overflow (due to rounding)
						- minor adjustments to the external prediction mode
1.05
July 20-21 2009
						- Overhauling prediction-mode
						- alterations of the output and default settings:
							NN_WT_MODE = KNN_WT_K_VOTE (was 0)
							SA_K = 0.8 (was 0.9)
1.04
July 17 2009			- modified Crossover() for variable-size GA
							any legal alignment of parents (before only head to head)					
						- minor fixes

1.00-1.03
July 14 2009			- SizeFnc4GA() added into the GA-based fitness functions
						  without size-based adjustments GA is biased toward hi-size solutions
						  (in NN_GA mutation->increase and crossover->averaging of sizes)
						- penalty factor buf for category kNN fixed

July 13 2009			- SA-scan of dimensions is made by flexible increment (def. 2)

0.99
~ June 2009				- Basic version with all functionalities
*/

#include "knn.h"
#include "feature_alg.h"

#define Version		"2.9"
#define VerDate		"2013"

#define SUBSEP		"@"
#define COMMENT		"#"

#define KNN_LOGFILE		"knn+.log"

#define	NEIBSFILE		".neibs"
#define	PREDSFILE		".preds"
#define	SUMMARYFILE		".tbl"

#define MODELFILE		".mod"
#define PREDFILE		".pred"

#define	SPLITFILE		".t2t"
#define	LISTFILE		".list"

#define	DATAFILE_SVM	".svm"
#define	DATAFILE_NEW	".xa"
#define	DATAFILE_OLD	".x"

//definitions for subspace selection methods
#define NN_STEPWISE1		5	//down
#define NN_STEPWISE			4	//up
#define NN_FULL_DIMSEARCH	1

#define NN_SA				0
#define NN_GA				2
#define NN_GA1				3

#define NN_AC				6	//ant colony
#define NN_PS				7	//particle swarm

//
#define	MIN_GROUP_N			4	//minimum 4 datapoints per class to make reliable stats, also used for overall min.size of test set in continuous models

dataset datasetNN, datasetTRN, datasetTST;
knn knnCore;
GenAlgorithm gaNN;
SimAnneal saNN;
ACO acNN;
PSO psNN;

STRING_TYPE EvalP1 = "q2", EvalP2 = "R2"; //names of evaluation parameters for logging
char sbuffr[256];
STRING_TYPE stOutputFile, stOutputPredFile, stOutputNeibFile, stInputFile, stScreenedFile;
FILETYPE_IN fiNNinput;
apvector<SIGNED_4B_TYPE> dmlist;	//service var to explicate set-content 
apvector<STRING_TYPE> seedims;		//seeded dimensions for a solution to try
UNSIGNED_4B_TYPE nKnnModelId = 0;	//count of generated kNN models

UNSIGNED_4B_TYPE srand_seed = 0;
bool NewInputFormat = true;
bool RunPredictionMode = false, OldFormatOutput = false;

// -------------  general parameters of kNN    ------------------------//
UNSIGNED_1B_TYPE NN_DIMSEL_METHOD = 0;	//GA, SA, full exploration, ets.
UNSIGNED_1B_TYPE NN_APPROX = 0;			//lossy searching of nearest neighbours
bool NN_REPORT_TERSE = true;			//suppresses excessive output from optimizations
bool NN_REPORT_MODSIM = false;			//to report analysis of models (inter-similarity)
bool NN_REPORT_NEIBS = false;			//to report neibors for each prediction

UNSIGNED_2B_TYPE NN_DIMS_MIN = 5, NN_DIMS_MAX = 50, NN_DIMS_INC = 2;//SOLUTION SIZE LIMITS
UNSIGNED_2B_TYPE NN_K_MIN = 1, NN_K_MAX = 9;
REALNUM_TYPE NN_R = 0.0;
UNSIGNED_2B_TYPE NN_WT_MODE = KNN_WT_K_VOTE | KNN_WT_F_MNK; 
REALNUM_TYPE NN_WT_K = 1.0;

//distance metrics
UNSIGNED_1B_TYPE METRIC_K = 0; //Euclidean
REALNUM_TYPE METRIC_V = 2.0;

//Applicability Domain
REALNUM_TYPE NN_AD_ZCUT = 0.5;
UNSIGNED_1B_TYPE NN_AD_MODE = 0; //corresponds to AD(Distance^2) and that at least 1 neighbor required for prediction

//prediction & evaluation
bool NN_PRED_DIM_STRICT = true; //control dimensions check in prediction, iweather labels and positions should match exactly
UNSIGNED_2B_TYPE NN_PRED_MODE = 0, NN_EVAL_MODE = KNN_EVAL_EXT_R2, NN_FILTER = 1;
GENERIC_POINTER NN_PRED_FUNC = NULL, NN_EVAL_FUNC = NULL;
apvector<REALNUM_TYPE> NN_FILTER_CCR;	//to store filtering criteria for CCR

UNSIGNED_2B_TYPE NN_EVAL_LGO = 1;	//#datapoints to leave out during cross-validation prediction
REALNUM_TYPE NN_EVAL_LGO_F = 0;		//to specify NN_EVAL_LGO by train-set fraction, not used by default

REALNUM_TYPE NN_EVAL_Q2 = 0.5, NN_EVAL_R2 = 0.6;
//extra parameters for discrete activity scale
apvector<REALNUM_TYPE> NN_EVAL_BPS(1, 0.5), NN_EVAL_WTS(2, 0.5);
REALNUM_TYPE NN_EVAL_PENALTY = 0.0;	//penalty for imbalanced accuracy
REALNUM_TYPE NN_CLASSSEP_ADJ = 0.0;	//weight of the class-separation term in the model evaluation

REALNUM_TYPE GA_SIZE_ADJ = 0.1;	//basis for the GA score size-related adjustements
//----/-------    parameters for the Genetic Algorithm    -------/-----//
 UNSIGNED_4B_TYPE POPULATION_SIZE			= 200;
 UNSIGNED_4B_TYPE MAX_N_GENERATIONS			= 300; 
 REALNUM_TYPE MUTATION_RATE					= 0.7;
 REALNUM_TYPE CROSSOVER_CHANCE				= 0.8;
 
 REALNUM_TYPE MIN_DIVERSION					= 0.001;
 UNSIGNED_4B_TYPE STABILITY_RANGE			= 10;			//If diversion is minimal throughout this range of generations then quit
 REALNUM_TYPE IDEAL_FIT						= 1.0;			//Ideal fitness 
 
//modes
 UNSIGNED_4B_TYPE CROSSOVER_MODE			= UNIFORM;
 UNSIGNED_4B_TYPE PARENTSELECTION_MODE		= TOURNAMENT;  
 bool ELITISM_MODE							= true;

//extra parameters for some modes
 UNSIGNED_4B_TYPE GROUP_SIZE				= 3;
 REALNUM_TYPE ELITEPART_SIZE				= 0.01;
// --------------------------------------------------------------------//

//----/-------   parameters for the Simulated Annealing only-------/-----//
 UNSIGNED_2B_TYPE SA_NRUNS = 3;					//number of SA runs "from scratch for each configuration and split (in a way it ~ to GA's population size!)
 UNSIGNED_2B_TYPE SA_NBEST2STORE = 3;			//number of best models to store
 REALNUM_TYPE SA_MUTATION_RATE_PER_GENE = 0.2;	//probability of mutation per gene of solution
 UNSIGNED_4B_TYPE SA_N_TRIALS = 100;			//number of mutation trials to do before changing temperature
 bool SA_FULL_ITER = false;						//mode to continue with trials after accepted mutation
 REALNUM_TYPE SA_T_0 = 2, SA_T_END = -5, SA_K =  0.75, SA_T_CONV = -2;
// --------------------------------------------------------------------//

//----/----------   parameters for the Ant Colony only ----------/-----//
 REALNUM_TYPE AC_P = 0.5, AC_A = 1.0, AC_B = 1.0, AC_TMIN = 0.001, AC_TMAX = 1.0;
 UNSIGNED_1B_TYPE AC_PH = 0;
 UNSIGNED_2B_TYPE AC_NBEST = 50;
 UNSIGNED_4B_TYPE AC_SIZE = 100;
// --------------------------------------------------------------------//

//----/----------   parameters for the Particle Swarm only ----------/-----//
UNSIGNED_4B_TYPE SWARM_SIZE			= 100;
REALNUM_TYPE PS_W = 0.7, PS_C1 = 1.5, PS_C2 = 1.5, PS_VMAX =10.0;
UNSIGNED_2B_TYPE PS_K = 0, PS_TOPN = 0; //PS_TOPN is internally used and initialized in a cycle

//predeclarations
void ProcessArgumentString(STRING_TYPE &S);
void PrintParametersInLOG();
void PrintMoreHelp(STRING_TYPE);
void PrintKNNHelp(bool = false);
void GetSeededModel( set&);
void PrintOutputHeader(STRING_TYPE, STRING_TYPE, bool = false);
void PrintOutPredictionsHeader(STRING_TYPE, STRING_TYPE = "");
void PrintOutPredictionsHeader(STRING_TYPE, dataset &, STRING_TYPE = "", bool = false);

void OutputOptProgress(STRING_TYPE);

void OutputkNNPredictions(STRING_TYPE, dataset &, apvector<REALNUM_TYPE>&, REALNUM_TYPE , UNSIGNED_4B_TYPE, bool = true, bool = false);
void OutputkNNPredictions(STRING_TYPE, knn &, apvector<REALNUM_TYPE> &, REALNUM_TYPE , UNSIGNED_4B_TYPE, bool = false); //version for external predictions

void OptimizekNNByK(knn &, UNSIGNED_2B_TYPE, UNSIGNED_2B_TYPE);
bool StepwiseKNN(knn &, UNSIGNED_2B_TYPE, UNSIGNED_2B_TYPE, bool = false);
void CopyDims2kNN(apvector<class Gene> &);
void CopyDims2Genes(set &see, apvector<class Gene> &, UNSIGNED_1B_TYPE = 0);
bool LoadDataset(dataset &, STRING_TYPE, STRING_TYPE, bool AfileNeeded = true);
bool LoadDatasetOLD(apvector<STRING_TYPE>&, dataset&, dataset&, dataset&);
bool LoadFromSplit(STRING_TYPE &,  STRING_TYPE &, dataset &, dataset &, dataset &);
void ReadModelFile(FILETYPE_IN &, REALNUM_TYPE &, UNSIGNED_2B_TYPE &, set &, apvector<STRING_TYPE> &);
STRING_TYPE SkipFileComments(FILETYPE_IN &);
bool VerifyDescriptors(apvector<STRING_TYPE> &, set &, dataset &, dataset &, bool = true);
void DescriptorAnalysis(STRING_TYPE, dataset &, apvector_set_type &);

bool knnFilter(STRING_TYPE);
REALNUM_TYPE SizeFnc4GA(REALNUM_TYPE, SIGNED_4B_TYPE);
REALNUM_TYPE CalcKNNFitGA(apvector<class Gene> &);
REALNUM_TYPE CalcKNNFitGA1(apvector<class Gene> &);
REALNUM_TYPE CalcKNNFitSA(apvector<class Gene> &);
REALNUM_TYPE CalcKNNFitAC(apvector<class Gene> &);
REALNUM_TYPE CalcKNNFitPS(apvector<REALNUM_TYPE> &);
// --------------------------------------------------------------------//

int main(int argc, char* argv[])
{
	srand_seed = (UNSIGNED_4B_TYPE)time(NULL);
	srand(srand_seed);	//default random seed

	if (argc < 2)
	{
		PrintKNNHelp();
		return 0;
	}

	stInputFile = argv[1];

	stOutputFile = stInputFile;
	stOutputFile.tolowercase();

	if (argc == 2)
	{//extended help
		if ( (stOutputFile.find("-o=") == 0) || (stOutputFile.find("-evl") == 0) )
		{
			PrintMoreHelp(stOutputFile);
			return 0;
		}
		else
			if (stOutputFile.find("/?") >= 0)
			{
				PrintKNNHelp( true );
				return 0;
			}
	}
	
	//service variables	
	SIGNED_4B_TYPE nArg = 1, u, z;
	STRING_TYPE stArg, stJ = "command line: knn+ ";
	stJ += argv[nArg];

	//---------------------------------------------------------
	while (argc > ++nArg)
	{ 
		stArg = argv[nArg];
		stJ += BLANK;
		stJ += stArg;
		ProcessArgumentString(stArg);
	}//while (argc > ++nArg) loop

	//---------------------------------------------------------
	//post-processing handling
	if (NN_PRED_MODE)
	{//class or category
		NN_EVAL_MODE -= (NN_EVAL_MODE & KNN_EVAL_EXT_R2);
		if (NN_PRED_MODE == KNN_PRED_CLASS) NN_EVAL_MODE  |= KNN_EVAL_CLASS_ERR;

		if ( NN_FILTER_CCR.length() < NN_EVAL_WTS.length() )
			NN_FILTER = 0;
		if (NN_FILTER == 0)	NN_FILTER_CCR.resize(0);
	}
	else
	{//continuous kNN, remove unnecessary settings
		NN_FILTER_CCR.resize(0);
		NN_EVAL_MODE -= (NN_EVAL_MODE & KNN_EVAL_AVERR);
		if ((NN_EVAL_MODE & KNN_EVAL_EXT_R2) == 0) 		
			NN_FILTER = 0; //do not do post-evaluation check for test-metrics other than R2
	}

	if (NN_DIMSEL_METHOD == NN_AC)
	{//remove redundant modes in ACO settings (for the log-reporting sake)
		if ((AC_PH & 3) == 3) AC_PH -= 2;
		if ((AC_PH & 6) == 6) AC_PH -= 4;
		if ((AC_PH & 11) > 8) AC_PH -= 8;
	}
	//---------------------------------------------------------


	//--------------------------------------------------------
	//prepare progress report labels for modeling mode:
	if (NN_PRED_MODE)
		EvalP1 = EvalP2 = "CCR";
	else
	{
		if ((NN_EVAL_MODE & KNN_EVAL_ERR) == KNN_EVAL_ERR)
			EvalP1 = "MAEq";

		if ((NN_EVAL_MODE & KNN_EVAL_ALT) == KNN_EVAL_ALT)
			EvalP1 += "(alt.)";

		if ((NN_EVAL_MODE & KNN_EVAL_EXT_R2) == 0)
			EvalP2 = EvalP1;
	}
	//--------------------------------------------------------


	if (!CheckStrEnding(stInputFile, LISTFILE) && 
		!CheckStrEnding(stInputFile, SPLITFILE) &&
		!CheckStrEnding(stInputFile, SUMMARYFILE))
	{//if no extension then by default assume a new format
		CutStrEnding(stInputFile);
		if (RunPredictionMode) stInputFile += SUMMARYFILE; else stInputFile += SPLITFILE;
		FILETYPE_IN fiTest(stInputFile.c_str());
		if (fiTest.eof() || fiTest.fail())
		{//revert to old format
			CutStrEnding(stInputFile);
			stInputFile += LISTFILE;
		}
		fiTest.close();
	}

	fiNNinput.open(stInputFile.c_str());
	if (fiNNinput.eof() || fiNNinput.fail())
	{
		cout << "Can not open input file: " << stInputFile << endl;
		return -1;
	}

	//explicit report (in case few files with the same core but different extension exist)
	cout << "Processing input file: '" << stInputFile << "' ..." << endl; 

	NewInputFormat = true;	
	if (CheckStrEnding(stInputFile, LISTFILE))	NewInputFormat = false;

	//-----handle prediction output file name------------------
	stOutputPredFile = stOutputFile;	
	if (stScreenedFile.length() > 0)	stOutputPredFile += "_vs_" + stScreenedFile;	
	stOutputPredFile.tolowercase();
	if (	CheckStrEnding(stOutputPredFile, DATAFILE_NEW) ||
			CheckStrEnding(stOutputPredFile, DATAFILE_OLD)	)
			CutStrEnding(stOutputPredFile);
	stOutputPredFile += PREDSFILE;	
	//-----end of handling prediction output file name------------------
	

	//------------------------------------------------------------------------
	//------------------------------------------------------------------------
	ModulePath		= "";			//use current folder	
	bool ModelInput, SummaryFileInput = CheckStrEnding(stInputFile, SUMMARYFILE);
	ModelInput = SummaryFileInput;	//the case of list-files will be resolved later

	//construct log-file and output-file names
	if ( CheckStrEnding(stOutputFile, SUMMARYFILE)	)	CutStrEnding(stOutputFile);
		
	LOG_FILENAME	= stOutputFile + "_" + KNN_LOGFILE;	
	if ((!RunPredictionMode) && ModelInput)
		if (stInputFile.substr(0, stInputFile.find(SUMMARYFILE)) == stOutputFile) 
			stOutputFile += "_mdf";

	stOutputFile += SUMMARYFILE;
	

	//begin reporting into the log-file
	PutInLogFile("\n\n***********************************");
	stArg = "Current knn+ version: ";
	stArg += Version;
	PutInLogFile(stArg);
	
	GetTimeStamp(stArg);
	PutInLogFile("Starting knn+ on " + stArg);
	PutInLogFile(stJ);

	PrintParametersInLOG();
	//----------- end of handling the settings
	
	REALNUM_TYPE CurrentModelQuality;	//will be used for the optimization of existing models
	SIGNED_2B_TYPE splitcount = 0;	//count of training sets
	apvector<STRING_TYPE> stSplits;
	bool splitloaded = true;
	stOutputNeibFile = stOutputPredFile + NEIBSFILE;

	if (RunPredictionMode)
	{//prediction mode
		PutInLogFile("\n\nPrediction Mode is ON...");
		datasetTST.dump();
		
		if (stScreenedFile.length())
		{//if prediction file is specified...
			if (!CheckStrEnding(stScreenedFile, DATAFILE_NEW) &&
				!CheckStrEnding(stScreenedFile, DATAFILE_OLD))
			{
				CutStrEnding(stScreenedFile);
				stScreenedFile += DATAFILE_NEW;
				FILETYPE_IN fiTest(stScreenedFile.c_str());
				if (fiTest.eof() || fiTest.fail())
				{//revert to old format
					CutStrEnding(stScreenedFile);
					stScreenedFile += DATAFILE_OLD;
				}
				fiTest.close();
			}

			if (CheckStrEnding(stScreenedFile, DATAFILE_OLD))
			{
				stArg = stScreenedFile;
				CutStrEnding(stArg);
				stArg += ".a";
			}

			if (!LoadDataset(datasetTST, stScreenedFile, stArg, false))
			{
				cout << "Can not open data-file for predicting: " << stScreenedFile << endl;
				return -1;
			}	
			
			PrintOutPredictionsHeader(stOutputPredFile, datasetTST, stJ);
			if (NN_REPORT_NEIBS)	PrintOutPredictionsHeader(stOutputNeibFile, datasetTST, stJ, true);
		}; //if (stScreenedFile.length())
	}//if (RunPredictionMode)..else
	else
	{//prepare for modeling
			PrintOutputHeader(stOutputFile, stJ, NN_PRED_MODE );
			PrintOutPredictionsHeader(stOutputPredFile, stJ);
	}

	stJ = SkipFileComments(fiNNinput);
	
	//-------------read and process data-split
	STRING_TYPE PreviousBase;		//to avoid base-reloading
	apvector_set_type setNNDims;	//models for prediction mode
	set setSeedDscrs, setDscrs;		//seeds and work variables
	Individual indSeed;
	STRING_TYPE stTrnDataset;

WORK_W_SPLIT:

	if (splitcount >= setNNDims.length()) setNNDims.resize(splitcount + 128 );
	CurrentModelQuality = 0;

	if (stJ.length() == 0)
	if (!SummaryFileInput)
	{
		stJ = SkipFileComments(fiNNinput);
		if (stJ.length() == 0) goto WORK_W_CYCLEND;; //probably empty lines at the file's end
	}

	if (NewInputFormat)
	{
		if (!ModelInput)
		{//.t2t-file
			if (stJ == PreviousBase)
			{
				splitloaded = true;
				stTrnDataset = stTrnDataset.substr(0, stTrnDataset.find(stInputFile)-1);
			}
			else
			{			
				SplitString(stJ, BLANK, stSplits);

				if (stSplits.length())	
				{//generate filename of the training dataset for later reporting
					stTrnDataset = stSplits[0];
					if ((stSplits.length() > 1) && (!CheckStrEnding(stSplits[0], DATAFILE_NEW)))
					{
						stTrnDataset += BLANK + stSplits[1];
						splitloaded = LoadDataset(datasetNN, stSplits[0], stSplits[1]);
					}
					else
						splitloaded = LoadDataset(datasetNN, stSplits[0], stSplits[0]);
				}
				else splitloaded = false;
			}

			if (splitloaded) 
			{
				LoadSetAsText(fiNNinput, datasetNN.train);
				LoadSetAsText(fiNNinput, datasetNN.test);
				datasetTRN = datasetNN.get_training_set();
				datasetTST = datasetNN.get_test_set();
			}
			goto WORK_W_LOADED;
		}
		
		//.tbl-file
		fiNNinput >> nKnnModelId;
		fiNNinput >> stJ;
		
		//recreate model's base-data				
		stArg.getline(fiNNinput); //stops at the next tab
		stArg.parse_string();
		stJ += BLANK + stArg;
		
		if (stArg.length() == 0) goto WORK_W_CYCLEND; //probably empty lines at the file's end
		stTrnDataset = stArg = stJ;
			
		ReadModelFile(fiNNinput, NN_R, NN_K_MIN, setNNDims[splitcount], stSplits);
		if (!LoadFromSplit(stArg, PreviousBase, datasetNN, datasetTRN, datasetTST))
		{
			PutInLogFile("Can't load model's base-data: " + stArg);					
			goto WORK_W_CYCLEND;
		}
	}
	else
	{//.list-file containing names of .mod-files or splits
		nKnnModelId++;
		STRING_TYPE stSpl;
		SplitString(stJ, BLANK, stSplits);
		u = stSplits.length();
		
		if (u > 5)
		{//split .list format: 'aa_trn0.x aa_trn0.a 42 aa_tst0.x aa_tst0.a 6'
			stTrnDataset = stSplits[0];
			if (!CheckStrEnding(stSplits[0], DATAFILE_NEW))	stTrnDataset += BLANK + stSplits[1];
			splitloaded = LoadDatasetOLD(stSplits, datasetNN, datasetTRN, datasetTST);
			goto WORK_W_LOADED;
		}

		ModelInput = true;
		if (u > 1) 
		{
			stArg = stSplits[--u]; 
			stSpl = stJ + " split#1"; //dummy
		}
		else stArg = stJ;
		
		//read mod-file
		FILETYPE_IN fiMod(stArg.c_str());
		if (fiMod.eof() || fiMod.fail())
		{
			PutInLogFile("Can't find mod-file: " + stArg);
			fiMod.close();
			goto WORK_W_CYCLEND;
		}
		else
		{
			ReadModelFile(fiMod, NN_R, NN_K_MIN, setNNDims[splitcount], stSplits);
			if (stSpl.length() == 0)
			{
				stSpl.getlinewithtabs(fiMod); //line with possible base-data info
				stSpl.parse_string();
			}
		}				
		fiMod.close();
		
		//for the mod-files training set datafile can be identified by: 				
		//a) from the service line after model's description (new format)		
		if (!LoadFromSplit(stSpl, PreviousBase, datasetNN, datasetTRN, datasetTST))
		{//now analyze the name e.g. '48t_a.0_0_2.mod' -> '48t_a.0 48t_a1.0' as the base
			z = u = stArg.find('.');			
			while (++z < stArg.length()) if (stArg[z] == '_') break;
			STRING_TYPE stPfx, stSfx;
			stPfx = stArg.substr(0, u);
			stSfx = stArg.substr(u, z - u) + BLANK + stPfx;
			stSpl = stPfx + stSfx + "1" + stSfx + BLANK + "split#1"; //dummy

			//b) analyzining its name (old format) 
			if (!LoadFromSplit(stSpl, PreviousBase, datasetNN, datasetTRN, datasetTST))
			{
				PutInLogFile("Can't load model's base-data: " + stArg);
				goto WORK_W_CYCLEND;
			}
		}
	}//else ..if (NewInputFormat)

WORK_W_LOADED:

	if (ModelInput)
	{
		CurrentModelQuality = knnCore.qualV;	//save loaded model's quality before reset/initialization
		sprintf(sbuffr, "Current model #%d; ID=%d Qual.=%6.3f", splitcount+1, nKnnModelId, CurrentModelQuality);
	}
	else
	{
		if (RunPredictionMode)
		{
			fiNNinput.close();
			stJ = "Invalid input file for the prediction mode";
			PutInLogFile(stJ);
			cout << stJ << endl;
			return -1;
		}
		sprintf(sbuffr, "Current split: #%d; %s", splitcount+1, stJ.c_str());
	}
	PutInLogFile(sbuffr);


	if (RunPredictionMode)
	{		
		//harmonize decriptor labels vs. specified positions
		if (!VerifyDescriptors(stSplits, setNNDims[splitcount], datasetTRN, datasetTST, NN_PRED_DIM_STRICT)) 
		{
			PutInLogFile("Descriptor labels are inconsistent.. model is skipped!");
			goto WORK_W_CYCLEND;
		}
		
		if (datasetTST.get_Ndatapoints()) //if datasetTST is ready
		{//initialize knn for predictions (dataset, metric, AD, and other prediction-related settings)
			if ((datasetTST.get_Ndimensions() == datasetTRN.get_Ndimensions()) || !NN_PRED_DIM_STRICT)
			{
				knnCore.configure( &(setNNDims[splitcount]), &datasetTRN, 
				NN_K_MIN, NN_R, NN_WT_MODE, NN_WT_K, NN_PRED_MODE, NN_PRED_FUNC, NN_EVAL_MODE, NN_EVAL_LGO, 
				NN_EVAL_FUNC,&NN_EVAL_BPS, &NN_EVAL_WTS, NN_EVAL_PENALTY, NN_CLASSSEP_ADJ,
				METRIC_K, METRIC_V,	NN_AD_ZCUT, NN_AD_MODE, NN_APPROX );

				apvector<REALNUM_TYPE> apClcs;
				apvector<SIGNED_4B_TYPE> apNbs;
				knnCore.evaluate_ext(datasetTST, apClcs, apNbs, false, !NN_PRED_DIM_STRICT);					
				OutputkNNPredictions(stOutputPredFile, knnCore, apClcs, datasetTRN.get_MaxAct(), nKnnModelId);
				if (NN_REPORT_NEIBS)
				{
					knnCore.apvNeibs = apNbs; //to be used inside OutputkNNPredictions()
					OutputkNNPredictions(stOutputNeibFile, knnCore, apClcs, datasetTRN.get_MaxAct(), nKnnModelId, true);
				}
			}
			else
				PutInLogFile(" Descriptor dimensions mismatch b/w test and model sets!");
		}

		splitcount++;
		goto WORK_W_CYCLEND;
	}//if (RunPredictionMode)

	splitloaded = splitloaded && datasetNN.test.Size() && datasetNN.train.Size();

	if (splitloaded)
		PreviousBase = stJ;
	else
	{
		PreviousBase = "";	//invalidate
		PutInLogFile("Loading failed: " + stJ);
		goto WORK_W_CYCLEND;
	}

	if (ModelInput)
	{
		setSeedDscrs = setNNDims[splitcount]; //replace the seed with the currently loaded model
		nKnnModelId *= 10;	 //in case of model optimization scale new ids
		if (NN_DIMSEL_METHOD == NN_SA) //override dimensions if SA will be used for post-optimization
			NN_DIMS_MIN = NN_DIMS_MAX = setSeedDscrs.Size();
	}
	else
	{//to report the split id in the input-file
		if (seedims.length() && setSeedDscrs.IsEmpty()) GetSeededModel(setSeedDscrs);

		sprintf(sbuffr, " %s split#%d", stInputFile.c_str(), splitcount+1);
		stTrnDataset += sbuffr; 
	}
		
	GetTimeStamp(stArg);
	PutInLogFile(stArg);

	if ( (NN_EVAL_LGO_F > 0) && (NN_EVAL_LGO == 1) )
		NN_EVAL_LGO  = max( NN_EVAL_LGO, UNSIGNED_2B_TYPE(NN_EVAL_LGO_F * datasetTRN.get_Ndatapoints()) );
		
	knnCore.configure(NULL, &datasetTRN, 
		NN_K_MIN, NN_R, 
		NN_WT_MODE, NN_WT_K,
		NN_PRED_MODE, NN_PRED_FUNC, 
		NN_EVAL_MODE, NN_EVAL_LGO, 
		NN_EVAL_FUNC,
		&NN_EVAL_BPS, &NN_EVAL_WTS, NN_EVAL_PENALTY, NN_CLASSSEP_ADJ,
		METRIC_K, METRIC_V,
		NN_AD_ZCUT, NN_AD_MODE,
		NN_APPROX);

	z = knnCore.dbase->get_Ndimensions();
	if (NN_DIMSEL_METHOD == NN_GA)
	{
		gaNN.Configure(true, GENE_INTEGER, false, 2,
			0, 0, z, z, z, 
			max(3, POPULATION_SIZE), MAX_N_GENERATIONS, 
			MUTATION_RATE, CROSSOVER_CHANCE, 
			MIN_DIVERSION, STABILITY_RANGE, SizeFnc4GA(IDEAL_FIT, NN_DIMS_MIN), 
			GENERIC_POINTER(&CalcKNNFitGA), CROSSOVER_MODE, PARENTSELECTION_MODE, 
			ELITISM_MODE, GROUP_SIZE, ELITEPART_SIZE);

		//seeding
		for (u = 0; u < gaNN.population.length(); u++)	gaNN.population[u].SetGenesTo(0);
		//above will result in smallest individuals, initially 

		if (!setSeedDscrs.IsEmpty())
		{
			indSeed = gaNN.population[0];
			CopyDims2Genes(setSeedDscrs, indSeed.body, 1);
			gaNN.InsertIndividual(indSeed);
		}

		gaNN.PerformOptimization(NULL, GENERIC_POINTER(&OutputOptProgress));

		for (u = 0; u < gaNN.population.length(); u++)
		{//copy current model into knnCore for final evaluation
			CalcKNNFitGA(gaNN.population[u].body);
			knnFilter(stTrnDataset);
		}//for u
	}//if (NN_DIMSEL_METHOD == NN_GA)

	if (NN_DIMSEL_METHOD == NN_GA1)
	{
		gaNN.Configure(false, GENE_INTEGER, true, z,
			0, 0, NN_DIMS_MAX, NN_DIMS_MIN, NN_DIMS_MAX, 
			max(3, POPULATION_SIZE), MAX_N_GENERATIONS, 
			MUTATION_RATE, CROSSOVER_CHANCE, 
			MIN_DIVERSION, STABILITY_RANGE, SizeFnc4GA(IDEAL_FIT, NN_DIMS_MIN), 
			GENERIC_POINTER(&CalcKNNFitGA1), CROSSOVER_MODE, PARENTSELECTION_MODE, 
			ELITISM_MODE, GROUP_SIZE, ELITEPART_SIZE);

		//seeding
		if (!setSeedDscrs.IsEmpty())
		{
			indSeed = gaNN.population[0]; //NB: setSeedDscrs may get truncated if exceeds NN_DIMS_MIN!!
			CopyDims2Genes(setSeedDscrs, indSeed.body);
			gaNN.InsertIndividual(indSeed);
		}

		gaNN.PerformOptimization(NULL, GENERIC_POINTER(&OutputOptProgress));

		for (u = 0; u < gaNN.population.length(); u++)
		{//copy current model into knnCore for final evaluation
			CalcKNNFitGA1(gaNN.population[u].body);
			knnFilter(stTrnDataset);
		}//for u
	}//if (NN_DIMSEL_METHOD == NN_GA1)

	if (NN_DIMSEL_METHOD == NN_AC)
	for ( z = NN_DIMS_MIN; z < NN_DIMS_MAX; z += max(1, min(NN_DIMS_INC, NN_DIMS_MAX - z)) )
	{//NB: it is expected that z < NN_DIMS_MIN * AC_SIZE
		u = min(z + NN_DIMS_INC - 1, NN_DIMS_MAX);			//to prevent overlap between (z,u) ranges
		if ( (1 + ((NN_DIMS_MAX - u) << 1)) < NN_DIMS_INC ) //to avoid separate sampling of the remaining sizes
			u = NN_DIMS_MAX;

		sprintf(sbuffr, "Ant Colony solutions of sizes #%d-%d", z, u);
		PutInLogFile(sbuffr);

		acNN.Configure(datasetTRN.get_Ndimensions(), AC_SIZE, z, u, MAX_N_GENERATIONS, 
			MIN_DIVERSION, STABILITY_RANGE, IDEAL_FIT, GENERIC_POINTER(&CalcKNNFitAC), 
			AC_P, AC_A, AC_B, AC_TMIN, AC_TMAX, AC_PH, AC_NBEST);
	
		//seeding
		if (!setSeedDscrs.IsEmpty())
		{
			indSeed = acNN.models[0];	//NB: setSeedDscrs may get truncated if exceeds NN_DIMS_MIN!!
			CopyDims2Genes(setSeedDscrs, indSeed.body);
			acNN.Insert(indSeed);
		}

		acNN.Run(GENERIC_POINTER(&OutputOptProgress));

		stJ = "Pheromarks:";
		for (u=0; u < acNN.pheromarks.length(); u++)
		{
			sprintf(sbuffr, "\t%6.4f", acNN.pheromarks[u]);
			stJ += sbuffr;
		}
		PutInLogFile(stJ);

		//process resulting models
		for (u = 0; u < min(AC_NBEST, acNN.models.length()); u++)
		{
			CalcKNNFitSA(acNN.models[u].body); //recalculates predictions
			knnFilter(stTrnDataset);
		}
	}//if (NN_DIMSEL_METHOD == NN_AC) for..

	if (NN_DIMSEL_METHOD == NN_PS)
	for ( z = NN_DIMS_MIN; z <= NN_DIMS_MAX; z += max(1, min(NN_DIMS_INC, NN_DIMS_MAX - z)) )
	{
		PS_TOPN = z;
		sprintf(sbuffr, "Particle Swarm solutions of size #%d", z);
		PutInLogFile(sbuffr);

		psNN.Configure(datasetTRN.get_Ndimensions(), SWARM_SIZE, MAX_N_GENERATIONS, 
			MIN_DIVERSION, STABILITY_RANGE, GENERIC_POINTER(&CalcKNNFitPS), 
			PS_W, PS_C1, PS_C2, PS_VMAX, PS_K);

		//seeding
		if (!setSeedDscrs.IsEmpty())
		{
			indSeed.body.resize(PS_TOPN);
			indSeed.SetGenesTo(); //binary genes will be used for seeding!
			CopyDims2Genes(setSeedDscrs, indSeed.body, 1);
			psNN.Insert(indSeed);
		}

		psNN.Run(GENERIC_POINTER(&OutputOptProgress));

		//process resulting models
		//scan adjacent ranges of descriptors as well		
		for (PS_TOPN = max(z - NN_DIMS_INC, NN_DIMS_MIN); PS_TOPN < min(z + NN_DIMS_INC, NN_DIMS_MAX); PS_TOPN++)
		{
			apvector_set_type sw( psNN.swarm.length() ); 
			for (u = 0; u < psNN.swarm.length(); u++)
			{
				CalcKNNFitPS(psNN.swarm[u].x0); //recalculates predictions
				
				//check for duplicate solutions (e.g. due to swarm collapse)
				SIGNED_4B_TYPE v = 0;
				for (; v < u; v++)									
					if ((sw[v] & knnCore.dims).Size() == PS_TOPN) 
						break;
				if (v == u) sw[u] = knnCore.dims; else continue;
				
				knnFilter(stTrnDataset);
			}//for u
		}//for PS_TOPN
	}//if (NN_DIMSEL_METHOD == NN_PS) for..


	if (NN_DIMSEL_METHOD == NN_SA)
		for ( z = NN_DIMS_MIN; z <= NN_DIMS_MAX; z += max(1, min(NN_DIMS_INC, NN_DIMS_MAX - z)) )
		{
			sprintf(sbuffr, "Simulated Annealing: #%d dimension solutions..", z);
			PutInLogFile(sbuffr);

			for (UNSIGNED_2B_TYPE rn = 0; rn < SA_NRUNS; rn++)
			{
				saNN.Configure(false, datasetTRN.get_Ndimensions(), 0, 1, 
				z, SA_FULL_ITER, SA_NBEST2STORE, GENERIC_POINTER(&CalcKNNFitSA), IDEAL_FIT, 
				SA_MUTATION_RATE_PER_GENE, SA_N_TRIALS,	SA_T_0, SA_T_END, SA_T_CONV, SA_K);
				if (!setSeedDscrs.IsEmpty())
				{
					indSeed = saNN.models[0];
					CopyDims2Genes(setSeedDscrs, indSeed.body);
					saNN.Insert(indSeed);
				}

				saNN.Run(GENERIC_POINTER(&OutputOptProgress));

				for (u = 0; u < min(SA_NBEST2STORE, saNN.models.length()); u++)
				{//copy current model into knnCore for final evaluation
					CalcKNNFitSA(saNN.models[u].body); //recalculates predictions
					knnFilter(stTrnDataset);
				}//for u
			}//for z, rn
	} //if (NN_DIMSEL_METHOD == NN_SA) for ...

	if ( (NN_DIMSEL_METHOD == NN_STEPWISE) || 
		(NN_DIMSEL_METHOD == NN_STEPWISE1) || 
		(NN_DIMSEL_METHOD == NN_FULL_DIMSEARCH) )
	{//define descriptor sets for combinatorial exploration
		if (ModelInput)	
			knnCore.qualV = CurrentModelQuality;	//restore the quality of the current model

		if (setSeedDscrs.IsEmpty())
			setDscrs = set(0, datasetTRN.get_Ndimensions());
		else
			setDscrs = setSeedDscrs;
	}

	if (NN_DIMSEL_METHOD == NN_FULL_DIMSEARCH)
	{
		for (z = NN_DIMS_MIN; z <= NN_DIMS_MAX; z++)		
		if (GetCombination(setDscrs, knnCore.dims, z))
		do
		{
			knnCore.kdist.SetSize(0, 0);
			OptimizekNNByK(knnCore, NN_K_MIN, NN_K_MAX);			
			knnFilter(stTrnDataset);
		} while (GetCombination(setDscrs, knnCore.dims));	
	}//if (NN_DIMSEL_METHOD == NN_FULL_DIMSEARCH)

	if ((NN_DIMSEL_METHOD == NN_STEPWISE) || (NN_DIMSEL_METHOD == NN_STEPWISE1))
	{
		if (NN_DIMSEL_METHOD == NN_STEPWISE)
			knnCore.dims = setSeedDscrs;
		else
			knnCore.dims = setDscrs;

		u = z =	knnCore.dims.Size();
		while ( ((NN_DIMSEL_METHOD == NN_STEPWISE) && (z < NN_DIMS_MAX)) ||
				((NN_DIMSEL_METHOD == NN_STEPWISE1) && (z > NN_DIMS_MIN)) )
		{
			if ( !StepwiseKNN(knnCore, NN_K_MIN, NN_K_MAX, (NN_DIMSEL_METHOD == NN_STEPWISE)) )	break;
			if (knnFilter(stTrnDataset)) u = 0;
			if (NN_DIMSEL_METHOD == NN_STEPWISE) z++; else z--;
		}

		if (ModelInput && (u > 0))
		{//if no improvement, attempt to report the original model
			knnCore.dims = setNNDims[splitcount];
			knnCore.kdist.SetSize(0, 0);
			OptimizekNNByK(knnCore, NN_K_MIN, NN_K_MAX);
			knnFilter(stTrnDataset);
		}
	}
	
	splitcount++;

WORK_W_CYCLEND:	
	stJ = ""; //should be kept empty for the next iteration

	//read another datasplit entry and repeat the work
	if (!fiNNinput.eof())	goto WORK_W_SPLIT;	

	if (RunPredictionMode && NN_REPORT_MODSIM)
	{
		setNNDims.resize(splitcount);		
		DescriptorAnalysis(stInputFile, datasetTRN, setNNDims);
	}

	cout << "Done." << endl;
	return 0;
}

//service functions ---------------------------------
STRING_TYPE SkipFileComments(FILETYPE_IN &fiTXT)
{
	STRING_TYPE stX;
	do 
	{//skip the comments in the file being read
		stX = "";
		stX.getlinewithtabs(fiTXT);
		stX.parse_string();
	} while ( ((stX.find(COMMENT) == 0) || (stX.length() == 0)) && !fiTXT.eof() );
	return stX;
}

void OutputOptProgress(STRING_TYPE S)
{
	if (NN_REPORT_TERSE) return;
	PutInLogFile(S);
}

void CopyDims2kNN(apvector<class Gene> &I)
{
	CopyGenes2Set(I, knnCore.dims);	
	knnCore.kdist.SetSize(0,0);
}

void CopyDims2Genes(set &seeD, apvector<class Gene> &I, UNSIGNED_1B_TYPE Md)
{
	apvector<SIGNED_4B_TYPE> dims;
	SIGNED_4B_TYPE dms;
	seeD.GetList(dims);
	if (Md)	
	{//GA,PS
		for (dms = 0; dms < I.length(); dms++) I[dms].GD = 0;
		for (dms = 0; dms < dims.length(); dms++)	I[dims[dms]].GD = 1;
		return;
	}

	//SA, AC, GA1
	for (dms = 0; dms < min(dims.length(), I.length()); dms++)	I[dms].GD = dims[dms];
}

void OptimizekNNByK(knn &K, UNSIGNED_2B_TYPE KMN, UNSIGNED_2B_TYPE KMX)
{
	REALNUM_TYPE qualB = -1000;
	UNSIGNED_2B_TYPE kB = K.k;	
	for (K.k = KMN; K.k <= KMX; K.k++)
	{
		K.evaluate();	
		if (K.qualV > qualB)	
		{
			qualB = K.qualV;
			kB = K.k;
		}
	}

	K.k = kB;
	K.evaluate(); //to update predictions matrix and qualV
}

REALNUM_TYPE SizeFnc4GA(REALNUM_TYPE rtV, SIGNED_4B_TYPE sz)
{//returns adjusted value for the GA fitness function based on model's size
	REALNUM_TYPE sf = GA_SIZE_ADJ;
	sf	/= 1 + sz - NN_DIMS_MIN;
	sf	+= 1;	//scaling should be better than a constant additive term!
	return (sf * rtV);
}

REALNUM_TYPE CalcKNNFitGA(apvector<class Gene> &I)
{//will randomly adjust the number of genes to comply with constraints

	if (knnCore.dbase == NULL) return 0;
	SIGNED_4B_TYPE i, j, Nm = I.length();
	
	knnCore.dims.Dump();
	for (j = i = 0; i < Nm; i++) 
	if (I[i].GD) 
	{
		knnCore.dims.PutInSet(i);
		j++;
	}
	
	while (j < NN_DIMS_MIN)
	{//randomly add some
		i = GetRandomNumber(Nm);
		if (!knnCore.dims.IsInSet(i))
		{
			knnCore.dims.PutInSet(i);
			I[i].GD = 1;
			j++;
		}
	}

	if (j > NN_DIMS_MAX)
	{//randomly remove some
		knnCore.dims.GetList(dmlist);
		while (j > NN_DIMS_MAX)
		{
			i = GetRandomNumber(dmlist.length());
			knnCore.dims.RemoveFromSet(dmlist[i]);
			I[i].GD = 0;
			j--;
		}
	}
	
	knnCore.kdist.SetSize(0,0); //to clear out previous distance-data!
	OptimizekNNByK(knnCore, NN_K_MIN, NN_K_MAX);
	return (SizeFnc4GA(knnCore.qualV, j)); //returns size-dependent fitness (relevant for optimization only)
}

REALNUM_TYPE CalcKNNFitGA1(apvector<class Gene> &I)
{//will randomly adjust the number of genes to comply with constraints
	if (knnCore.dbase == NULL) return 0;
	SIGNED_4B_TYPE i, j, Nm = I.length();
	
	knnCore.dims.Dump();
	for (j = i = 0; i < Nm; i++)
	{
		if (j < i) I[j].GD = I[i].GD;
		if (knnCore.dims.IsInSet(I[i].GD)) continue;		
		knnCore.dims.PutInSet(I[i].GD);
		j++;		
	}

	if (j < Nm)
	{
		I.resize(max(j, NN_DIMS_MIN));
		while (j < NN_DIMS_MIN)
		{//randomly add some
			while (knnCore.dims.IsInSet(I[j].GD)) I[j].MutateGene();
			knnCore.dims.PutInSet(I[j].GD);
			j++;
		}
	}
	
	knnCore.kdist.SetSize(0,0); //to clear out previous distance-data!
	OptimizekNNByK(knnCore, NN_K_MIN, NN_K_MAX);
	return (SizeFnc4GA(knnCore.qualV, j)); //returns size-dependent fitness (relevant for optimization only)
}

REALNUM_TYPE CalcKNNFitSA(apvector<class Gene> &I)
{
	if (knnCore.dbase == NULL) return 0;
	CopyDims2kNN(I);
	OptimizekNNByK(knnCore, NN_K_MIN, NN_K_MAX);
	return knnCore.qualV;
}

REALNUM_TYPE CalcKNNFitAC(apvector<class Gene> &I)
{//just a wrapper
	//return non-negative only, otherwise pheromone level will be reduced
	return max(CalcKNNFitSA(I), 0);	//NB: to subtract random performance results would require extra handling of prediction modes
}

REALNUM_TYPE CalcKNNFitPS(apvector<REALNUM_TYPE> &X)
//note, for binary use, PSO has redundant sets of solutions: because only relative values matter;
//i.e. all (x + v) are the same, if v is const*x. These redundant solutions can be detected by Spearman rank correlation.
{
	if (knnCore.dbase == NULL) return 0;
	
	knnCore.kdist.SetSize(0,0);
	apvector<SIGNED_4B_TYPE> psd(PS_TOPN + 1,0); //#descriptors to actually use
	UNSIGNED_2B_TYPE i = 1, j = 0;
	for (; i < PS_TOPN; i++) 
	{
		j = i;
		if ( X[ psd[j - 1] ] > X[i] ) 
			psd[j] = i;
		else
			while ( X[ psd[--j] ] < X[i] )
			{
				psd[j+1] = psd[j];
				psd[j] = i;	
				if (j == 0) break;
			}
	}

	for (j = PS_TOPN; j < X.length(); j++)
	{
		i = PS_TOPN;	
		while (X[psd[--i]] < X[j])
		{
			psd[i+1] = psd[i];
			psd[i] = j;
			if (i == 0) break;
		}
	}
	psd.resize(PS_TOPN);
	knnCore.dims = psd;
	
	OptimizekNNByK(knnCore, NN_K_MIN, NN_K_MAX);
	return knnCore.qualV;
}

bool LoadDataset(dataset &D, STRING_TYPE X, STRING_TYPE A, bool AfileNeeded)
{
	bool r = false;
	FILETYPE_IN fiX(X.c_str());
	if (fiX.fail()) return false;

	if (CheckStrEnding(X, DATAFILE_NEW))
		r = D.load(fiX);
	else
	{
		FILETYPE_IN fiA(A.c_str());
		if (fiA.fail())
		{
			if (!AfileNeeded) //load x-file only
				r = D.load(fiX, NULL, 10);
		}
		else
		{//Explicit check; Dec 1 2009
			r = D.load(fiX, &fiA, 1);
			fiA.close();
		}			
	}
	fiX.close();
	return r;
}

bool StepwiseKNN(knn &K, UNSIGNED_2B_TYPE kmin, UNSIGNED_2B_TYPE kmax, bool Up)
/*
description:	performs stepwise exploration of dimensions for kNN
				Up specifies addition of dimensions, !Up - reduction

postcondition:	updated model is in K and returns true
				if no change was possible returns false
*/
{
	SIGNED_4B_TYPE dm, dmB = INVALID;
	UNSIGNED_2B_TYPE kkB = K.k;
	REALNUM_TYPE qualB = K.qualV; //worst imaginable quality
	set D(0, K.dbase->get_Ndimensions());
	if (Up) D -= K.dims; else D = K.dims;

	while (D.GetElement(dm))
	{
		D.RemoveFromSet(dm);
		if (Up)	K.dims.PutInSet(dm); else K.dims.RemoveFromSet(dm);

		K.kdist.SetSize(0, 0);
		OptimizekNNByK(K, kmin, kmax);
		if (qualB < K.qualV) 
		{
			qualB = K.qualV;
			dmB = dm;
			kkB = K.k;
		}
		if (Up)	K.dims.RemoveFromSet(dm); else K.dims.PutInSet(dm);
	};

	K.k = kkB;
	K.qualV = qualB;
	K.kdist.SetSize(0, 0);
	if (dmB == INVALID) return false;	
	if (Up)	K.dims.PutInSet(dmB); else K.dims.RemoveFromSet(dmB);
	K.evaluate();
	return true;
}

void PrintOutPredictionsHeader(STRING_TYPE PF, dataset &dtstA, STRING_TYPE ARG, bool ifNeib)
{//for external predictions
	FILETYPE_OUT foNNpred(PF.c_str());
	SIGNED_4B_TYPE i, pdN = dtstA.get_Ndatapoints();
	
	foNNpred << COMMENT << ARG << endl;
	foNNpred << TAB << "dataID->";
	for (i = 0; i < pdN; i++) foNNpred << TAB << dtstA.get_sid(i);
	foNNpred << endl;	
	
	foNNpred << TAB << "Experim.ACT->";
	for (i = 0; i < pdN; i++) foNNpred << TAB << dtstA.get_Act(i);
	foNNpred << endl;

	foNNpred << "MODEL_ID" << TAB << "AD(distance)" << TAB;
	if (ifNeib) foNNpred << "TST_Neibs->"; else foNNpred << "TST_ClcACT->";
	foNNpred << endl;
	foNNpred.close();
}

void PrintOutPredictionsHeader(STRING_TYPE PF, STRING_TYPE ARG)
{//to output predictions during model training
	FILETYPE_OUT foNNpred(PF.c_str());
	foNNpred << COMMENT << ARG << endl;
	foNNpred << "MODEL_ID" << TAB;
	foNNpred << "TRN_dataID" << TAB << "TRN_ExpACT" << TAB << "TRN_ClcACT" << TAB;
	foNNpred << "TST_dataID" << TAB << "TST_ExpACT" << TAB << "TST_ClcACT" << endl;
	foNNpred.close();
}

void PrintOutputHeader(STRING_TYPE TABFILE, STRING_TYPE ARG, bool knnCateg)
{
	FILETYPE_OUT foNNres(TABFILE.c_str());
	SIGNED_1B_TYPE rep = 2;

	foNNres << COMMENT << ARG << endl;
	foNNres << COMMENT << "| MODEL-INFO || STATISTICS OF TRAINING-SET || STATISTICS OF TEST-SET |" << endl;
	if (knnCateg)
	{
		  foNNres << COMMENT << "STATISTICS OF CLASSIFICATION ACCURACY (CATEGORIES)" << endl;
		  foNNres << COMMENT << "QualityLimit-> i.e. one of CCR-functions;" << endl;		  
	}
	else
	{
		foNNres << COMMENT << "| STATISTICS OF REGRESSION LINES Y = b11*X + b01 AND X = b12*y + b02";
		foNNres << " || STATISTICS OF REGRESSIONS THROUGH ORIGIN Y = k1*X AND X = k2*Y |" << endl;
		foNNres << COMMENT << "QualityLimit-> i.e. q2(training) / R2(test);" << endl;
	}

	foNNres << "ModelID" << TAB << "Datafile" << TAB;
	foNNres << "Ndims" << TAB << "DimsIDs" << TAB << "DimsNames" << TAB;
	foNNres << "k(R)" << TAB;

	while (rep--)
	if (knnCateg)
	{
		foNNres << "QualityLimit" << TAB << "Ndatapoints" << TAB;
		foNNres << "Accuracy" << TAB << "CCR(Normalized Accuracy)" << TAB;
		foNNres << "Accuracy(with group-weights)" << TAB << "CCR(with group-weights)" << TAB;
		foNNres << "Accuracy(Max.Err.based)" << TAB << "CCR(Max.Err.based)" << TAB;
		foNNres << "Accuracy(Av.Err.based)" << TAB << "CCR(Av.Err.based)" << TAB;

		for (SIGNED_4B_TYPE c = 1; c <= NN_EVAL_WTS.length(); c++)
			foNNres << "Ndatapoints_Group" << c << TAB << "Accuracy_Group" << c << TAB;
	}
	else
	{
		foNNres << "QualityLimit" << TAB << "Ndatapoints" << TAB;
		foNNres << "st.dev(Act.)" << TAB << "st.dev(Act.calc.)" << TAB;

		foNNres << "b01" << TAB << "b11" << TAB << "b02" << TAB << "b12" << TAB;
		foNNres << "R" << TAB << "R^2" << TAB;
		foNNres << "MSE1" << TAB << "MSE2" << TAB << "F1" << TAB << "F2" << TAB;

		foNNres << "k1" << TAB << "k2" << TAB;
		foNNres << "R0^2" << TAB << "R01^2" << TAB;
		foNNres << "MSE01" << TAB << "MSE02" << TAB "F01" << TAB << "F02" << TAB;

		foNNres << "q^2" << TAB << "q'^2" << TAB << "MAEq" << TAB << "MAEq'" << TAB;
		foNNres << "MSE" << TAB << "MAE" << TAB;		
	}

	foNNres << endl;
	foNNres.close();
}

void OutputkNNPredictions(STRING_TYPE predFile, knn &knnM, apvector<REALNUM_TYPE> &preds, REALNUM_TYPE maxAct, UNSIGNED_4B_TYPE ModelId, bool ifNeibs)
{//to output series of external predictions
	REALNUM_TYPE rtAppDmn = knnM.calc_AD();

	FILETYPE_OUT foNNpreds(predFile.c_str(), ios_base::app | ios_base::ate);	
	SIGNED_4B_TYPE nn, i, pdN = preds.length();
	foNNpreds << ModelId << TAB << rtAppDmn;

	for (nn = i = 0; i < pdN; i++) 
	{
		foNNpreds  << TAB;
		if (preds[i] > maxAct) foNNpreds << "NA"; 
		else 
			if (!ifNeibs) foNNpreds << preds[i];
			else
				//report nearest neighbors
				while ( (nn + 1 < knnM.apvNeibs.length()) && (knnM.apvNeibs[nn] == i) )
				{
					foNNpreds << (knnM.apvNeibs[++nn] + 1);
					foNNpreds << "{" << (knnM.dbase->get_sid(knnM.apvNeibs[nn])) << "} ";
					nn++;
				}
	}
	foNNpreds << endl;
	foNNpreds.close();
}

void OutputkNNPredictions(STRING_TYPE predFile, dataset &dtstA, apvector<REALNUM_TYPE> &preds, REALNUM_TYPE maxAct, UNSIGNED_4B_TYPE ModelId, bool HeadLine, bool EndLine)
{//to output preditions during model training
	FILETYPE_OUT foNNpreds(predFile.c_str(), ios_base::app | ios_base::ate);	
	SIGNED_4B_TYPE i, pdN = preds.length();
	if (HeadLine) foNNpreds << ModelId << TAB;
	for (i = 0; i < pdN; i++) foNNpreds << dtstA.get_sid(i) << BLANK;
	foNNpreds << TAB;
	for (i = 0; i < pdN; i++) foNNpreds << dtstA.get_Act(i)  << BLANK;
	foNNpreds << TAB;	
	for (i = 0; i < pdN; i++) 
	{
		if (preds[i] > maxAct) foNNpreds << "NA"; else foNNpreds << preds[i];
		foNNpreds << BLANK;
	}

	if (EndLine)  foNNpreds << endl; else foNNpreds << TAB;
	foNNpreds.close();
}

void OutputkNNResults(STRING_TYPE MasterFile, STRING_TYPE DataFile, knn &Model, UNSIGNED_4B_TYPE ModelId, bool ifModFile = false, bool SavePreds = false)
//if input was a .t2t-file then DataFile should have at the end the split number referring to that file.
{
	FILETYPE_OUT foNNres(MasterFile.c_str(), ios_base::app | ios_base::ate);
	if (ifModFile)
	{//save the model in an individual file
		STRING_TYPE stModelname = DataFile.substr(0, DataFile.find(BLANK));
		CutStrEnding(stModelname);
		sprintf(sbuffr, "MID%.4d_%s%s", ModelId, stModelname.c_str(), MODELFILE);
		stModelname  = sbuffr;
		foNNres.close();
		foNNres.open(stModelname.c_str());
	}
	else
		foNNres << ModelId << TAB << DataFile << TAB;

	//save the rest of the model - data, including k/R etc.
	Model.dims.GetList(dmlist);
	foNNres << dmlist.length();
	if (ifModFile)	foNNres << endl; else	foNNres << TAB;

	SIGNED_4B_TYPE i;
	foNNres << dmlist[0];
	for (i = 1; i < dmlist.length(); i++) foNNres << BLANK << dmlist[i];
	if (ifModFile)	foNNres << endl; else	foNNres << TAB;

	foNNres << Model.dbase->get_dscr(dmlist[0]);
	for (i = 1; i < dmlist.length(); i++) foNNres << BLANK << Model.dbase->get_dscr(dmlist[i]);
	if (ifModFile)	foNNres << endl; else	foNNres << TAB;

	if (Model.k == 0)
	{
		sprintf(sbuffr, "%4.3f", Model.r);
		foNNres << sbuffr;
	}
	else
		foNNres << Model.k;

	if (!ifModFile)	
	{
		foNNres << TAB;	
		foNNres.close();
		return;
	}

	foNNres << BLANK << Model.qualV << endl;	
	foNNres << DataFile << endl;

	if (SavePreds)
	{//print-out self-predictions, needed for compatibility
		foNNres << "id" << BLANK << "Calc.Act." << BLANK << "Act." << endl;
		for (SIGNED_4B_TYPE i = 0; i < Model.pred_data.length(); i++)
		{
			foNNres << Model.dbase->get_sid(i) << BLANK;
			foNNres << Model.pred_data[i] << BLANK;
			foNNres << Model.dbase->get_Act(i)  << endl;
		}
	}

	foNNres.close();
}

void printStats4knnContinuous(STRING_TYPE MasterFile, knn &K, apvector<REALNUM_TYPE> &dpExp, apvector<REALNUM_TYPE> &dpClc, bool EndLine = false)
{
	FILETYPE_OUT foNNres(MasterFile.c_str(), ios_base::app | ios_base::ate);
	
	//quality-index
	foNNres << K.qualV << TAB;
	
	//#points
	SIGNED_4B_TYPE it, ndp = dpClc.length();
	foNNres << ndp << TAB;
	
	//standard deviation values for Experimental (Y) and Calculated (X) Activities
	REALNUM_TYPE sE =K.qsarBLOCK.stdev(dpExp), sC = K.qsarBLOCK.stdev(dpClc);
	foNNres << sE << TAB << sC << TAB;

	//Y = A*X + B;	X = A1*Y + B1	
	REALNUM_TYPE A, B, A1, B1, R;
	R = K.qsarBLOCK.correl(dpExp, dpClc);
	A = K.qsarBLOCK.trendline(dpExp, dpClc, B);
	A1 = K.qsarBLOCK.trendline(dpClc, dpExp, B1);
	
	//4 coefficients, R and R^2
	foNNres << B << TAB << A << TAB << B1 << TAB << A1 << TAB;
	foNNres << R << TAB << sqr(R) << TAB;
		
	apvector<REALNUM_TYPE> Yx (ndp), Xy (ndp);
	for (it = 0; it < ndp; it++)
	{
		Yx[it] = A*dpClc[it] + B;
		Xy[it] = A1*dpExp[it] + B1;
	};
	
	REALNUM_TYPE F, F1, MSE, MSE1, MSS, MSS1;
	MSS = K.qsarBLOCK.varianceV(dpExp, true);
	MSE = K.qsarBLOCK.MSE(Yx, dpExp);
	
	F = K.qsarBLOCK.Fvalue(MSS, MSE, 1, ndp - 2); 

	MSS1 = K.qsarBLOCK.varianceV(dpClc, true);
	MSE1 = K.qsarBLOCK.MSE(Xy, dpClc);
	F1 = K.qsarBLOCK.Fvalue(MSS1, MSE1, 1, ndp - 2); 

	//2 MSEs and 2 F-values
	foNNres << MSE << TAB << MSE1 << TAB << F << TAB << F1 << TAB;

	//trendline through the origin
	//Y = A*X;	X = A1*Y
	REALNUM_TYPE kR2, k1R2;
	A	 = K.qsarBLOCK.trendline0(dpExp, dpClc);
	A1	 = K.qsarBLOCK.trendline0(dpClc, dpExp);	

	kR2	 = K.qsarBLOCK.sqrR0(dpExp, dpClc);
	k1R2 = K.qsarBLOCK.sqrR0(dpClc, dpExp);
	
	//2 cofficientss and 2 R^2
	foNNres << A << TAB << A1 << TAB << kR2 << TAB << k1R2 << TAB;

	for (it = 0; it < ndp; it++)
	{
		Yx[it] = A*dpClc[it];
		Xy[it] = A1*dpExp[it];
	};

	MSE = K.qsarBLOCK.MSE(Yx, dpExp);
	MSS = K.qsarBLOCK.varianceVext(Yx, K.qsarBLOCK.meanV(dpExp), true);
	MSS += MSE;
	F = K.qsarBLOCK.Fvalue(MSS, MSE, 1, ndp - 2); 

	MSE1 = K.qsarBLOCK.MSE(Xy, dpClc);
	MSS1 = K.qsarBLOCK.varianceVext(Xy, K.qsarBLOCK.meanV(dpClc), true);
	MSS1 += MSE1;
	F1 = K.qsarBLOCK.Fvalue(MSS1, MSE1, 1,  ndp - 2);

	//2 MSEs and 2 F-values
	foNNres << MSE << TAB << MSE1 << TAB << F << TAB << F1 << TAB;

	REALNUM_TYPE q2, q2_1, MAEq, MAEq_1, MAE;
	q2	 = K.qsarBLOCK.q2etc(dpExp, dpClc);
	q2_1 = K.qsarBLOCK.q2etc(dpClc, dpExp);
	MAEq = K.qsarBLOCK.q2etc(dpExp, dpClc, 2);
	MAEq_1 = K.qsarBLOCK.q2etc(dpClc, dpExp, 2);
		
	MSE = K.qsarBLOCK.MSE(dpExp, dpClc);
	MAE = K.qsarBLOCK.RSS(dpExp, dpClc, true) / ndp;

	//4 q2-like parameters; MSE and MAE
	foNNres << q2 << TAB << q2_1 << TAB << MAEq << TAB << MAEq_1 << TAB;
	foNNres << MSE << TAB << MAE;

	if (EndLine)  foNNres << endl; else foNNres << TAB;
	foNNres.close();
}

void printStats4knnCategory(STRING_TYPE MasterFile, knn &K, matrix<SIGNED_4B_TYPE> &ConfMat, bool EndLine = false)
{
	FILETYPE_OUT foNNres(MasterFile.c_str(), ios_base::app | ios_base::ate);

	foNNres << K.qualV << TAB << K.qsarBLOCK.get_groupN(ConfMat) << TAB;
	
	foNNres << K.qsarBLOCK.get_ccri(ConfMat, -1, true) << TAB;
	foNNres << K.qsarBLOCK.get_ccri(ConfMat) << TAB;
	
	foNNres << K.qsarBLOCK.get_ccr(ConfMat, K.c_wts, 0, 1) << TAB;
	foNNres << K.qsarBLOCK.get_ccr(ConfMat, K.c_wts, 0, 0) << TAB;
	
	foNNres << K.qsarBLOCK.get_ccr(ConfMat, K.c_wts, 0, 3) << TAB;
	foNNres << K.qsarBLOCK.get_ccr(ConfMat, K.c_wts, 0, 2) << TAB;
	
	foNNres << K.qsarBLOCK.get_ccr(ConfMat, K.c_wts, 0, 7) << TAB;
	foNNres << K.qsarBLOCK.get_ccr(ConfMat, K.c_wts, 0, 6);

	for (SIGNED_4B_TYPE ci = 0; ci < K.c_wts.length(); ci++)
		foNNres  << TAB << K.qsarBLOCK.get_groupN(ConfMat, ci) << TAB << K.qsarBLOCK.get_ccri(ConfMat, ci);

	if (EndLine)  foNNres << endl; else foNNres << TAB;
	foNNres.close();
}

bool knnFilter(STRING_TYPE stTrain)
{
	REALNUM_TYPE q2index = knnCore.qualV, r2index = -1;
	
	STRING_TYPE stM = EvalP1;
	
	sprintf(sbuffr, "= %6.4f", q2index);
	stM += sbuffr; 
	stM.parse_string();
	PutInLogFile(stM); //report current's model q2/CCR

	if (q2index < NN_EVAL_Q2) return false;
	
	//check external prediction, it must have >50% coverage (at least 3 datapoints) unless "external-Q2" is used
	apvector<REALNUM_TYPE> predTST;
	SIGNED_4B_TYPE snPred = knnCore.evaluate_ext(datasetTST, predTST);
	
	//mdf. 10.6.2011
	if	( snPred < MIN_GROUP_N) //( 2*snPred < datasetTST.get_Ndatapoints() ) 
	if ( ((NN_EVAL_MODE & KNN_EVAL_EXT_Q2F) == 0) || (snPred == 0) )
	{//if not external-Q2 (which can predict even 1 datapoint) then quit
		PutInLogFile("Test set coverage is too low");
		return false;
	}

	if ((NN_EVAL_MODE & KNN_EVAL_EXT_R2) == KNN_EVAL_EXT_R2)
	{//calculate R2 for continuous kNN
		if (knnCore.qualV < 0)
			r2index = -sqr(knnCore.qualV);
		else
			r2index = sqr(knnCore.qualV);
	}
	else
		r2index = knnCore.qualV;	//otherwise, use quality parameter directly

	knnCore.qualV = q2index; //restore training-set based quality index

	stM = EvalP2;
	sprintf(sbuffr, "(test)= %6.4f", r2index );
	stM += sbuffr;
	stM.parse_string();
	PutInLogFile(stM); //report R2/CCR

	if (r2index < NN_EVAL_R2)	return false;
	
	//simple check is complete, prepare for reporting/extensive check (if needed)
	apvector<REALNUM_TYPE> expTST, expTRN;
	datasetTST.get_ActValues(expTST);
	knnCore.dbase->get_ActValues(expTRN);

	apvector<REALNUM_TYPE> p__TST(predTST), e__TST(expTST);
	knnCore.qsarBLOCK.remove_unpredicted(e__TST, p__TST, knnCore.dbase->get_MaxAct());
	
	//vars for discrete-scale kNN
	SIGNED_4B_TYPE ngroups = 0; 
	apvector<SIGNED_4B_TYPE> pn__, en__;
	matrix<SIGNED_4B_TYPE> cmtx_trn, cmtx_tst;
	
	if (NN_PRED_MODE)
	{//get confusion matrices for test and training sets
		ngroups = NN_EVAL_WTS.length();
		cmtx_tst.SetSize(ngroups, ngroups);
		knnCore.classify_predictions(e__TST, en__);
		knnCore.classify_predictions(p__TST, pn__, (NN_PRED_MODE & KNN_PRED_CLASS) );
		knnCore.qsarBLOCK.get_conf_mtx(en__, pn__, cmtx_tst);
		
		for (SIGNED_4B_TYPE ni = 0; ni < ngroups; ni++)
		if (knnCore.qsarBLOCK.get_groupN(cmtx_tst, ni) < MIN_GROUP_N)
		{
			sprintf(sbuffr, "'Act.group = %d' has < %d datapoints in test set.", ni, MIN_GROUP_N);
			PutInLogFile(sbuffr);
			return false;
		}

		cmtx_trn.SetSize(ngroups, ngroups);
		knnCore.classify_predictions(expTRN, en__);
		knnCore.classify_predictions(knnCore.pred_data, pn__, (NN_PRED_MODE & KNN_PRED_CLASS) );
		knnCore.qsarBLOCK.get_conf_mtx(en__, pn__, cmtx_trn);
	}

	if (NN_FILTER)
	{//extensive check
		if (NN_PRED_MODE)
		{//ccr-based filtering			
			if (NN_FILTER_CCR.length() == ngroups)
			{
				SIGNED_4B_TYPE i;
				//check the test set first
				for (i = 0; i < ngroups; i++)
				if (knnCore.qsarBLOCK.get_ccri(cmtx_tst, i) < NN_FILTER_CCR[i])	return false;
				
				for (i = 0; i < ngroups; i++)
				if (knnCore.qsarBLOCK.get_ccri(cmtx_trn, i) < NN_FILTER_CCR[i])	return false;
			}
		}
		else
		{//continuous
			//R, R2, slope based filtering
			REALNUM_TYPE kR2, k1R2;
			kR2 = knnCore.qsarBLOCK.sqrR0(e__TST, p__TST);
			k1R2 = knnCore.qsarBLOCK.sqrR0(p__TST, e__TST);
			if ( fabs(kR2 - k1R2) > 0.3) return false;
			REALNUM_TYPE k, k1;
			k = knnCore.qsarBLOCK.trendline0(e__TST, p__TST);
			k1 = knnCore.qsarBLOCK.trendline0(p__TST, e__TST);

			//if ( ( (k < 0.85) || (k > 1.15) ) && ( (k1 < 0.85) || (k1 > 1.15) ) ) //simplified filter
			if (((k < 0.85) || (k > 1.15) || (fabs(1 - kR2/r2index) > 0.1)) &&
				((k1 < 0.85) || (k1 > 1.15) || (fabs(1 - k1R2/r2index) > 0.1)))			
				return false;			
		}
	}//if (NN_FILTER)

	//report
	if (OldFormatOutput)	//save individual file
		OutputkNNResults(stOutputFile, stTrain, knnCore, nKnnModelId+1, true, true);
	
	OutputkNNResults(stOutputFile, stTrain, knnCore, ++nKnnModelId);

	OutputkNNPredictions(stOutputPredFile, *(knnCore.dbase), knnCore.pred_data, knnCore.dbase->get_MaxAct(), nKnnModelId);
	OutputkNNPredictions(stOutputPredFile, datasetTST, predTST, knnCore.dbase->get_MaxAct(), nKnnModelId, false, true);

	//print out statistics 	
	if (NN_PRED_MODE)
	{//for category-kNN!
		printStats4knnCategory(stOutputFile, knnCore, cmtx_trn);
		knnCore.qualV = r2index;
		printStats4knnCategory(stOutputFile, knnCore, cmtx_tst, true);
	}
	else
	{//for continuous-kNN!
		printStats4knnContinuous(stOutputFile, knnCore, expTRN, knnCore.pred_data);
		knnCore.qualV = r2index;
		printStats4knnContinuous(stOutputFile, knnCore, e__TST, p__TST, true);
	}

	knnCore.qualV = q2index; //restore training-set based quality index (e.g. used in stepwise comparison)
	return true;
}

void PrintParametersInLOG()
{
	char bfr[1024];
	FILETYPE_OUT fLOG(LOG_FILENAME.c_str(), ios_base::app | ios_base::ate); //appending
	fLOG << "The randomizer seed used: " << srand_seed << endl;

	fLOG << endl << "---------------kNN settings ---------------" << endl;
	if (NN_DIMSEL_METHOD == NN_FULL_DIMSEARCH)
		fLOG << "Exhaustive combinatorial exploration of dimensions." << endl;

	if (NN_DIMSEL_METHOD == NN_STEPWISE)
		fLOG << "Stepwise selection of dimensions." << endl;
	
	if (NN_DIMSEL_METHOD == NN_STEPWISE1)
		fLOG << "Stepwise reduction of dimensions." << endl;

	fLOG << "Range of dimensions: " <<  NN_DIMS_MIN	<< ".." << NN_DIMS_MAX << endl;

	if (NN_R > 0)		
	{
		sprintf(bfr, "%4.2f", NN_R);
		fLOG << "Distance cut-off R: " <<  bfr << endl;
	}
	else
	{
		fLOG << "range of k: " <<  NN_K_MIN	<< ".." << NN_K_MAX << endl << endl;
		fLOG << "Applicability Domain: max.R";
		if ((NN_AD_MODE & 1) == 0)	fLOG << "^2 = <d^2>"; else fLOG << "= <d>";
		fLOG << " + " << NN_AD_ZCUT << "*st.dev" << endl;
		fLOG << "Prediction is within AD if ";
		if (NN_AD_MODE < 2) 
			fLOG << "at least 1 neighbor within AD"; //lax mode, default
		else
		{
			if (NN_AD_MODE & 2) fLOG << "av.dist to k-neighbors within AD"; //compatibility mode
			if (NN_AD_MODE & 4) fLOG << "at least half of k-neighbors within AD";
			if (NN_AD_MODE & 8) fLOG << "all k-neighbors within AD"; //strict mode
		}
		fLOG << endl;
		// { <d> - av.dist. b/w k-nearest neibors }
	}

	fLOG << "Distance metric is ";
	switch (METRIC_K)
	{		
		case 1: fLOG << "Cosine"; break;
		case 2: fLOG << "Correlation-based"; break;
		case 3: fLOG << "Tanimoto"; break;
		default:
		case 0:	fLOG << "Euclidean/Minkowski"; break;
	};
	fLOG << ", power-coff=" << METRIC_V << endl;
	if (NN_APPROX) fLOG << "NB! Fast&lossy search-mode for nearest neighbors is ON!" << endl;

	fLOG << endl << "Distance based weighting scheme for nearest neighbours: " << endl;
	
	if (NN_WT_MODE & KNN_WT_K_FULLD) fLOG << "check found neighbours in space of all dimensions" << endl;
	if (NN_WT_MODE & (KNN_WT_K_VOTE | KNN_WT_K_ALL)) fLOG << "include all extra neighbours" << endl;
	if (NN_WT_MODE & KNN_WT_K_VOTE) fLOG << "take all but only k votes" << endl;
	
	if ((NN_WT_MODE & KNN_WT_F_NO)  == 0)
		fLOG << "all nearest neighbours have equal weights" << endl;
	else
	{
		if (NN_WT_MODE & KNN_WT_F_EXP)  fLOG << "Exponential";
		if (NN_WT_MODE & KNN_WT_F_MNK) 	fLOG << "Minkowski-kind (default)";
		if (NN_WT_MODE & KNN_WT_F_HYP) 	fLOG << "Hyperbolic";

		fLOG << " distance based weight-func()," << endl;

		if ((NN_WT_MODE & KNN_WT_D1) == 0) 	fLOG << "squared ";
		if (NN_WT_MODE & KNN_WT_ABS_D) 	fLOG << "absolute"; else fLOG << "relative";
		fLOG << " distances are used in the weight-func()," << endl;

		sprintf(bfr, "%4.2f", NN_WT_K);
		fLOG << "Weight-func() power-parameter: " <<  bfr << endl;
	}
	fLOG << endl;

	//prediction
	fLOG << "Activity type prediction: ";
	if (NN_PRED_MODE == 0) fLOG << "Continuous" << endl;
	if (NN_PRED_MODE == KNN_PRED_CATEG) fLOG << "Category" << endl;
	if (NN_PRED_MODE == KNN_PRED_CLASS) fLOG << "Classes" << endl;
	if (NN_PRED_MODE == KNN_PRED_EXTERNAL_F) fLOG << "External-prediction func()" << endl;
	fLOG << endl;
		
	if (RunPredictionMode) 
	{
		fLOG.close();
		return;
	}

	//evaluation
	fLOG << "Evaluation:" << endl;	
	fLOG << "Leave-group-out size: ";
	if (NN_EVAL_LGO_F > 0) 
		fLOG << UNSIGNED_4B_TYPE(100*NN_EVAL_LGO_F) << '%';
	else 
		fLOG << NN_EVAL_LGO << " datapoints";
	fLOG << endl;
	
	fLOG << "Min. " << EvalP1 << " and " << EvalP2 << "(test) are";	
	sprintf(bfr, " %4.2f and %4.2f", NN_EVAL_Q2, NN_EVAL_R2);
	fLOG << bfr << endl;

	if ((NN_EVAL_MODE & KNN_EVAL_EXTERNAL_F) == KNN_EVAL_EXTERNAL_F) 
		fLOG << "External evaluation func()" << endl;

	if ((NN_EVAL_MODE & KNN_EVAL_ALT) == KNN_EVAL_ALT)			fLOG << "Alternative evaluation func() is used" << endl;
	if ((NN_EVAL_MODE & KNN_EVAL_ERR) == KNN_EVAL_ERR)			fLOG << "Error-based evaluation func()" << endl;
	if ((NN_EVAL_MODE & KNN_EVAL_AVERR) == KNN_EVAL_AVERR)		fLOG << "based on average-error instead of max-error" << endl; //should be used only for group kNNs!
	if ((NN_EVAL_MODE & KNN_EVAL_EXT_R2) == KNN_EVAL_EXT_R2)	fLOG << "R2 used for test set." << endl;

	if (NN_FILTER)
	{
		if (NN_PRED_MODE)
		{		
			fLOG << "Additional requirements on CCR for activity-groups:";
			for (SIGNED_4B_TYPE i = 0; i < NN_FILTER_CCR.length(); i++)
			{
				sprintf(bfr, "%4.2f", NN_FILTER_CCR[i]);
				fLOG << BLANK << bfr;
			}
			fLOG << endl;
		}
		else
		{
			fLOG << "Additional requirements on correlation:" << endl;
			fLOG << "For y=k*x and x=k1*y (y - observed, x - predicted):" << endl;
			fLOG << "1) corr. coffs should comply with: |R2(k) - R2(k1)| <= 0.3" << endl;
			fLOG << "2) slope coffs should comply with: 0.85 < (k or k1) < 1.15" << endl;
			fLOG << "3) and 0.1 > | 1 - {R2(k) or R2(k1)} /R2 |" << endl << endl;
		}
	}

	if (NN_PRED_MODE)
	{//extra parameters for discrete activity scale
		
		if (NN_PRED_MODE == KNN_PRED_CLASS)
		{
			sprintf(bfr, "%4.2f", NN_CLASSSEP_ADJ);			
			fLOG << "Weight of the class-separation term in the model training fitness f(): " << bfr << endl;
		}

		sprintf(bfr, "%4.2f", NN_EVAL_PENALTY);
		fLOG << "Penalty for imbalanced predictions of activity-groups: " << bfr << endl;
		SIGNED_4B_TYPE i, n = NN_EVAL_WTS.length();
		if (n)
		{
			fLOG << "Weights for activity-groups:";
			for (i = 0; i < n; i++)
			{
				sprintf(bfr, "%4.2f", NN_EVAL_WTS[i]);
				fLOG << BLANK << bfr;
			}
			fLOG << endl;
		}
		
		n = NN_EVAL_BPS.length();
		if (n && (NN_PRED_MODE == KNN_PRED_CATEG))
		{
			fLOG << "Break-points for activity-categories:";
			for (i = 0; i < n; i++)
			{
				sprintf(bfr, "%4.2f", NN_EVAL_BPS[i]);
				fLOG << BLANK << bfr;
			}
			fLOG << endl;
		}
	}

	fLOG << endl;
	if (NN_DIMSEL_METHOD == NN_AC)
	{
		fLOG << "---------------Ant Colony (AC) parameters ---------------" << endl;
		fLOG << "Population size			: " <<  AC_SIZE	<< endl;
		fLOG << "Max. #cycles			: " <<  MAX_N_GENERATIONS	<< endl;
		fLOG << "#best models to store		: " <<  AC_NBEST	<< endl;
		fLOG << "Pheromone volatility rate	: " <<  AC_P	<< endl;
		fLOG << "Min.& max. pheromones		: " <<  AC_TMIN << " .. " << AC_TMAX	<< endl;
		fLOG << "Pheromone impact (alpha)	: " <<  AC_A	<< "(power-factor)" << endl;
		fLOG << "model size impact (beta)	: " <<  AC_B	<< "(power-factor)" << endl;
		if (AC_PH & 1)	fLOG << "Pheromone upate by iteration best model."	<< endl;
		if (AC_PH & 2)	fLOG << "Pheromone upate by global best model only."	<< endl;
		if (AC_PH & 4)	fLOG << "Extra pheromone upate by global best model."	<< endl;
		if (AC_PH & 8)	fLOG << "ASRank mode (models are ranked by quality)."	<< endl;
		if (AC_PH & 16)	fLOG << "Post-optimization mode is ON."	<< endl;
	}

	if (NN_DIMSEL_METHOD == NN_PS)
	{
		fLOG << "---------------Particle Swarms (PS)parameters ---------------" << endl;
		fLOG << "Population size		: " <<  SWARM_SIZE	<< endl;
		fLOG << "Max. #cycles		: " <<  MAX_N_GENERATIONS	<< endl;
		
		fLOG << "Particle inertia weight	: " <<  PS_W	<< endl;
		fLOG << "Cognitive learning rate	: " <<  PS_C1	<< endl;
		fLOG << "Social learning rate	: " <<  PS_C2	<< endl;
		fLOG << "Max.coord. for Velocity	: " <<  PS_VMAX	 << endl;
		fLOG << "Social learning based on: ";
		if (PS_K) fLOG <<  PS_K	<< " neighbors" << endl; else fLOG << "global best" << endl;
	}

	if (NN_DIMSEL_METHOD == NN_SA)
	{
		fLOG << "---------------Simulated Annealing (SA) parameters ---------------" << endl;
		fLOG << "Dimensions increment: " << NN_DIMS_INC << endl;
		fLOG << "#SA-runs for each set of dimensions: " << SA_NRUNS << endl;
		fLOG << "#best models to store             : " <<  SA_NBEST2STORE	<< endl;
		fLOG << "Max. #mutation trials w/o T-change: " <<  SA_N_TRIALS	<< endl;
		
		sprintf(bfr, "%4.2f%%", 100*SA_MUTATION_RATE_PER_GENE);
		fLOG  << "Mutation probability per gene     : " << bfr << endl;

		fLOG  << "Temperature change from 10^" << SA_T_0 << " till 10^" << SA_T_END << endl;
		sprintf(bfr, "%4.2f", SA_K);
		fLOG  << "K-coefficient of T-change is       : " << bfr << endl;
		fLOG  << "Convergence span of T is 10^" << SA_T_CONV << endl;
	}

	if  ( (NN_DIMSEL_METHOD == NN_GA) || (NN_DIMSEL_METHOD == NN_GA1) )
	{
		fLOG << "---------------Genetic Algorithm parameters ---------------" << endl;

		if  (NN_DIMSEL_METHOD == NN_GA)
			fLOG << "Standard GA: #genes in solution is equal #descriptors, each gene is 0 or 1." << endl;
		else
			fLOG << "Alternative GA: gene values are dim ids, solution size is " <<  NN_DIMS_MIN	<< ".." << NN_DIMS_MAX  << endl;
		
		fLOG << "Solution size penalty base: " << GA_SIZE_ADJ << endl;

		fLOG << "Population size         : " <<  POPULATION_SIZE	<< endl;
		fLOG << "Max. #generations: " <<  MAX_N_GENERATIONS	<< endl;
		sprintf(bfr, "%4.2f%%", 100*CROSSOVER_CHANCE);
		fLOG  << "Crossover rate          : " << bfr << endl;
		sprintf(bfr, "%4.2f%%", 100*MUTATION_RATE);
		fLOG  << "Mutation rate           : " << bfr << endl;
		fLOG  << "Parent selectiond mode  : ";
		
		switch (PARENTSELECTION_MODE)
		{
			case RANK_BASED:	fLOG << "RANK BASED";	break;					
			case ROULETTE:		fLOG << "ROULETTE";		break;
			case TOURNAMENT:	fLOG << "TOURNAMENT";	break;
		};//PARENTSELECTION_MODE

		fLOG << endl << "Crossover mode          : ";
		switch (CROSSOVER_MODE)
		{
			case ONE_POINT:			fLOG << "SINGLE POINT";	break;
			case TWO_POINT:			fLOG << "DOUBLE POINT";	break;
			case TWO_OR_ONE_POINT:	fLOG << "SINGLE/DOUBLE POINT";	break;
			case UNIFORM:			fLOG << "UNIFORM";	break;
		};//CROSSOVER_MODE

		fLOG << endl << "Elitism mode:           : ";
		if (ELITISM_MODE)	fLOG << "ON";	else	fLOG << "OFF";
		fLOG << endl;

		if (ELITISM_MODE || (PARENTSELECTION_MODE == RANK_BASED))
		{
			sprintf(bfr, "%4.2f%%", 100*ELITEPART_SIZE);
			fLOG  << "Elite individuals ratio : " << bfr << endl;
		};

		if (PARENTSELECTION_MODE == TOURNAMENT)
			fLOG  << "Tournament group size   : " << GROUP_SIZE << endl;
	}

	if (MIN_DIVERSION != INVALID)
	if  ( (NN_DIMSEL_METHOD == NN_PS)||(NN_DIMSEL_METHOD == NN_AC)||
		(NN_DIMSEL_METHOD == NN_GA)||(NN_DIMSEL_METHOD == NN_GA1) )	
	{
		sprintf(bfr, "Min.diversion: 10^%3d", SIGNED_4B_TYPE(log10(MIN_DIVERSION)) );
		fLOG  << bfr << endl;
		fLOG  << "Convergence condition: " << STABILITY_RANGE << " generations with <= Min.diversion." << endl;					
	};


	fLOG << endl << endl;

	fLOG.close();
}

void PrintMoreHelp(STRING_TYPE ARG)
{
	if (ARG.find("-evl") == 0)
	{
		cout << "'-EVL=[flags]X@X' sets thresholds for model selection" << endl;
		cout << "by their training <def. X=0.5> and test set <def. X=0.6> performance." << endl;
		cout << "For continuous: TRN@TST are Q2@R2<def.> with correlation slope checking." << endl;
		cout << "For class/category: TRN@TST are both CCR<def.> or '-EVL=TRN@TST_CL1@CL2@..CLn'" << endl;
		cout << "where TRN@TST is overall performance and CLi - limits for individual categories." << endl << endl;

		cout << "Additional one-letter flags can change the performance parameters used." << endl;
		cout << ">>>For continuous models:" << endl;
		cout << "E - error-based index: (abs. instead of sqr.diff. for RSS and variation)" << endl;
		cout << "A - alternative index (normalized by variance of predicted values)" << endl;		
		cout << "F - external-Q2 for test set {=1-PRESS(TST)/RSS(TRN) }" << endl;
		cout << "Q - the same index for training and test sets evaluation" << endl;
		cout << "R - R2 as a check for test set <def.>" << endl;
		cout << ">>>For class/category models:" << endl;
		cout << "A - alternative index (accuracy instead of CCR)" << endl;	
		cout << "E - error-based index: (misclassifications as indicator)" << endl;					
		cout << "V - aver.error based (normalized by mean instead of max.misclassification)" << endl << endl;
		
		cout << "S - skips post-evaluation for continuous/category models," << endl;
		cout << "i.e. no check of correlation slope or CCR limits for individual categories." << endl << endl;
	}

	if (ARG.find("-o=ps") == 0)
	{
		cout << "'-PS@...' - Particle Swarm settings: e.g. '-PS@N=200@D=100@S=5@V=-1'" << endl;
		cout << "'..@N=' - swarm size <def." << SWARM_SIZE << ">; '..@D=' - max.#cycles <def." << MAX_N_GENERATIONS << ">" << endl;
		cout << "'..@S=' - max.#stable cycles to stop; '..@V=' - convergence thershold (10^-x)" << endl;		
		cout << "'..@W=' - inertia <def."<< PS_W << ">; '..@VMAX=' - max.velocity <def." << PS_VMAX << ">" << endl;
		cout << "'..@B=' - #neighbors to socialize, 0 is for entire swarm <def." << PS_K << ">" << endl;
		cout << "'..@C1=,..@C2=' - cognitive and social learning rates  <def. " << PS_C1 << " and " << PS_C2 << ">" << endl << endl;
	}

	if (ARG.find("-o=ac") == 0)
	{
		cout << "'-AC@...' - Ant Colony settings: e.g. '-AC@N=300@D=600@S=20@V=-4@E=3@P=0.1'" << endl;
		cout << "'..@N=' - population size <def." << AC_SIZE << ">; '..@D=' - max.#cycles <def." << MAX_N_GENERATIONS << ">" << endl;
		cout << "'..@I=' - ideal fit; '..@V=' - min. significant fitness difference (10^-x)" << endl;
		cout << "'..@S=' - max.#stable cycles to stop; '..@E=' - #best models to store <def." << AC_NBEST << ">" << endl;
		cout << "'..@P=' - pheromone volatility <def."<< AC_P << ">; '..@A=' - pheromone impact <def." << AC_A << ">" << endl;
		cout << "'..@B=' - model size impact on pheromone update <def." << AC_B << ">" << endl;
		cout << "'..@TMIN=,..@TMAX=' - min and max pheromone levels <def. " << AC_TMIN << " and " << AC_TMAX << ">" << endl << endl;

		cout << "Pheromone update modes <def. by all ants in proportion to model fitness>:" << endl;
		cout << "'..@GLB' - by global best, '..@LOC' - by iter best'" << endl;
		cout << "'..@ADG' - extra update by global best at each iter" << endl;
		cout << "'..@ASR' - update based on model rank; '..'@PST' - post-optimize" << endl << endl;
	}

	if (ARG.find("-o=sa") == 0)
	{
		cout << "'-SA@...' - Simulated Annealing settings: e.g. '-SA@B=3@TE=-2@K=0.6@DT=-3@ET=-5'" << endl;
		cout << "'..@N=' - #SA runs to repeat; '..@D=' - #mutations at each T" << endl;
		cout << "'..@FULL' - to do all mutations at each T; '..@B=' - #best models to store" << endl;
		cout << "'..@T0=x' - start T (10^x); '..@TE=' - final T; '..@DT=' - convergence range of T" << endl;
		cout << "'..@K=' - T decreasing coeff.; '..@M=' - mutation probability per dimension" << endl << endl;
	}
	
	if (ARG.find("-o=ga") == 0)
	{
		cout << "GA: #genes in solution is #descriptors, each gene is 0 or 1." << endl;		
		cout << "GA1: gene values are dim ids, solution size is the allowed range of dimensions."  << endl << endl;

		cout << "'-GA@...' - Genetic Algorithm settings: e.g. '-GA@N=500@D=1000@S=20@V=-4@G=7'" << endl;
		cout << "'..@N=' - population size; '..@D=' - max.#generations; '..@I=' - ideal fit" << endl;
		cout << "'..@S=' - #stable generations to stop; '..@V=' - min signif fitness diff (10^x)" << endl;
		cout << "'..@X=' - crossover rate; '..@M=' - mutation rate" << endl;
		cout << "'..@C=' - crossover mode: 1P, 2P, UN, 12" << endl;
		cout << "(i.e. ONE_POINT, TWO_POINT, UNIFORM, TWO_OR_ONE_POINT modes)" << endl;
		cout << "'..@P=' - parent selection mode: RANK, TOUR, RLTT" << endl;
		cout << "'..@G=' - group size for tournament ('TOUR') selection of parents" << endl;
		cout << "'..@E=' - to retain best solutions; e.g. '@E=0.01' (population portion)" << endl; 
		cout << "or '@E=7' (#solutions). Use '@E=OFF' to disable, default is ON with 0.01" << endl;
		cout << "'..@Z=' - to set a penalty term for the solution size <default is 0.1>" << endl << endl;
	}

	if (ARG.find("-o=xx") == 0)
	{
		cout << "Tries all possible combinations of descriptors." << endl;
		cout << "'Use -SEED=.. to limit the set of descriptors to be used." << endl << endl;
	}

	if ( (ARG.find("-o=++") == 0) && (ARG.find("-o=--") == 0) )
	{
		cout << "Gradually adds/removes descriptors until the best model is found." << endl;
		cout << "'Use -SEED=.. to specify a starting set of descriptors." << endl << endl;
	}
}

void PrintKNNHelp(bool detailed)
{//print help
	cout << endl << "kNN+ V" << Version << " - Integrated family of kNN algorithms. " << VerDate << endl;
	cout << "Usage: 'knn+ inputfile [flags]'" << endl;
	cout << "Allowed input files: " << SPLITFILE << ", " << LISTFILE << " (for modeling)" << endl;
	cout << " and " << SUMMARYFILE << ", " << LISTFILE << " (for prediction)" << endl << endl;

	cout << "Use '-4PRED=testfile' for prediction mode, otherwise it is model building mode." << endl;
	cout << "'-OUT=...' - user-defined output file, otherwise inputname is used" << endl << endl;

	cout << "'-O=...' - descriptor selection methods: 'SA'<def.>,'GA','GA1','AC','PS','++','--','XX'" << endl;
	cout << "'GA','GA1' - genetic algorithms, 'SA' - simulated annealing" << endl;
	cout << "'AC' - ant colony, 'PS' - particle swarm" << endl;
	cout << "'++','--' - stepwise kNN, 'XX' - combinatorial search" << endl;
	if (detailed)
		cout << "For more info type 'knn+ -O=...' with the method of interest." << endl;
	cout << endl;
	
	cout << "'-D=...' - range of dimensions to use, e.g. '-D=5@50' or" << endl;
	cout << "'-D=5@50@2' to scan with steps of 2 <def.>" << endl << endl;

	cout << "'-KR=...' sets knn range (#nearest neighbors) or R(Radial cut-off)" << endl;
	cout << "e.g. '-KR=1@9' <def.>(uses 1-9 neighbors) or '-KR=0.5'" << endl << endl;

	cout << "'-M=...' - model type: 'CNT' continuous <def.>,'CTG' category,'CLS' classes" << endl << endl;
	
	if (detailed)
	{
		cout << "'-WT=...' - choosing and weighting neighbors; e.g. -WT=AE2.0" << endl;
		cout << "'2.0' is power-coefficient of weight-f(); E - exponential weight-f()" << endl;
		cout << "'H' - hyperbolic, 'M' - Minkowski-kind weight-f() <def.>," << endl;
		cout << "'N' - no weight-f(), all neighbors are equal;" << endl;
		cout << "'B' - use absolute distances in weight-f(), 'R' - relative <def.>," << endl;
		cout << "'S' - use direct-distances, 'Q' - squared distances <def.>;" << endl;
		cout << "'A' - include all qualified neighbors, 'K'- strictly k only," << endl;
		cout << "'V' - extra neighbors share votes <def.>," << endl;
		cout << "'F' - check extra neighbors in full-D space" << endl << endl;	
	}

	cout << "'-LGO=' - size or fraction of Leave-Group-Out self-prediction <def. 1>" << endl;
	cout << "'-EVL=' - models filter by train/test results: -EVL=0.5@0.6 <def.>" << endl;
	if (detailed)
	{
		cout << "For continuous kNN it means q2 >0.5 and R2>0.6" << endl;
		cout << "A - alternative fit-index; E - error-based fit-index" << endl;
		cout << "V - aver.error based (only for discrete activity modeling)" << endl;
		cout << "R - use R2 as a check for test set <def.>" << endl;
		cout << "F - use external-Q2 for test set" << endl;
		cout << "Q - use the same index for training and test sets" << endl;
		cout << "S - skips post-evaluation (no check of correlation slope, etc.)." << endl;	
		cout << "For more info type 'knn+ -EVL='." << endl << endl;	
	}

	cout << "'-AD=' - applicability domain: e.g. -AD=0.5, -AD=0.5d1_mxk" << endl;
	cout << "'0.5' is z-cutoff <def.>; d1 - direct-distance based AD <def. is dist^2>" << endl;

	if (detailed)
	{
		cout << "Additional options of AD-checking before making prediction:" << endl;
		cout << "'_avd' - av.dist to k neighbors should be within AD (traditional)" << endl;
		cout << "'_mxk' - all k neighbors should be within AD" << endl;
		cout << "'_avk' - k/2 neighbors within AD, '_mnk' - at least 1 within AD <def.>" << endl << endl;

		cout << "'-Z=zx' - metric f() to use; x - is a power coff. <def. 2.0>" << endl;
		cout << "z: E -Euclidean <def.>, T -Tanimoto, R -Corr., C -Cosine" << endl << endl;
	
		cout << "'-SEED=' to start from given dimensions (e.g. =dscr1@dscr2@... or =filename.mod)" << endl << endl;
		cout << "'-SRND=..': to seed randomizer <default is by time>" << endl;

		cout << "Specific settings for categories/classes kNN:" << endl;
		cout << "'-N=' - sets number of groups and their weights and break-points" << endl;
		cout << "e.g. '-N=3W0.6@0.2@0.2B0.5@1.5' sets 3 groups with 0.5 and 1.5 breakpoints and some weights" << endl;
		cout << "'-PNL=' - penalty factor for unbalanced prediction <def. 0>" << endl;
		cout << "'-SEP=' - weight of class-separation term in model evaluation <def. 0>" << endl << endl;

		cout << "'-LAXDIMS' - predicted file's descriptors can be in any order, but labels must match" << endl;
		cout << "'-2OLD' - save individual models into .mod files (old format)" << endl;
		cout << "'-LOGALL' - report progress of each optimization cycle" << endl;
		cout << "'-DOMODS' - report comparison of individual models" << endl;
		cout << "'-DONEIBS' - report neighbors for each prediction" << endl;
		//cout << "'-FAST=' - turns on fast&lossy search for nearest neighbors" << endl;	
	}
	else 
		cout << "type 'knn+ /?' to see additional settings." << endl;
}

void ProcessArgumentString(STRING_TYPE &S)
{
	REALNUM_TYPE rtX;
	
	S.parse_string();
	STRING_TYPE stX = S;
	S.touppercase();
		
	SIGNED_4B_TYPE intX, intU = S.find("-OUT="); 
	if (intU == 0)
	{//output file
		stOutputFile = stX.substr(intU + 5, stX.length());
		return;
	}

	
	if (S.find("-LAXDIMS") == 0)
	{//output file
		NN_PRED_DIM_STRICT = false;
		return;
	}

	intU = S.find("-SRND=");
	if (intU >= 0)
	{		
		intU = atoi( S.substr(intU + 6, S.length()).c_str() );
		srand_seed = UNSIGNED_4B_TYPE(intU);
		srand(srand_seed);
		return;
	}

	intU = S.find("-4PRED=");
	if (intU == 0)
	{//output file
		stScreenedFile = stX.substr(intU + 7, stX.length());
		stScreenedFile.parse_string();
		RunPredictionMode = true;
		return;
	}

	if (S.find("-2OLD") == 0)
	{//output file
		OldFormatOutput = true;
		return;
	}

	intU = S.find("-Z="); //distance metric to use
	if (intU >= 0)
	{
		if (S[intU + 3] == 'C') METRIC_K = 1;
		if (S[intU + 3] == 'R') METRIC_K = 2;
		if (S[intU + 3] == 'T') METRIC_K = 3;
		if (S.length() > intU + 4)
		{
			rtX = atof( S.substr(intU + 4, S.length()).c_str() );
			if (rtX > 0)	METRIC_V = rtX;
		}
		return;
	}

	intU = S.find("-O="); //optimization method
	if (intU == 0)
	{
		if (S.find("=GA") > 0) NN_DIMSEL_METHOD = NN_GA;
		if (S.find("=GA1") > 0) NN_DIMSEL_METHOD = NN_GA1;
		if (S.find("=SA") > 0) NN_DIMSEL_METHOD = NN_SA;
		if (S.find("=AC") > 0) NN_DIMSEL_METHOD = NN_AC;
		if (S.find("=PS") > 0) NN_DIMSEL_METHOD = NN_PS;
		if (S.find("=++") > 0) NN_DIMSEL_METHOD = NN_STEPWISE;
		if (S.find("=--") > 0) NN_DIMSEL_METHOD = NN_STEPWISE1;
		if (S.find("=XX") > 0) NN_DIMSEL_METHOD = NN_FULL_DIMSEARCH;
		return;
	}

	intU = S.find("-SEED="); //seeded solution to try, e.g. -SEED=MW@FR2@C108 or filename.mod
	if (intU == 0)
	{
		S = stX; //use case-sensitive version
		stX = S.substr(intU + 6, S.length());
		stX += SUBSEP;
		intU = stX.find(SUBSEP);
		intX = 0;
		seedims.resize(stX.length());
		while (intU > 0)
		{
			seedims[intX] = stX.substr(0, intU);
			S = stX.substr(++intU, stX.length());
			stX = S;
			intU = stX.find(SUBSEP);
			intX++;
		}
		seedims.resize(intX);
		return;
	}

	intU = S.find("-D="); //range of dimensions to scan and increment; e.g. -D=5@50@2
	if (intU == 0)
	{
		stX = S.substr(intU + 3, S.length());
		stX += SUBSEP;
		intU = stX.find(SUBSEP);
		intX = atoi( stX.c_str() );
		if (intX > 0) NN_DIMS_MIN = intX;
		
		stX = stX.substr(++intU, stX.length());
		intX = atoi ( stX.c_str() );
		if ( intX >= NN_DIMS_MIN ) NN_DIMS_MAX = intX;
		
		intU = stX.find(SUBSEP);
		intX = atoi ( stX.substr(++intU, stX.length()).c_str() );
		if (intX > 0) NN_DIMS_INC = max(1, min(intX, NN_DIMS_MAX - NN_DIMS_MIN));
		return;
	}

	intU = S.find("-KR="); //sets k or R; e.g. -KR=1@9 -KR=0.5
	if (intU == 0)
	{
		stX = S.substr(intU + 4, S.length());
		if (stX.find('.') >= 0)
		{//R cut-off!
			rtX = atof( stX.c_str() );
			if (rtX > 0)
			{
				NN_R = rtX;
				NN_K_MIN = NN_K_MAX = 0;	//force using radial cut-off
			}
			return;
		}

		stX += SUBSEP;
		intU = stX.find(SUBSEP);
		intX = atoi( stX.c_str() );
		if (intX > 0) NN_K_MIN = intX;
		intX = atoi ( stX.substr(++intU, stX.length()).c_str() );
		if ( intX >= NN_K_MIN ) NN_K_MAX = intX;
		
		return;
	}

	intU = S.find("-M="); //model kind
	if (intU == 0)
	{
		stX = S.substr(intU + 3, S.length());
		if (stX == "CTG") NN_PRED_MODE = KNN_PRED_CATEG;
		if (stX == "CLS") NN_PRED_MODE = KNN_PRED_CLASS;
		if (stX == "EXT") NN_PRED_MODE = KNN_PRED_EXTERNAL_F;
		//default is: "CNT" - KNN_PRED_CONTIN		
	}

	intU = S.find("-WT="); //kNN weighting scheme; e.g. -WT=AE2.0
	if (intU == 0)
	{		
		if (strpbrk(stX.c_str(), "AVFK") != NULL)
		{//remove default!
			NN_WT_MODE |= KNN_WT_K_ONLY; 
			NN_WT_MODE -= KNN_WT_K_ONLY;
		}

		stX = S.substr(intU + 4, S.length());
		stX += BLANK;

		intU = 0;
		while ( isalpha(stX[intU]) )
		{
			if (stX[intU] == 'A')	NN_WT_MODE |= KNN_WT_K_ALL;
			if (stX[intU] == 'V')	NN_WT_MODE |= KNN_WT_K_VOTE;
			if (stX[intU] == 'F')	NN_WT_MODE |= KNN_WT_K_FULLD;
						
			if (strchr("NEHM", stX[intU]) != NULL)
			{//overwrite all bits the code weight-functions to ensure only 1 is selected.
				NN_WT_MODE |= KNN_WT_F_NO; 
				NN_WT_MODE -= KNN_WT_F_NO; 
			}

			if (stX[intU] == 'E')	NN_WT_MODE |= KNN_WT_F_EXP;
			if (stX[intU] == 'H')	NN_WT_MODE |= KNN_WT_F_HYP;
			if (stX[intU] == 'M')	NN_WT_MODE |= KNN_WT_F_MNK; //default

			if (stX[intU] == 'B')	NN_WT_MODE |= KNN_WT_ABS_D;	//absolute distances (B), not relative (R)
			if (stX[intU] == 'S')	NN_WT_MODE |= KNN_WT_D1;	//normal distance (S), not squared (Q)

			stX[intU++] = ' ';
		}

		rtX = atof( stX.c_str() );
		if (rtX > 0)	NN_WT_K = rtX;
		return;
	}

	intU = S.find("-EVL="); //models' evaluation e.g. -EVL=SE0.5@0.6  -EVL=AEV0.5@0.7_0.6@0.6 for q2/R2 and CCR filters
	if (intU == 0)
	{
		stX = S.substr(intU + 5, S.length());
		stX += SUBSEP;
		intU = 0;
		while ( isalpha(stX[intU]) )
		{
			if (stX[intU] == 'A') NN_EVAL_MODE |= KNN_EVAL_ALT;		//alternative equation (q2 & accuracy)
			if (stX[intU] == 'E') NN_EVAL_MODE |= KNN_EVAL_ERR;		//error-based equation
			if (stX[intU] == 'V') NN_EVAL_MODE |= KNN_EVAL_AVERR;	//based on average-error, should be used only for group kNNs!
			if (stX[intU] == 'X') NN_EVAL_MODE |= KNN_EVAL_EXTERNAL_F;
			
			if (stX[intU] == 'R') 
			{//use R2 for test set, may require post-evaluation of slope, etc.
				NN_EVAL_MODE -= (NN_EVAL_MODE & KNN_EVAL_EXT_Q2F);
				NN_EVAL_MODE |= KNN_EVAL_EXT_R2;	
			}

			if (stX[intU] == 'F') 
			{//use extQ2 for test set, can be calculated even for 1 datapoint in the test set!
				NN_EVAL_MODE -= (NN_EVAL_MODE & KNN_EVAL_EXT_R2);
				NN_EVAL_MODE |= KNN_EVAL_EXT_Q2F;
			}

			if (stX[intU] == 'Q') //use the same for train and test set
				NN_EVAL_MODE -= (NN_EVAL_MODE & (KNN_EVAL_EXT_R2 | KNN_EVAL_EXT_Q2F));			

			if (stX[intU] == 'S') NN_FILTER		= 0;				//default filtering mode is 1, 0 is Simple (no extra check)
			stX[intU++] = ' ';
		}

		rtX = atof( stX.c_str() );
		if ( fabs(rtX) > 0)	NN_EVAL_Q2 = rtX;
		intU = stX.find(SUBSEP);
		rtX = atof( stX.substr(++intU, stX.length()).c_str() );
		if ( fabs(rtX) > 0)	NN_EVAL_R2 = rtX;
		
		intU = stX.find('_');
		if (intU > 0)
		{//more options for filtering by CCR(i)
			NN_FILTER_CCR.resize(100);
			intU++;
			stX = stX.substr(intU, stX.length());
			intX = 0;
			do
			{
				rtX = atof( stX.c_str() );
				if (rtX > 0)
				{
					if (intX == NN_FILTER_CCR.length()) NN_FILTER_CCR.resize((intX << 1) + 1);
					NN_FILTER_CCR[intX++] = rtX;
				}
				intU = stX.find(SUBSEP);
				stX = stX.substr(++intU, stX.length());
			} while (intU > 0);
			NN_FILTER_CCR.resize(intX);
		}
		return;
	}

	intU = S.find("-LGO="); //leave group out - size
	if (intU == 0)
	{//also add using fraction: e.g. -LGO=0.1
		stX = S.substr(intU + 5, S.length());
		if (stX.find('.') >=0)
		{
			rtX = atof( stX.c_str() );
			if ( (rtX > 0) && (rtX < 1) ) NN_EVAL_LGO_F = rtX;
		}
		else
		{
			intX = atoi( stX.c_str() );
			if (intX > 0)	NN_EVAL_LGO = intX;
		}
		return;
	}

	intU = S.find("-N="); //e.g. -N=3W0.7@0.2@0.2B0.5@1.5 //3 groups with weights and breakpoints
	if (intU == 0)
	{//#categories
		stX = S.substr(intU + 3, S.length());
		intX = atoi( stX.c_str() );
		if (intX <= 0) return;

		//extra parameters for discrete activity scale
		NN_EVAL_BPS.resize(--intX);
		for (intU = 0; intU < intX; intU++) NN_EVAL_BPS[intU] = intU + 0.5;
		NN_EVAL_WTS.resize(++intX);
		for (intU = 0; intU < intX; intU++) NN_EVAL_WTS[intU] = 1.0/intX;

		intX = S.find('W');
		if (intX++ > 0)
		{
			intU = S.find('B');
			stX = S.substr(intX, ( (intU > intX) ? (intU - intX) : S.length() ) );
			intX = intU = 0;
			while ( (intU < NN_EVAL_WTS.length()) && (intX < stX.length()) )
			{
				rtX = atof( stX.c_str() + intX );
				if (rtX > 0) NN_EVAL_WTS[intU] = rtX;
				intX = stX.find(SUBSEP);
				if (intX > 0)	stX[intX] = ' '; else intX = stX.length();
				intU++;
			}
		}

		intX = S.find('B');
		if (intX++ > 0)
		{
			stX = S.substr(intX, S.length() );
			intX = intU = 0;
			while ( (intU < NN_EVAL_BPS.length()) && (intX < stX.length()) )
			{
				rtX = atof( stX.c_str() + intX );
				if ( (rtX > floor(NN_EVAL_BPS[intU])) && (rtX < ceil(NN_EVAL_BPS[intU])) ) 
					NN_EVAL_BPS[intU] = rtX;
				intX = stX.find(SUBSEP);
				if (intX > 0)	stX[intX] = ' '; else intX = stX.length();
				intU++;
			}
		}

		return;
	}

	intU = S.find("-PNL="); //penalty factor for discrete-activity scale kNN
	if (intU == 0)
	{//penalty

		stX = S.substr(intU + 5, S.length());
		rtX = atof( stX.c_str() );
		NN_EVAL_PENALTY = rtX; //hehe, can be negative, anyone cares to encourage disbalance? : )
		return;
	}

	intU = S.find("-SEP="); //weight for class separation term in the model evaluation
	if (intU == 0)
	{
		stX = S.substr(intU + 5, S.length());		
		rtX = atof(stX.c_str());
		if (rtX >= 0) NN_CLASSSEP_ADJ = rtX;		
		return;
	}
	
	intU = S.find("-AD="); //applicability domain settings, e.g. -AD=0.5 or -AD=0.5d1 or -AD=0.5_avk
	if (intU == 0)
	{
		stX = S.substr(intU + 4, S.length());
		intX = stX.find("D1");
		if (intX > 0)	NN_AD_MODE |= 1; else intX = stX.find("D2");
		
		//extra handling needed, otherwise 0.5D1 is interpreted as 0.5*10^1
		if (intX > 0) stX[intX] = '_';
		rtX = atof( stX.c_str() );
		NN_AD_ZCUT  = rtX;
		//"_mnk" is default
		if (stX.find("_MXK") > 0)	NN_AD_MODE |= 8;
		if (stX.find("_AVK") > 0)	NN_AD_MODE |= 4;
		if (stX.find("_AVD") > 0)	NN_AD_MODE |= 2;
		return;
	}

	intU = S.find("-FAST=");
	if (intU == 0)
	{//mode for fast but lossy search of nearest neighbors
		stX = S.substr(intU + 6, S.length());
		NN_APPROX = atoi( stX.c_str() );
	}

	if (S.find("-LOGALL") == 0) NN_REPORT_TERSE = false;
	if (S.find("-DOMODS") == 0) NN_REPORT_MODSIM = true;
	if (S.find("-DONEIBS") == 0) NN_REPORT_NEIBS = true;
	
	//--------------------------------------------------------------------------------------
	intU = S.find("-SA"); //one line for all SA-specific parameters
	if (intU == 0)
	{//Simulated annealing settings

		intU = S.find("N=");
		if (intU > 0)
		{//number of runs to do
			stX = S.substr(intU + 2, S.length());
			intX = atoi(stX.c_str());
			if (intX > 0) SA_NRUNS = intX;
		}
		
		intU = S.find("B=");
		if (intU > 0)
		{//number of best models to store
			stX = S.substr(intU + 2, S.length());
			intX = atoi(stX.c_str());
			if (intX > 0) SA_NBEST2STORE = intX;
		}

		intU = S.find("D=");
		if (intU > 0)
		{//number of iterations (mutation trials) to do at each temperatue
			stX = S.substr(intU + 2, S.length());
			intX = atoi(stX.c_str());
			if (intX > 0) SA_N_TRIALS = intX;
		}

		intU = S.find("FULL"); //mode to go on with iterations after accepted mutation
		if (intU > 0)	SA_FULL_ITER = true;

		intU = S.find("T0=");
		if (intU > 0)
		{//start temperature
			stX = S.substr(intU + 3, S.length());
			rtX = atof( stX.c_str() );
			SA_T_0 = rtX;
		}
		
		intU = S.find("TE=");
		if (intU > 0)
		{//end temperature
			stX = S.substr(intU + 3, S.length());
			rtX = atof( stX.c_str() );
			if (rtX < 0) SA_T_END = rtX;
		}

		intU = S.find("DT=");
		if (intU > 0)
		{//convergence span for the temperature scale
			stX = S.substr(intU + 3, S.length());
			rtX = atof( stX.c_str() );
			if (rtX < 0) SA_T_CONV = rtX;
		}

	 	intU = S.find("M=");
		if (intU > 0)
		{//probability of mutation per gene of solution
			stX = S.substr(intU + 2, S.length());
			rtX = atof( stX.c_str() );
			if ( (rtX >= 0) && (rtX < 1.0) ) SA_MUTATION_RATE_PER_GENE = rtX;
		}

		intU = S.find("K=");
		if (intU > 0)
		{//K coefficient for simulated annealing
			stX = S.substr(intU + 2, S.length());
			rtX = atof( stX.c_str() );
			if ( (rtX >= 0) && (rtX < 1.0) ) SA_K = rtX;
		}
		return;
	}

	//--------------------------------------------------------------------------------------
	intU = S.find("-GA"); //one line for all GA-specific parameters; -GA&N=500&D=1000&...
	if (intU == 0)
	{//Genetic algorithm settings
		intU = S.find("Z=");
		if (intU > 0)
		{//solution size penalty term
			stX = S.substr(intU + 2, S.length());
			if (stX.find('.') >= 0)
			{
				rtX = atof(stX.c_str());
				if ((rtX >= 0) && (rtX < 1)) GA_SIZE_ADJ = rtX;
			}
		}
		
		intU = S.find("N=");
		if (intU > 0)
		{//population size
			stX = S.substr(intU + 2, S.length());
			intX = atoi(stX.c_str());
			if (intX > 0) POPULATION_SIZE = intX;
		}

		intU = S.find("D=");
		if (intU > 0)
		{//maximum #generations
			stX = S.substr(intU + 2, S.length());
			intX = atoi(stX.c_str());
			if (intX > 0) MAX_N_GENERATIONS = intX;
		}

		intU = S.find("S=");
		if (intU > 0)
		{//stability range (#generations)
			stX = S.substr(intU + 2, S.length());
			intX = atoi(stX.c_str());
			if (intX > 0) STABILITY_RANGE = intX;
		}

		intU = S.find("V=");
		if (intU > 0)
		{//minimum variation
			stX = S.substr(intU + 2, S.length());
			intX = atoi(stX.c_str());
			if (intX) MIN_DIVERSION = pow(10.0, REALNUM_TYPE(intX) );
		}

		intU = S.find("G=");
		if (intU > 0)
		{//group size (for tournament selection of parents)
			stX = S.substr(intU + 2, S.length());
			intX = atoi(stX.c_str());
			if (intX > 0) GROUP_SIZE = intX;
		}

		intU = S.find("E=");
		if (intU > 0)
		{//elitism mode (always to store some best individuals)
			ELITISM_MODE = true;
			stX = S.substr(intU + 2, S.length());
			stX.parse_string();
			
			rtX = atof(stX.c_str());
			if (stX.find('.') < 0) rtX /= POPULATION_SIZE;
			if ( (rtX >= 0) && (rtX < 1.0) ) ELITEPART_SIZE = rtX; 			
			
			stX.tolowercase();
			if ( (stX.find("off") == 0) || (stX.length() == 0) || (ELITEPART_SIZE == 0) ) 
				ELITISM_MODE = false;
		}

		intU = S.find("P=");
		if (intU > 0)
		{//mode for selection of parents
			if ( S.find("P=RANK") > 0) PARENTSELECTION_MODE= RANK_BASED;
			if ( S.find("P=TOUR") > 0) PARENTSELECTION_MODE= TOURNAMENT;
			if ( S.find("P=RLTT") > 0) PARENTSELECTION_MODE= ROULETTE;
		}

		intU = S.find("C=");
		if (intU > 0)
		{//crossover
			if ( S.find("C=1P") > 0) CROSSOVER_MODE= ONE_POINT;
			if ( S.find("C=2P") > 0) CROSSOVER_MODE= TWO_POINT;
			if ( S.find("C=UN") > 0) CROSSOVER_MODE= UNIFORM;
			if ( S.find("C=12") > 0) CROSSOVER_MODE= TWO_OR_ONE_POINT;
		}

		intU = S.find("X=");
		if (intU > 0)
		{//crossover rate
			stX = S.substr(intU + 2, S.length());
			rtX = atof( stX.c_str() );
			if ( (rtX >= 0) && (rtX < 1.0) )	CROSSOVER_CHANCE = rtX;
		}

		intU = S.find("I=");
		if (intU > 0)
		{//ideal fit
			stX = S.substr(intU + 2, S.length());
			rtX = atof( stX.c_str() );
			IDEAL_FIT = rtX;
		}

		intU = S.find("M=");
		if (intU > 0)
		{//mutation rate
			stX = S.substr(intU + 2, S.length());
			rtX = atof( stX.c_str() );
			if ( (rtX >= 0) && (rtX < 1.0) ) MUTATION_RATE = rtX;
		}
	}  

//--------------------------------------------------------------------------------------
	intU = S.find("-AC"); //one line for all AC-specific parameters
	if (intU == 0)
	{//Ant colony optimization algorithm settings
		
		intU = S.find("N=");
		if (intU > 0)
		{//anthill size
			stX = S.substr(intU + 2, S.length());
			intX = atoi(stX.c_str());
			if (intX > 0) AC_SIZE = intX;
		}

		intU = S.find("D=");
		if (intU > 0)
		{//maximum #cycles
			stX = S.substr(intU + 2, S.length());
			intX = atoi(stX.c_str());
			if (intX > 0) MAX_N_GENERATIONS = intX;
		}

		intU = S.find("S=");
		if (intU > 0)
		{//stability range (#generations)
			stX = S.substr(intU + 2, S.length());
			intX = atoi(stX.c_str());
			if (intX > 0) STABILITY_RANGE = intX;
		}

		intU = S.find("V=");
		if (intU > 0)
		{//minimum variation
			stX = S.substr(intU + 2, S.length());
			intX = atoi(stX.c_str());
			if (intX) MIN_DIVERSION = pow(10.0, REALNUM_TYPE(intX) );
		}

		intU = S.find("I=");
		if (intU > 0)
		{//ideal fit
			stX = S.substr(intU + 2, S.length());
			rtX = atof( stX.c_str() );
			IDEAL_FIT = rtX;
		}

		intU = S.find("TMIN=");
		if (intU > 0)
		{
			stX = S.substr(intU + 5, S.length());			
			rtX = atof(stX.c_str());
			if (rtX >= 0) AC_TMIN = rtX;			
		}

		intU = S.find("TMAX=");
		if (intU > 0)
		{
			stX = S.substr(intU + 5, S.length());			
			rtX = atof(stX.c_str());
			if (rtX >= 0) AC_TMAX = rtX;			
		}

		intU = S.find("B=");
		if (intU > 0)
		{
			stX = S.substr(intU + 2, S.length());
			AC_B = atof(stX.c_str());			
		}

		intU = S.find("A=");
		if (intU > 0)
		{
			stX = S.substr(intU + 2, S.length());			
			AC_A = atof(stX.c_str());
		}

		intU = S.find("P=");
		if (intU > 0)
		{//volatiliy rate
			stX = S.substr(intU + 2, S.length());
			if (stX.find('.') >= 0)
			{
				rtX = atof(stX.c_str());
				if ((rtX >= 0) && (rtX < 1.0)) AC_P = rtX;
			}
		}

		intU = S.find("E=");		
		if (intU > 0)
		{//number of best models to store
			stX = S.substr(intU + 2, S.length());
			intX = atoi(stX.c_str());
			if (intX > 0) AC_NBEST = intX;
		}

		if (S.find("LOC") > 0)	AC_PH |= 1;
		if (S.find("GLB") > 0)	AC_PH |= 2;
		if (S.find("ADG") > 0)	AC_PH |= 4;
		if (S.find("ASR") > 0)	AC_PH |= 8;
		if (S.find("PST") > 0)	AC_PH |= 16;
	}//ACO
		
	intU = S.find("-PS"); //one line with all parameters, e.g. '-PS@N=200@D=100@S=5@V=-1'
	if (intU == 0)
	{//particle swarm settings
		
		intU = S.find("N=");
		if (intU > 0)
		{//swarm size
			stX = S.substr(intU + 2, S.length());
			intX = atoi(stX.c_str());
			if (intX > 0) SWARM_SIZE = intX;
		}

		intU = S.find("D=");
		if (intU > 0)
		{//maximum #cycles
			stX = S.substr(intU + 2, S.length());
			intX = atoi(stX.c_str());
			if (intX > 0) MAX_N_GENERATIONS = intX;
		}

		intU = S.find("S=");
		if (intU > 0)
		{//stability range (#generations)
			stX = S.substr(intU + 2, S.length());
			intX = atoi(stX.c_str());
			if (intX > 0) STABILITY_RANGE = intX;
		}

		intU = S.find("V=");
		if (intU > 0)
		{//minimum variation
			stX = S.substr(intU + 2, S.length());
			rtX = atof(stX.c_str());
			MIN_DIVERSION = pow(10.0, rtX );
		}		

		//cognitive and social rates
		intU = S.find("C1=");
		if (intU > 0) PS_C1 = atof(S.substr(intU + 3, S.length()).c_str());		
		intU = S.find("C2=");
		if (intU > 0) PS_C2 = atof(S.substr(intU + 3, S.length()).c_str());

		intU = S.find("W="); //inertia weight, can be negative
		if (intU > 0)
		{
			rtX = atof(S.substr(intU + 2, S.length()).c_str());
			if ((rtX < 1.0) && (rtX > -1.0)) //otherwise no convergence
				PS_W = rtX;
		}

		intU = S.find("VMAX="); //maximum for any coordinate of the velocity vector
		if (intU > 0)	PS_VMAX = atof(S.substr(intU + 5, S.length()).c_str());

		intU = S.find("B="); //#neighbors to socialize
		if (intU > 0)
		{
			stX = S.substr(intU + 2, S.length());
			intX = atoi(stX.c_str());
			if (intX >= 0) PS_K = intX;
		}
	} //PSO
}

void GetSeededModel(set &setTD)
{//read specified seed and returns it as setAvDm
	SIGNED_4B_TYPE si = 0, di, nD, dscr, sdL = seedims.length();
	STRING_TYPE stD;
	FILETYPE_IN fiMOD;	
	setTD.Dump();
	while (si < sdL)
	{		
		seedims[si].parse_string();
		stD = seedims[si];
		seedims[si].tolowercase();		
		if (CheckStrEnding(seedims[si], MODELFILE))
		{//read descriptors from the file
			fiMOD.open(stD.c_str());
			if ( fiMOD.eof() || !fiMOD.fail() )
			{
				fiMOD >> nD;			
				seedims.resize(sdL + nD);
				for (di = 0; di < nD; di++)	fiMOD >> dscr; //skip the numbers
				for (di = 0; di < nD; di++)	fiMOD >> seedims[sdL + di];
				sdL += nD;
			}
			fiMOD.close();			
		}
		else
		{
			dscr = knnCore.dbase->get_dscr_pos(stD);
			if (dscr < 0) 
			{
				char *stend;
				dscr = strtol(stD.c_str(), &stend, 10);
				if (stend[0] != char(0)) dscr = 0;
				dscr--;
				if (knnCore.dbase->get_dscr(dscr).length()) setTD.PutInSet(dscr);
			}
			else setTD.PutInSet(dscr);
		}
		si++;
	}//while si
}

bool LoadDatasetOLD(apvector<STRING_TYPE>&NS, dataset&FL, dataset&TR, dataset&TS)
{//loads a dadtaset from an old format list-file split, e.g. 'aa_trn0.x aa_trn0.a 42 aa_tst0.x aa_tst0.a 6'
	if (NS.length() < 6) return false;
	if (!LoadDataset(TR, NS[0], NS[1])) return false;
	if (!LoadDataset(TS, NS[3], NS[4])) return false;
			
	//reconstruct full dataset
	FL = TR;
	if (!FL.expandby(TS, true)) return false;
	
	set XXX(0, FL.get_Ndatapoints());
	FL.train = XXX - FL.test;
	return true;
}

SIGNED_4B_TYPE GetSplitNum(STRING_TYPE &S)
{//identify split number
	SIGNED_4B_TYPE u = S.length();
	while (u > 0) if ( !isdigit(S[--u]) ) break;
	return atoi( S.substr(++u, S.length()).c_str() );
}

bool LoadFromSplit(STRING_TYPE &Args, STRING_TYPE &oldArgs, dataset &datasetT, dataset &dtBase, dataset &dtTest)
{//loads training/test datasets depending on the global work mode (RunPredictionMode)
	if (oldArgs == Args) return (Args.length() > 0); //skip loading, if dtBase is going to be the same!

	oldArgs = "";	//erase to invalidate (by default)

	SIGNED_4B_TYPE cnt = INVALID;
	apvector<STRING_TYPE> ap_strs;	
	SplitString(Args, BLANK, ap_strs);

	datasetT.dump();
	dtBase.dump();
	if (!RunPredictionMode) dtTest.dump();

	if (ap_strs.length() < 3) return false;
	if (ap_strs.length() == 3) 
	{
		ap_strs.resize(4);
		ap_strs[3] = ap_strs[2];
		ap_strs[2] = ap_strs[1];
		ap_strs[1] = ap_strs[0];
	}

	if (!CheckStrEnding(ap_strs[2], SPLITFILE)) 
	{
		if (CheckStrEnding(ap_strs[2], LISTFILE))
		{
			cnt = GetSplitNum(ap_strs[3]);
			if (cnt > 0)
			{
				FILETYPE_IN fiList(ap_strs[2].c_str());
				if (!(fiList.eof() || fiList.fail()))
				{
					STRING_TYPE stX = SkipFileComments(fiList);
					while (--cnt) stX.getlinewithtabs(fiList);
					SplitString(stX, BLANK, ap_strs);
					dataset dtTmpTest;
					if (LoadDatasetOLD(ap_strs, datasetT, dtBase, dtTmpTest))
					{//success, initialize						
						fiList.close();
						if (!RunPredictionMode) dtTest = dtTmpTest;
						oldArgs = Args;
						return true;
					}
				}
				fiList.close();
			}
		}
		//try to load at least the training set (still useful for prediction)
		if (LoadDataset(dtBase, ap_strs[0], ap_strs[1]))
		{
			oldArgs = Args;		//initialize
			return true;
		}
		return false;
	}

	cnt += GetSplitNum(ap_strs[3]);	
	if (cnt < 0) return false;

	FILETYPE_IN fiSplit(ap_strs[2].c_str());
	if (fiSplit.eof() || fiSplit.fail())
	{
		fiSplit.close();
		return false;
	}

	if (!LoadDataset(datasetT, ap_strs[0], ap_strs[1]) ) return false;
	
	SkipFileComments(fiSplit);

	cnt *= 3; //#lines to skip
	STRING_TYPE stX;
	while (cnt--) stX.getlinewithtabs(fiSplit);

	LoadSetAsText(fiSplit, datasetT.train);
	LoadSetAsText(fiSplit, datasetT.test);

	dtBase = datasetT.get_training_set();
	if (!RunPredictionMode)
		dtTest = datasetT.get_test_set();
	fiSplit.close();

	oldArgs = Args;
	return true;
}

void ReadModelFile(FILETYPE_IN &fiMod, REALNUM_TYPE & rtR, UNSIGNED_2B_TYPE &unK, set &setDm, apvector<STRING_TYPE> &DimLabels)
{
	SIGNED_4B_TYPE uu, zz;
	uu = LoadSetAsText(fiMod, setDm, 0);
	DimLabels.resize(uu);
	for (zz = 0; zz < uu; zz++)	fiMod >> DimLabels[zz];
	unK = 0; //if 0 then R-based cut-off will be used!
	STRING_TYPE stmpx;
	fiMod >> stmpx;
	fiMod >> knnCore.qualV;
	stmpx.parse_string();
	rtR = unK = 0;
	if (stmpx.find('.') > 0) rtR = atof(stmpx.c_str()); else unK = atoi(stmpx.c_str());
	stmpx.getlinewithtabs(fiMod); //to get the rest of the line	
}

bool VerifyDescriptors(apvector<STRING_TYPE> &Lbls, set &Dims, dataset &TRN, dataset &TST, bool ifStrict)
//checks if a descriptor set  is present in TRN and TST datasets:
//Lbls - descriptor names, Dims - their positions
{
	set setX;
	SIGNED_4B_TYPE i, pz, tpz, dn = Lbls.length(), dpn = TST.get_Ndatapoints();
	for (i = 0; i < dn; i++)
	{
		pz = TRN.get_dscr_pos(Lbls[i]);
		if (pz < 0)	return false;		
		
		if (dpn)
		{//if test and modeling set are compared
			tpz = TST.get_dscr_pos(Lbls[i]);
			if (tpz < 0)
			{
				PutInLogFile("Some descriptors are absent in the test set!");
				return false;
			}
			if (ifStrict && ( pz != tpz ))	return false; //position mismatch b/w test and modeling set
		}

		setX.PutInSet(pz);
	}

	set setXX = Dims & setX;
	if (setXX.Size() != dn)
	{//Inconsistency between labels and positions in the modeling set itself
		PutInLogFile("NB! Descriptor labels do not match with original positions..");
		Dims = setX;
	}
	return true;
}

void DescriptorAnalysis(STRING_TYPE In, dataset &B, apvector_set_type &D)
{
	STRING_TYPE stF(In), stSpace, stSim;
	CutStrEnding(stF);
	stSpace = stF + "_modspace.txt";
	stSim = stF + "_modsimil.txt";
	FILETYPE_OUT foSpace(stSpace.c_str()), foSim(stSim.c_str());

	set setX;
	REALNUM_TYPE rtX;
	SIGNED_4B_TYPE i, j, nm = D.length(), nd = B.get_Ndimensions();
	apvector<SIGNED_4B_TYPE> apDims;
	matrix<REALNUM_TYPE> Msim, Mspace;
	
	Msim.SetSize(nm, nm);
	Mspace.SetSize(nm, nd);
	Mspace.Null();
	
	//space of models matrix (based on descriptors)
	foSpace << "Model#";	
	for (i= 0; i < nd; i++) foSpace  << TAB << B.get_dscr(i);
	foSpace << endl;

	//model similarity matrix
	foSim << "Model#";
	for (i= 0; i < nm; i++) foSim  << TAB << (i+1);
	foSim << endl;
		
	for (i = 0; i < nm; i++)
	{
		Msim(i,i) = 1;
		for (j = i + 1; j < nm; j++)
		{
			setX = D[i] & D[j];
			rtX = setX.Size();
			if (rtX > 0)	
			{
				setX = D[i] | D[j];
				rtX /= setX.Size();
			}
			Msim(i, j) = Msim(j, i) = rtX;
		}

		foSim << (i + 1);
		for (j= 0; j < nm; j++) foSim  << TAB << Msim(i, j);
		foSim << endl;

		D[i].GetList(apDims);
		for (j = 0; j < apDims.length(); j++) Mspace(i,apDims[j]) = 1;

		foSpace << (i + 1);
		for (j= 0; j < nd; j++) foSpace  << TAB << Mspace(i, j);
		foSpace << endl;
	} //for i

	foSpace.close();
	foSim.close();
}
