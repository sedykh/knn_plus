/*
kNN - Class-library for the family of nearest-neighbour algorithms
May. 2008 - Jan. 2011

Evaluation of model's quality: evaluate():
penalty and weights are applied in 2-tier way for the same purpose: 
to compensate for imbalanced datasets. 
Since weights are not optimized, overall penalty can do the final adjustment.

Coordinators: Aleksandr Sedykh, Alexander Golbraikh, Christopher Grulke
*/

#if !defined(KNN_CLASSES)
#define KNN_CLASSES

#include "dataset.h"

//definitions for kNN settings

//weighting scheme and dealing with extra # nearest neighbours
//default is to use k nearest neighbours only!
#define KNN_WT_K_ALL			1	//use all satisfactory nearest neighbours found
#define KNN_WT_K_VOTE			2	//use all but the votes of extra neighbours will be shared over the remaining seats 
#define KNN_WT_K_FULLD			4	//use full-dimensions to see which neighbors are closest
#define KNN_WT_K_ONLY			7	//mask for KNN_WT_K_FULLD | KNN_WT_K_VOTE | KNN_WT_K_ALL

//parameters and equations to calculate the weights of neighbors in the k-neighborhood: 
#define KNN_WT_D1				8	//conventional distances, if OFF then squared distances are used (similar to 0bit in ADmode-variable in knn class)
#define KNN_WT_ABS_D			16	//absolute distances (i.e. X = Di or Di^2), if OFF - relative distances are used (i.e. X = Di/Sum(Dk) or Di^2/Sum(Dk^2))
#define KNN_WT_F_HYP			32	//Hyperbolic equation:		Weight = (1 + X)^-wt_k; (inverse of distance); hyperbolic for wt_k = 1
#define KNN_WT_F_MNK			64	//Minkowski\Ferma-kind:		Weight = (1 - X^wt_k)^1/wt_k; linear for wt_k = 1; only for relative distances (KNN_WT_ABS_D bit is ignored);
#define KNN_WT_F_EXP			128	//Exponential (Gaussian):	Weight = exp(-wt_k*X^2); only for conventional distances (i.e. KNN_WT_D1 bit is ignored)
#define KNN_WT_F_NO				224 //mask for KNN_WT_F_HYP | KNN_WT_F_MNK | KNN_WT_F_EXP
//NB: default is equal weights for all neighbors, if no equation is specified!!! 
//NB: one adjustable parameter is wt_k, wt_k should be > 0; 

//prediction mode
#define KNN_PRED_CONTIN			0	//continuous scale
#define KNN_PRED_CATEG			1	//descrete scale (categories)
#define KNN_PRED_CLASS			2	//multi-class
#define KNN_PRED_EXTERNAL_F		4	//external prediction function

//evaluation mode
#define KNN_EVAL_DEF			0	//q2 for continous, Accuracy for group KNNs
#define KNN_EVAL_ALT			1	//alternative equations for q2 and Accuracy
#define KNN_EVAL_ERR			2	//MAEq for continous, Error-based for group kNNs
#define KNN_EVAL_AVERR			4	//only for group kNNs: calc.Error using av.error / max.error
#define KNN_EVAL_CLASS_ERR		8	//to adjust error-based functions for classes
#define KNN_EVAL_EXT_R2			16	//for continous only: always R2 for the test set
#define KNN_EVAL_EXT_Q2F		32	//for continous only: always external-Q2 for the test set (Consonni et al 2009)
#define KNN_EVAL_EXTERNAL_F		256	//external evaluation function

typedef REALNUM_TYPE (*kNNExternPredFunc)(apvector<SIGNED_4B_TYPE> &, apvector<REALNUM_TYPE> &);
typedef REALNUM_TYPE (*kNNExternEvalFunc)(apvector<REALNUM_TYPE> &, apvector<REALNUM_TYPE> &);

class knn
{
public:
//initialization
	knn();
	knn(const knn &knn_in);
	virtual ~knn();
	void dump();

	void configure(set *setDims = NULL, dataset *dsDbase = NULL,
	UNSIGNED_2B_TYPE nK = 1, REALNUM_TYPE rtR = 0.0, 
	UNSIGNED_2B_TYPE nWt_mode = KNN_WT_F_MNK, REALNUM_TYPE rtWt_k = 1.0,
	//prediction
	UNSIGNED_2B_TYPE nPr_mode = 0, GENERIC_POINTER pPr_f = NULL,
	//evaluation
	UNSIGNED_2B_TYPE nEvl_mode = 0, UNSIGNED_2B_TYPE nEv_lgo = 1,  
	GENERIC_POINTER pEv_f = NULL,
	GENERIC_POINTER rtEv_pred_bps = NULL, GENERIC_POINTER rtEv_pred_wts = NULL, 
	REALNUM_TYPE rtEv_pred_pn = 0, REALNUM_TYPE rtCl_sep_cf = 0,
	UNSIGNED_1B_TYPE nMetrK = 0, REALNUM_TYPE rtMetrV = 2.0,
	REALNUM_TYPE rtZcutoff = 0.5, UNSIGNED_1B_TYPE nADmode = 0, 
	UNSIGNED_1B_TYPE nAprx_mode = 0);

	//to pick from b{} Nearest Neighbours for a datapoint el
	UNSIGNED_1B_TYPE pick(SIGNED_4B_TYPE el, set &b, REALNUM_TYPE ADCutOff = INVALID);
	void calc_dist();	//recalculates kdist matrix for currently selected dimensions
	REALNUM_TYPE calc_AD();
	
	void calc_class_sim();	//recalculates interclass similarities for class-based kNN mode
	REALNUM_TYPE calc_mean_class_separation();

	//activity prediction
	void predict(set &b, set &t, REALNUM_TYPE AD = INVALID);
	REALNUM_TYPE predict_ext(apvector<REALNUM_TYPE> &dp, REALNUM_TYPE AD = INVALID);

	void classify_predictions(apvector<REALNUM_TYPE> &p_data, apvector<SIGNED_4B_TYPE> &p_ndata, bool ifSimple = true); //service function for categorial kNNs
	//quality evaluation
	REALNUM_TYPE compare(apvector<REALNUM_TYPE> &e_data, apvector<REALNUM_TYPE> &p_data, bool extrnVal = false);
	void evaluate();
	SIGNED_4B_TYPE evaluate_ext(dataset &ext_ds, apvector<REALNUM_TYPE> &pr_act, bool calcqual = true, bool bylabel = false);
	SIGNED_4B_TYPE evaluate_ext(dataset &ext_ds, apvector<REALNUM_TYPE> &pr_act, apvector<SIGNED_4B_TYPE> &pr_neibs, bool calcqual = true, bool bylabel = false);

//analysis tools
	QSAR qsarBLOCK;				
//data & nearest neighbour settings
	dataset * dbase;	//modeling dataset, should have train & test sets initialized
	set dims;			//dimensions to use

	UNSIGNED_1B_TYPE aprx_mod;	//defines with neighbours-search is done by lossy approximation!
	//apvector<lneib> dimneib, dimneib1;	//to store potential neighbors for each dimension

	UNSIGNED_1B_TYPE metrK;		//0 - Euclidean metric
	REALNUM_TYPE metrV;			//2.0 - metric coefficient (power for Euclidean metric)
	//AD - Applicability Domain: (R^2)max = <d^2> + zs; <..> is aver., s is st.dev., z is z-cutoff
	//NB: for RNN AD is not needed since the actual specified R is used
	REALNUM_TYPE ADzCutoff;		//0.5 - z-cutoff 
	UNSIGNED_1B_TYPE ADmode;	//Applicability domain options
								//0bit: 0 - uses squared distances to calculate AD, 1 - uses direct distances
								//1bit: uses "k-average distance to compare with AD" (compatibility option)
								//2bit: uses "#neighbors within AD should >= k/2" k-average mode
								//3bit: uses "#neighbors within AD should >= k" strict mode
								//default is "at least 1 neighbor within AD" lax mode

	UNSIGNED_2B_TYPE k;			//default is 1; if 0 then radial cut-off is used
	REALNUM_TYPE r;				//radial cut-off

	UNSIGNED_2B_TYPE wt_mode;	//weighting mode, also describes what to do if there are extra neighbours
	REALNUM_TYPE wt_k;			//weighting parameter (i.e. smoothing factor K, etc.)

//prediction
	UNSIGNED_2B_TYPE pred_mode;				//prediction mode: continuous, discrete, multi-class, external function
	GENERIC_POINTER pred_f;					//arbitrary prediction function	

//self-prediction, modeling training
	matrix<REALNUM_TYPE> kdist;			//distance matrix for the current subspace
	
	matrix<REALNUM_TYPE> iclass_sim;	//interclass similarities for class-based kNN
	REALNUM_TYPE class_sep_cf;			//weight of the class-separation term in the model fitness evaluation

	UNSIGNED_2B_TYPE lgo;				//size of leave-group-out (cross-validation), default 1
	lneib lgo_sets;						//lgo sets, if lgo > 1
	apvector<REALNUM_TYPE> pred_data;	//predicted values for datapoints in dataset

//evaluation 
	UNSIGNED_2B_TYPE eval_mode;	//evaluation mode: various functions to use
	GENERIC_POINTER eval_f;		//arbitrary evaluation function

//variables for class or category scale only:
	apvector<REALNUM_TYPE> c_bps, c_wts; //break points and weights
	REALNUM_TYPE pn;			//penalty for unbalanced predictions (when one group is predicted much better)

//model's quality descriptors
	REALNUM_TYPE qualV; //accuracy and ability to give prediction	

//model's work variables to store nearest neighbors and their distances (for 1 current prediction only!)
	apvector<SIGNED_4B_TYPE> apvNeibs;
	apvector<REALNUM_TYPE> apvDists;
};

#endif // !defined(KNN_CLASSES)
