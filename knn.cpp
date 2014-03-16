/*
Implementation of nearest neighbour algorithms.

NB:
1) min & max values for the continuous range of activity 
should be always added as starting points in sphere exclusion for kNN!
2) it is faster to do LGO cross-validation in kNN instead of LOO
3) kd-tree algorithm can help in cases where #datapoints >> #dimensions!!

TODO: (Add here problems and suggestions to be worked upon)
1)  Implement lossy but fast versions of finding nearest neighbours in predict()
	dimneib structures (nearest neighbours by individual dimensions)
	is only helpful for Euclid-like distances, not those based on similarity etc.! 

DONE:
Oct-Nov 2011	class-kNN implemented: calc_class_sim(), calc_mean_class_separation(), etc.

Oct 2011		ext-Q2 added for external set evaluation (QSAR::q2F1F3())

January 2011	weight-f() overhauled, now 4 basic options (no weights, hyperbolic, exp, minkowski-kind)

April 2010		fix for reporting neibors, a version of evaluate_ext() now can report pr_neibs()
				in a format of {#id #neibor #id #neibor #id #neibor}, where #id can repeat.

Jan 2010		special mode introduced: KNN_EVAL_EXT_R2, to evaluate
				model quality for the test set (as part of model building);

				minor fix in the call of qsarBLOCK.get_ccr() for class kNN

Dec 20 2009		expansion of evaluate_ext(), descriptor label check-mode is added
				apvCands[], apvDists[] made as knn public members (for reporting purposes)
				consequently, pick() input was adjusted, etc.

June 30 2009	small corrections
June 22 2009	initial version of predict()
May 27 2009		initial implementation
Dec 16 2008		knn:evaluate()
*/

#include "knn.h"

knn::knn()
{
	dump();
}

knn::knn(const knn &knn_in)
: dims(knn_in.dims), kdist(knn_in.kdist), iclass_sim(knn_in.iclass_sim), lgo_sets(knn_in.lgo_sets), pred_data(knn_in.pred_data), c_bps(knn_in.c_bps), c_wts(knn_in.c_wts), apvNeibs(knn_in.apvNeibs), apvDists(knn_in.apvDists) //, dimneib(knn_in.dimneib), dimneib1(knn_in.dimneib1)
{
	dbase = knn_in.dbase;
	k = knn_in.k;
	r = knn_in.r;
	
	ADmode = knn_in.ADmode;
	ADzCutoff = knn_in.ADzCutoff;
	
	metrK = knn_in.metrK;
	metrV = knn_in.metrV;
	
	aprx_mod = knn_in.aprx_mod;

	wt_mode = knn_in.wt_mode;
	wt_k = knn_in.wt_k;

	//prediction
	pred_mode = knn_in.pred_mode;
	pred_f = knn_in.pred_f;

	//evaluation
	eval_mode = knn_in.eval_mode;
	lgo = knn_in.lgo;
	
	eval_f = knn_in.eval_f;
	
	pn = knn_in.pn;

	class_sep_cf = knn_in.class_sep_cf;

	//etc
	qualV = knn_in.qualV;
}

knn::~knn()
{
	dump();
}

void knn::dump()
{
	configure();
}

void knn::configure(set *setDims, dataset *dsDbase, 
UNSIGNED_2B_TYPE nK, REALNUM_TYPE rtR,
UNSIGNED_2B_TYPE nWt_mode, REALNUM_TYPE rtWt_k,
//prediction
UNSIGNED_2B_TYPE nPr_mode, GENERIC_POINTER pPr_f,
//evaluation
UNSIGNED_2B_TYPE nEvl_mode, UNSIGNED_2B_TYPE nEv_lgo,  
GENERIC_POINTER pEv_f, 
GENERIC_POINTER rtEv_pred_bps, GENERIC_POINTER rtEv_pred_wts, REALNUM_TYPE rtEv_pred_pn,  REALNUM_TYPE rtCl_sep_cf, 
UNSIGNED_1B_TYPE nMetrK, REALNUM_TYPE rtMetrV,
REALNUM_TYPE rtZcutoff, UNSIGNED_1B_TYPE nADmode, 
UNSIGNED_1B_TYPE nAprx_mode)
//NB: kdist and iclass_sim matrices are set to 0 dimensions at every configure() call!
{
	apvNeibs.resize(0); 
	apvDists.resize(0);

	iclass_sim.SetSize(0,0);
	kdist.SetSize(0,0);
	dbase = dsDbase;
	if (dbase == NULL)
	{
		pred_data.resize(0);
//		dimneib.resize(0);
//		dimneib1.resize(0);
	}
	else
	{
		SIGNED_4B_TYPE np = dbase->get_Ndatapoints();

		pred_data.resize(np);
		REALNUM_TYPE x = 2*dbase->get_MaxAct() - dbase->get_MinAct(); //to use as a dummy value
		while (np) pred_data[--np] = x;
		
		//dbase->calc_allDneibs(k + lgo - 1, dimneib, dimneib1);
	}

	if (setDims == NULL) dims.Dump(); else	dims = (*setDims);

	k = nK;
	r = rtR;
	
	//Applicability Domain, distance-metrics, etc.
	ADmode = nADmode;
	ADzCutoff = rtZcutoff;
	
	metrK = nMetrK;
	metrV = rtMetrV;
	aprx_mod = nAprx_mode;
	

	wt_mode = nWt_mode;
	wt_k = rtWt_k;
	if ((wt_mode & KNN_WT_F_MNK) && (wt_mode & KNN_WT_ABS_D)) wt_mode -= KNN_WT_ABS_D;
	if (wt_mode & KNN_WT_F_EXP) wt_mode |= KNN_WT_D1;

//config prediction
	pred_mode = nPr_mode;

	if (pred_mode & KNN_PRED_EXTERNAL_F)
		pred_f = pPr_f;
	else
		pred_f = NULL;

//config evaluation
	eval_mode = nEvl_mode;
	lgo = nEv_lgo;

	if ( (lgo > 1) && (dbase != NULL) )
		dbase->lgo_split(lgo_sets, lgo); //should be set only once for the dataset!

	if (eval_mode & KNN_EVAL_EXTERNAL_F)
		eval_f = pEv_f;
	else
		eval_f = NULL;

	qualV = pn = class_sep_cf = 0;
	c_bps.resize(0);
	c_wts.resize(0);
	if ( (pred_mode & KNN_PRED_CLASS) || (pred_mode & KNN_PRED_CATEG) )
	{	
		apvector<REALNUM_TYPE> * pApvS = NULL;
		pApvS = (apvector<REALNUM_TYPE> *)rtEv_pred_bps;
		if (pApvS == NULL)
		{
			c_bps.resize(1); 
			c_bps[0] = 0.5;
		}
		else
			c_bps = (*pApvS);

		pApvS = (apvector<REALNUM_TYPE> *)rtEv_pred_wts;
		SIGNED_4B_TYPE nc = c_bps.length() + 1;
		if (pApvS) c_wts = (*pApvS);
		if ( c_wts.length() != nc )
		{
			c_wts.resize(nc);
			while (nc--) c_wts[nc] = 1;
		}
		for (pn = nc = 0; nc < c_wts.length(); nc++) 
		{
			c_wts[nc] = fabs(c_wts[nc]);
			pn += c_wts[nc];
		}

		if (pn > SMALL)
			for (nc = 0; nc < c_wts.length(); nc++) c_wts[nc] /= pn;
		else
			for (nc = 0; nc < c_wts.length(); nc++) c_wts[nc] = 1 / c_wts.length();

		pn = rtEv_pred_pn;
		class_sep_cf = rtCl_sep_cf;
	}	
}

void knn::calc_class_sim()
//precondition: kdist() should be initialized at least with INVALID distances so that pick() can populate it
{	
	if ((pred_mode & KNN_PRED_CLASS) == 0) return; //quit if not in a class mode	
	
	SIGNED_4B_TYPE N = c_wts.length();
	if (iclass_sim.ColNO() == UNSIGNED_2B_TYPE(N)) return; //no change, no need to recalculate
	
	if (dbase == NULL) return;
	SIGNED_4B_TYPE i, j, np = dbase->get_Ndatapoints();

	iclass_sim.SetSize(N, N);
	apvector<set> class_members(N);		

	for (i = 0; i < np; i++) class_members[ SIGNED_4B_TYPE(dbase->get_Act(i)) ].PutInSet(i);

	//reset neighbor-searching settings for the interclass purpose
	UNSIGNED_2B_TYPE old_k = k, old_wt_mode = wt_mode;
	k = 1; 
	wt_mode -= (wt_mode & KNN_WT_K_ONLY);

	//fill up nearest neighbor intra and inter-class distances
	matrix<REALNUM_TYPE> near_dists;	
	near_dists.SetSize(np, N);
	for (i = 0; i < np; i++)
	for (j = 0; j < N; j++)	
	{
		bool sameclass = (SIGNED_4B_TYPE(dbase->get_Act(i)) == j);
		if (sameclass) class_members[j].RemoveFromSet(i);
		pick(i, class_members[j]);
		near_dists(i,j) = apvDists[0];
		if (sameclass) class_members[j].PutInSet(i);
	}

	//restore valid settings
	k = old_k;
	wt_mode = old_wt_mode;

	//calculate class-by-class similarities
	apvector<REALNUM_TYPE> cl_dists;
	apvector<SIGNED_4B_TYPE> cl_lists;
	REALNUM_TYPE av, std;
	for (i = 0; i < N; i++)
	{//first, calculate AD(k=1) distribution for each class
		iclass_sim(i, i) = 1.0;
		class_members[i].GetList(cl_lists);	
		cl_dists.resize(cl_lists.length());
		for (j = 0; j < cl_lists.length(); j++)	
			cl_dists[j] = near_dists(cl_lists[j], i);		
		av = qsarBLOCK.meanV(cl_dists);
		std = qsarBLOCK.stdev(cl_dists);
		
		for (j = 0; j < N; j++)
		{//based on i-th class AD, get individual similarities of j-th class compounds to it, then average
			if (i == j) continue;
			class_members[j].GetList(cl_lists);	
			cl_dists.resize(cl_lists.length());
			for (SIGNED_4B_TYPE p = 0; p < cl_lists.length(); p++)	
				cl_dists[p] = qsarBLOCK.norm_cdf(near_dists(cl_lists[p], i), av, std);
			
			iclass_sim(j,i) = qsarBLOCK.meanV(cl_dists);
		}
	}
}


REALNUM_TYPE knn::calc_mean_class_separation()
{
	SIGNED_4B_TYPE N = c_wts.length();
	if (iclass_sim.ColNO() != UNSIGNED_2B_TYPE(N)) return 0; //not valid, skip	
	SIGNED_4B_TYPE nV = N*(N-1);
	apvector<REALNUM_TYPE> simv(nV,0);
	
	for (SIGNED_4B_TYPE i = 0; i < N; i++)
	for (SIGNED_4B_TYPE j = 0; j < N; j++)
	{
		if (i == j) continue;
		simv[--nV] = iclass_sim(i,j);
	}

	return (1 - qsarBLOCK.meanV(simv));
}

void knn::calc_dist()
{
	if (dbase == NULL) return;
	SIGNED_4B_TYPE i, j, N = dbase->get_Ndatapoints();
	kdist.SetSize(N, N);
	for (i = 0; i < N; i++)
	{
		kdist(i, i) = 0;
		for (j = i+1; j < N; j++)
		kdist(i, j) = kdist(j, i) = dbase->get_indDistance(i, j, &dims, metrV, metrK); //metric etc
	}	
}

REALNUM_TYPE knn::calc_AD()
{
	UNSIGNED_1B_TYPE kmode = 0;	
	bool ADbySquareD = ((ADmode & 1) == 0);
	if (ADbySquareD) kmode++;
	if ( (wt_mode & KNN_WT_K_VOTE) || ((wt_mode & KNN_WT_K_ALL) == 0) ) kmode += 2; //use strictly k only

	SIGNED_4B_TYPE fdbn = dbase->get_Ndatapoints();
	if ( SIGNED_4B_TYPE(kdist.ColNO()) != fdbn )
		calc_dist();
	else
	if (lgo > 1) //fill kdist values for the members of the same leave-out-group
	{
		for (SIGNED_4B_TYPE lli = 0; lli < fdbn-1; lli++)
		for (SIGNED_4B_TYPE llj = lli+1; llj < fdbn; llj++)
		if ( kdist(lli, llj) < 0 )
			kdist(lli, llj) = kdist(llj, lli) = dbase->get_indDistance(lli, llj, &dims, metrV, metrK);
	}

	//REALNUM_TYPE dbMx, dbMn, dbMn0, dbAv;
	matrix<REALNUM_TYPE> dsdist;
	dbase->get_distMatr(dsdist);
		
	dbase->set_distMatr(kdist);
	dbase->calc_dist_pars();

	apvector<REALNUM_TYPE> statz;
	dbase->get_NearNeibDistances(statz, k, r, kmode);

	//restore (quicker)
	dbase->set_distMatr(dsdist);
	
	//Rmax or R^2max calculation (max.allowed distance-cutoff)
	statz[0] += ADzCutoff * statz[2];
	if (ADbySquareD) return sqrt(statz[0]);

	return (statz[0]);
}

UNSIGNED_1B_TYPE knn::pick(SIGNED_4B_TYPE el, set &b, REALNUM_TYPE ADCutOff)
/*
//		picks candidates for el, stores them in apvNeibs[], apvDists[]
//		aprx_mod defines if lossy approximation for search of nearest neighbours should be used
//		0 - no
//		1 - cut-off by distance of full dimensions (for #dimensions > 20)
//		2 - cut-off by distance of individual dimensions (#dimensions < 20)
*/
{
	apvNeibs.resize(0); 
	apvDists.resize(0);
	if ( b.IsEmpty() || (dbase == NULL) || b.IsInSet(el)) 	return 0;

	UNSIGNED_1B_TYPE shared_kvotes = 0; //only for specific case of k+ candidates and KNN_WT_K_VOTE
	
	REALNUM_TYPE *FF = NULL;
	SIGNED_4B_TYPE d, c, cndN, *AA = NULL;
	apvector<SIGNED_4B_TYPE> apvCands;		
	
	if (aprx_mod == 0) //full recalculation of distance matrix
			b.GetList(apvCands);
	else
	{/* //implement lossy and faster version of finding nearest neighbours! //??!!??
		set setC, setC1;
		setC.Dump();
		for (d = 0; d < apvDims.length(); d++)	setC |= dimneib[el][d];
		setC1 = dimneib1[el][0];
		for (d = 1; d < apvDims.length(); d++)	setC1 &= dimneib1[el][d];
		setC |= setC1;
		setC.GetList(apvCands);		
	*/
	}
	cndN = apvCands.length();

	AA = GRAB_MEM_BLOCKS(SIGNED_4B_TYPE, cndN);
	FF	= GRAB_MEM_BLOCKS(REALNUM_TYPE, cndN);
	QSortScore = FF;

	for (c = 0; c < cndN; c++)
	{
		if ( kdist(el, apvCands[c]) < 0)
			kdist(el, apvCands[c]) = kdist(apvCands[c], el) = dbase->get_indDistance(el, apvCands[c], &dims, metrV, metrK); //metric etc

		FF[c] = kdist(el, apvCands[c]);
		AA[c] = c;
	}
	qsort(AA, (size_t)cndN, sizeof(SIGNED_4B_TYPE), QSortCompareGreater);
	QSortScore = NULL;
	d = k;
	if (k)		
	{
		while ( FF[ AA[d-1] ] == FF[ AA[d] ] ) { d++; if (d == cndN) break; };
		if ( (wt_mode & KNN_WT_K_ONLY) == 0 ) d = k;
		if (d > k)
		{
			c = k;
			while (--c) if ( FF[ AA[c] ] > FF[ AA[c-1] ] ) break;
			shared_kvotes = k - c;

			if (wt_mode & KNN_WT_K_FULLD)
			{//first, analyze extra-neighbors in full-D space (after that other modes may stil apply)
				apvDists.resize(d - c);
				apvNeibs.resize(d - c);					
				while (d-- > c)	
				{					
					//apvDists[d - c]	 = dbase->get_Distance(el, apvCands[AA[d]]);
					apvDists[d - c]	 = dbase->get_indDistance(el, apvCands[AA[d]], NULL, metrV, metrK);
					apvNeibs[d - c]	 = AA[d];
				}
				BubbleSort(apvDists, apvNeibs); //descending order
				d = apvDists.length() - shared_kvotes;
				while (d) if ( apvDists[d] < apvDists[d - 1] ) break; else d--;
				while (d < apvNeibs.length())	AA[c++] = apvNeibs[d++];
				d = c;
				if (d == k) shared_kvotes = 0; //no need to vote
			}//if ( wt_mode & KNN_WT_K_FULLD )

			if ( (wt_mode & KNN_WT_K_VOTE) == 0 )
			{
				shared_kvotes = 0;
				if ( (wt_mode & KNN_WT_K_ALL) == 0 ) d = k; //remove extra neighbors
			}
		} //if (d > k)
	}
	else	//radial cut-off
		while ( FF[ AA[d] ] < r ) d++;

	if (ADCutOff > 0)
	{//apply applicability domain cut-off to filter out some of the "nearest neighbours"
		
		if ((k > 1) && (ADmode & 2))
		{//compatibility mode, checks average k-dist vs. AD
			REALNUM_TYPE rtAVD = 0;
			bool SqrAD = ((ADmode & 1) == 0);
			for (c = 0; c < k; c++)				
				if (SqrAD)	rtAVD += sqr(FF[ AA[c] ]); else rtAVD += FF[ AA[c] ];
			rtAVD /= k;
			if (SqrAD) rtAVD = sqrt(rtAVD);
			if (rtAVD > ADCutOff) d = 0; //no prediction
		}
		else
		{//normal filtering by AD, final handling will be done on upper levels
			while (d) 
				if ( FF[ AA[d - 1] ] > ADCutOff ) d--; else break;
		}
	}//if (ADCutOff > 0)

	apvNeibs.resize(d);
	apvDists.resize(d);
	while (d--)
	{
		apvNeibs[d]	 = apvCands[AA[d]];
		apvDists[d]	 = kdist(el, apvNeibs[d]);
	}

	DROP_MEM_BLOCKS(FF);
	DROP_MEM_BLOCKS(AA);

	if (shared_kvotes && (apvNeibs.length() <= k)) shared_kvotes = 0; //to remove irrelevant ties (e.g. may happen due to AD-check)
	return (shared_kvotes);
}


void knn::predict(set &b, set &t, REALNUM_TYPE AD)
//		activity prediction of datapoints from t on the basis of datapoints in b
//		t and b must not overlap!
//
//NB:	it is inefficient to perform dataset reduction (by dimensions and by b subset) every time,
//		because there will be quite many such calls
//
{
	SIGNED_4B_TYPE d, c, np = pred_data.length();
	if ( (np == 0) || t.IsEmpty() || b.IsEmpty() || (dbase == NULL) ) 
		return;
		
	if ( SIGNED_4B_TYPE(kdist.ColNO()) != np) 
	{
		kdist.SetSize(np, np);
		for (d = 0; d < np; d++)
		{
			kdist(d,d) = 0;
			for (c = d + 1; c < np; c++)
				kdist(c,d) = kdist(d,c) = DUMMY_VAL;
		}
		if (pred_mode & KNN_PRED_CLASS)	iclass_sim.SetSize(0, 0);	//force recalculation of inter-class similarities	
	}

	calc_class_sim();	//recalculates inter-class similarities, if applicable

	UNSIGNED_1B_TYPE shared_kvotes; //only for specific case of k+ candidates and KNN_WT_K_VOTE
	SIGNED_4B_TYPE el;

	set setT(t);
	while (!setT.IsEmpty())
	{
		setT.GetElement(el);
		setT.RemoveFromSet(el);

		shared_kvotes = pick(el, b, AD); //no Applicability Domain when training model!

		if ((k > 1) && (AD > 0) && apvNeibs.length())
		{//applicability domain options: default is lax mode (even 1 kneib is enough)
			if ( ((ADmode & 8) && (apvNeibs.length() < k)) || //strict k mode
				((ADmode & 4) && (2*apvNeibs.length() < k)) ) //av. k mode
			 apvNeibs.resize(0);
		}

		if (apvNeibs.length() == 0) //no prediction
			continue;

		if ( (pred_mode & KNN_PRED_EXTERNAL_F) && (pred_f != NULL) )
		{//supply nearest neighbours and their distances to the external function
			kNNExternPredFunc pf = (kNNExternPredFunc)pred_f;
			pred_data[el] = pf(apvNeibs, apvDists);
			continue;
		}

		REALNUM_TYPE rtx, sumd = 0;

		//apvDists will be used to store weights
		if ( ((wt_mode & KNN_WT_F_NO) == 0) || (apvDists.length() == 1) || (k == 1) )
		{//no need to calculate weights
			for (d = 0; d < apvDists.length(); d++)	apvDists[d] = 1;
		}
		else
		{//we have d nearest neighbours, let us apply weights
			if (wt_mode & KNN_WT_ABS_D) 
				sumd = 1;
			else
			{
				for (sumd = d = 0; d < apvDists.length(); d++)	
				{
					if (k && shared_kvotes && (d == k)) break; //use only first k weights if votes are shared
					if (wt_mode & KNN_WT_D1) sumd += apvDists[d]; else sumd += sqr(apvDists[d]);
				}				
				if (sumd == 0.0) sumd = 1; //in case of all neighbors at "zero-distance"
			}

			if (wt_mode & KNN_WT_F_EXP)		//exponential weighing scheme
			for (d = 0; d < apvDists.length(); d++)
			{//wt_k here works as a smoothing coefficient, it is equiv to 1/(2c^2) of Gaussian
				rtx = apvDists[d]/sumd;
				apvDists[d] = exp( -wt_k * sqr(rtx) );
			}
			
			if (wt_mode & KNN_WT_F_HYP)
			for (d = 0; d < apvDists.length(); d++)
			{
				if ( wt_mode & KNN_WT_D1) rtx = apvDists[d]; else rtx = sqr(apvDists[d]);
				rtx /= sumd;
				rtx += 1;
				apvDists[d] = pow(rtx, -wt_k);
			}

			if (wt_mode & KNN_WT_F_MNK)		//default weighing scheme (linear, if wt_k == 1)
			for (d = 0; d < apvDists.length(); d++)
			{//NB: in case of 2 neighbors, where 1st is at zero-distance, the weight of 2nd will be 0!!
				if ( wt_mode & KNN_WT_D1) rtx = apvDists[d]; else rtx = sqr(apvDists[d]);
				rtx /= sumd;
				rtx = 1 - pow(rtx, wt_k);
				apvDists[d] = pow(rtx, 1/wt_k);
			}
		}//else..if ( (wt_mode & KNN_WT_F_NO) || (apvDists.length() == 1) )

		if (k && shared_kvotes)
		{//adjust weights of extra neighbours
			d	 = k - shared_kvotes;
			rtx	 = shared_kvotes;
			rtx	/= apvDists.length() - d;
			for (; d < apvDists.length(); d++)	apvDists[d] *= rtx;
		}

		pred_data[el] = 0;
		apvector<REALNUM_TYPE> repr(c_wts.length(), 0); //number of groups in KNN_PRED_CLASS
		
		//calculate sum of weights		
		sumd = SMALL_NUMBER; //important to prevent rare cases of pred.act to be > max.act due to rounding July 23 2009
		for (d = 0; d < apvDists.length(); d++)	sumd += apvDists[d];
		for (d = 0; d < apvDists.length(); d++)
		{
			rtx = apvDists[d] / sumd;
			if ( pred_mode & KNN_PRED_CLASS )
			{//also see Kovatcheva et al, J.Chem.Inf.Comput.Sci. [now JCIM] 2004, 44(2), 582-595
				for (c = 0; c < repr.length(); c++)
					repr[c] += iclass_sim(SIGNED_4B_TYPE(dbase->get_Act(apvNeibs[d])), c) * rtx;				
			}
			else
				pred_data[el] += dbase->get_Act(apvNeibs[d]) * rtx;
		}

		//finalize class memberships		
		if ( pred_mode & KNN_PRED_CLASS)	//assign to the maximum represented class			
		for (rtx = d = 0; d < repr.length(); d++)
		if ( repr[d] > rtx )
		{//if equal weights to first class with maximum population will be used
			rtx = repr[d];
			pred_data[el] = d;
		}

		//NB: for KNN_PRED_CATEG rounding will be done during evaluation

	}//while (!setT.IsEmpty())
}

REALNUM_TYPE knn::predict_ext(apvector<REALNUM_TYPE> &dp, REALNUM_TYPE AD)
/*description:	predicts one datapoint, AD is aplicability domain (max. distance)
				if no prediction is possible (out of AD) then return value is set 
				to >>max.act of the dataset
 precondition:	dp[] is array of dimensions for a datapoint being predicted
 postcondition:	predicted activity is returned
*/
{
	UNSIGNED_4B_TYPE i, N = dbase->get_Ndatapoints();
	
	kdist.SetSize(N + 1, N + 1);
	for (i = 0; i <= N; i++)	kdist(i, N) = kdist(N, i) = DUMMY_VAL;

	pred_data.resize(N + 1);
	pred_data[N] = 2*dbase->get_MaxAct() - dbase->get_MinAct();

	dbase->add_dp(dp, pred_data[N]);

	set bs(0, N), ts;
	ts.PutInSet(N);
	predict(bs, ts, AD);
	REALNUM_TYPE pp = pred_data[N];
	dbase->remove_dp(N);

	pred_data.resize(N);
	kdist.SetSize(N, N);
	return pp;
}

void knn::classify_predictions(apvector<REALNUM_TYPE> &p_data, apvector<SIGNED_4B_TYPE> &p_ndata, bool ifSimple)
//NB: only for category or class kNN!
{
	SIGNED_4B_TYPE i, j, np = p_data.length();
	p_ndata.resize(np);
	
	if (ifSimple)
	{//classes, etc.
		for (i = 0; i < np; i++) 
			p_ndata[i] = (SIGNED_4B_TYPE) p_data[i];
		return;
	}
	
	//rounding by break points
	for (i = 0; i < np; i++)
	{
		for (j = 0; j < c_bps.length(); j++)
		if (p_data[i] < c_bps[j]) 
		{
			p_ndata[i] = (SIGNED_4B_TYPE) floor(c_bps[j]);
			break;
		}
		if (j == c_bps.length()) p_ndata[i] = (SIGNED_4B_TYPE) ceil(c_bps[c_bps.length()-1]);
	}
}

REALNUM_TYPE knn::compare(apvector<REALNUM_TYPE> &e_data, apvector<REALNUM_TYPE> &p_data, bool extrnVal)
{
	SIGNED_4B_TYPE np = p_data.length();
	if (np != e_data.length()) return INVALID;

	if ( (eval_mode & KNN_PRED_EXTERNAL_F) && (eval_f != NULL) )
	{//perform external evaluation of the model
		kNNExternEvalFunc evf = (kNNExternEvalFunc)eval_f;
		return evf(e_data, p_data);		
	}

	if ( (pred_mode & KNN_PRED_CATEG) || (pred_mode & KNN_PRED_CLASS) )
	{
		apvector<SIGNED_4B_TYPE> exp_data_d, pred_data_d;
		classify_predictions(e_data, exp_data_d);
		classify_predictions(p_data, pred_data_d, (pred_mode & KNN_PRED_CLASS) > 0);

		SIGNED_4B_TYPE ngroups = c_wts.length();
		matrix<SIGNED_4B_TYPE> conf_mtx(ngroups, ngroups);
		qsarBLOCK.get_conf_mtx(exp_data_d, pred_data_d, conf_mtx);
		
		//calculate quality factors
		if (pred_mode & KNN_PRED_CLASS)
			return qsarBLOCK.get_ccr(conf_mtx, c_wts, pn, (eval_mode | KNN_EVAL_CLASS_ERR) & 0xFF);
		else
			return qsarBLOCK.get_ccr(conf_mtx, c_wts, pn, eval_mode & 0xFF );
	}

	if (extrnVal)
	{//it is assumed that the min.number of datapoints required for evaluation has been already checked
		if ((eval_mode & KNN_EVAL_EXT_R2) == KNN_EVAL_EXT_R2)
			return qsarBLOCK.correl(e_data, p_data);

		if ( (eval_mode & KNN_EVAL_EXT_Q2F) == KNN_EVAL_EXT_Q2F)
		{
			apvector<REALNUM_TYPE> trn_var;
			dbase->get_ActValues(trn_var);
			return qsarBLOCK.q2F13(e_data, p_data, qsarBLOCK.varianceV(trn_var, true), 0);
		}
	}

	return qsarBLOCK.q2etc(e_data, p_data, (eval_mode & 3) );
}

void knn::evaluate()
//		kNN-model's quality evaluation, based on pred_data[]
//		pred_data[] is filled by predict() calls
//
//NB:	kdist() matrix is used if non-empty and non-negative! 
//		Hence, kdist() should be reset before predictions every time dimensions are changed!
//		same applies for iclass_sim() matrix

//NB:	test-set evaluation is done separately!
{
	SIGNED_4B_TYPE i, np = pred_data.length();

	//prediction
	set setb(0, np);
	if ( (lgo == 1) || (lgo_sets.length() == 0) )
	{
		set sett, seti(setb);
		while (!seti.IsEmpty())
		{
			seti.GetElement(i);
			seti.RemoveFromSet(i);

			sett.PutInSet(i);
			setb.RemoveFromSet(i);

			predict(setb, sett);

			sett.RemoveFromSet(i);
			setb.PutInSet(i);
		}
	}
	else
	{
		set blgo;
		for (i = 0; i < lgo_sets.length(); i++)
		{
			blgo = setb - lgo_sets[i];
			predict(blgo, lgo_sets[i]);
		}
	}

	apvector<REALNUM_TYPE> exp_data;
	dbase->get_ActValues(exp_data);

	//evaluation
	qualV = compare(exp_data, pred_data);
	if (pred_mode & KNN_PRED_CLASS)
		if (class_sep_cf > 0)
			qualV += class_sep_cf*calc_mean_class_separation();
}

SIGNED_4B_TYPE knn::evaluate_ext(dataset &ext_ds, apvector<REALNUM_TYPE> &pr_act, bool calcqual, bool bylabel)
{
	apvector<SIGNED_4B_TYPE> neibs;
	return evaluate_ext(ext_ds, pr_act, neibs, calcqual, bylabel);
}

SIGNED_4B_TYPE knn::evaluate_ext(dataset &ext_ds, apvector<REALNUM_TYPE> &pr_act, apvector<SIGNED_4B_TYPE> &pr_neibs, bool calcqual, bool bylabel)
/*
	description:	predicts datapoints of external dataset ext_ds against current model
					calculates Applicability Domain based on current knn settings

	precondition:	ext_ds data should have same dimensions and 
					normalized the same way as the training dataset.

	postcondition:	calculated quality of external data prediction is stored in qualV
					returns actual number of datapoints that received prediction
					returns predictions in pr_act[] array, values > training set's MaxAct() mean that prediction is N/A
*/
{
	SIGNED_4B_TYPE i, dpn = ext_ds.get_Ndatapoints();	
	REALNUM_TYPE rtAD = INVALID;

	UNSIGNED_2B_TYPE nn = 0; //counter of nearest neighbor entries
	pr_neibs.resize(dpn * max(k, 2)); //array to store nearest neighbors
	
	if (k) rtAD = calc_AD();
	calc_class_sim();	//recalculates inter-class similarities, if applicable
	
	apvector<REALNUM_TYPE> dscrs, ex, pr(dpn);
	ext_ds.get_ActValues(ex);

	qualV = 0;

	for (i = 0; i < dpn; i++)
	{
		ext_ds.get_DimValues(i, dscrs);	
		
		//----------------------------------------------------
		if (bylabel)
		{//validate by labels and reorganize values -> helpful when descriptor order is different
			apvector<SIGNED_4B_TYPE> tdim;			
			dims.GetList(tdim);
			SIGNED_4B_TYPE ndims = tdim.length();
			apvector<REALNUM_TYPE> xdim(ndims, 0);
			SIGNED_4B_TYPE ii;
			for (ii = 0; ii < ndims; ii++)	xdim[ii] = dscrs[ ext_ds.get_dscr_pos(dbase->get_dscr(tdim[ii])) ];
			dscrs.resize(dbase->get_Ndimensions());
			for (ii = 0; ii < ndims; ii++) dscrs[tdim[ii]] = xdim[ii];			
		}
		//----------------------------------------------------

		pr[i] = predict_ext(dscrs, rtAD);
		
		SIGNED_4B_TYPE tb, napn = apvNeibs.length();
		if ( pr_neibs.length() < (nn + (napn << 1)) ) pr_neibs.resize(nn + dpn + (napn << 1) );
		for (tb = 0; tb < napn; tb++)	{	pr_neibs[nn++] = i;	pr_neibs[nn++] = apvNeibs[tb];	};
	} //for i
	
	pr_neibs.resize(nn);
	pr_act = pr;	//to return predicted values

	if (calcqual)
	{
		qsarBLOCK.remove_unpredicted(ex, pr, dbase->get_MaxAct());
		dpn = pr.length();
		if ((eval_mode & KNN_EVAL_EXT_Q2F) == KNN_EVAL_EXT_Q2F)
		{
			if (dpn == 0)	return 0;		
		}
		else
			if (dpn < 3)	return dpn;	//requires minimum 3 data points
			
		qualV = compare(ex, pr, true);
	}

	return dpn;
}
