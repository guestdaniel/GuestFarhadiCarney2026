#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> 
#include <mex.h>
#include <time.h>

#include "complex.hpp"

#define MAXSPIKES 1000000
#ifndef TWOPI
#define TWOPI 6.28318530717959
#endif

#ifndef __max
#define __max(a,b) (((a) > (b))? (a): (b))
#endif

#ifndef __min
#define __min(a,b) (((a) < (b))? (a): (b))
#endif

/**
 *
 * Mex function for `sim_efferent_model`
 *
 * Mex function that parses MATLAB inputs, appropriately converts them into a form 
 * suitable for passing to the model C code, calls the model on the converted inputs, and 
 * then appropriately converts the outputs and passes them back to MATLAB.
 *
 * Note that there are no safety checks built into this function --- all checking should
 * happen in MATLAB *before* passing inputs to this function.
 *
 * @param nlhs Number of return values (i.e., [n]umber [l]eft [h]and [s]ide values)
 * @param plhs mxArray of pointers to output variables (i.e., [p]ointers to [l]eft [h]and
 * [s]ide values). Obtain a pointer with `mexGetPr`.
 * @param nrhs Number of return values (i.e., [n]umber [l]eft [h]and [s]ide values)
 * @param prhs mxArray of pointers to input variables (i.e., [p]ointers to [r]eft [h]and
 * [s]ide values). Obtain a pointer with `mexGetPr`.
 *
 * The inputs are as follows on the RHS:
 * - [1] px double *
 * - [2] randNums_hsr double **
 * - [3] randNums_lsr double **
 * - [4] cf double *
 * - [5] n_chan int
 * - [6] tdres double
 * - [7] coch double *
 * - [8] cihc double *
 * - [9] species double
 * - [10] moc_cutoff double
 * - [11] moc_beta double * 
 * - [12] moc_offset double *
 * - [13] moc_weight double *
 * - [14] moc_width double
 * - [15] powerlaw_mode double
 * - [16] moc_minrate double
 * - [17] moc_maxrate double
 * - [18] dur_settle double
 * - [19] moc_delay double
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* Declare signature for `model`, the C function for the efferent model. */
	void model(
		// Principal inputs
		double *,   // px
		double **,  // randNums_hsr
		double **,  // randNums_lsr
		double *,   // cf
		// IHC/AN parameters
		int,        // n_chan
		double,     // tdres
		int,        // totalstim
		double *,   // cohc
		double *,   // cihc
		int,        // species
		double,     // spont
		int,        // powerlaw_mode
		// CN parameters
		double,     // cn_tau_e
		double,     // cn_tau_i
		double,     // cn_delay
		double,     // cn_amp
		double,     // cn_inh
		// IC parameters
		double,     // ic_tau_e
		double,     // ic_tau_i
		double,     // ic_delay
		double,     // ic_amp
		double,     // ic_inh
		// MOC parameters
		double,     // moc_cutoff
		double *,   // moc_beta
		double *,   // moc_offset
		double,     // moc_minrate
		double,     // moc_maxrate
		double *,   // moc_weight
		double,     // moc_width
		double,     // dur_settle
		double,     // moc_delay
		int,        // moc_fix_gain
		// BM/IHC outputs
		double **,  // controlout
		double **,  // c1out
		double **,  // c2out
		double **,  // ihcout
		// HSR outputs
		double **,  // expout_hsr
		double **,  // sout1_hsr
		double **,  // sout2_hsr
		double **,  // synout_hsr
		// LSR outputs
		double **,  // expout_lsr
		double **,  // sout1_lsr
		double **,  // sout2_lsr
		double **,  // synout_lsr
		// AN outputs
		double **,  // anrateout_hsr
		double **,  // anrateout_lsr
		// CN/IC outputs
		double **,  // cnout
		double **,  // icout
		// MOC outputs
		double **,  // mocwdr
		double **,  // mocic
		double **,  // gain
		double **   // gain_postmix
	);
							
	/* Check for proper number of arguments */
	if (nrhs != 19) 
	{
		mexErrMsgTxt("model requires 19 input arguments.");
	}; 

	if (nlhs != 4)  
	{
		mexErrMsgTxt("model requires 4 output argument.");
	};

	/* Get sizes of inputs */
	int totalstim = mxGetN(prhs[0]);

	/* De-reference (and, where needed, cast to int) scalar input values */
	int n_chan = (int) *mxGetPr(prhs[4]);
    double tdres = *mxGetPr(prhs[5]);
	int species = (int) *mxGetPr(prhs[8]);
	double moc_cutoff = *mxGetPr(prhs[9]);
	double moc_width = *mxGetPr(prhs[13]);
	int powerlaw_mode = (int) *mxGetPr(prhs[14]);
	double moc_minrate = *mxGetPr(prhs[15]);
	double moc_maxrate = *mxGetPr(prhs[16]);
	double dur_settle = *mxGetPr(prhs[17]);
	double moc_delay = *mxGetPr(prhs[18]);

	/* 
	 * Handle pressure vector by copying data from MATLAB mxArray pointer to dynamically
	 * allocated array in C
	 */
	double *px_mex = mxGetPr(prhs[0]);
	double *px = (double*) calloc(totalstim, sizeof(double));
	for (int i = 0; i < totalstim; i++) {
		px[i] = px_mex[i];
	}

	/* 
	 * Handle CF vector by copying data from MATLAB mxArray pointer to dynamically
	 * allocated array in C
	 */
	double *cf_mex = mxGetPr(prhs[3]);
	double *cf = (double*) calloc(n_chan, sizeof(double));
	for (int i = 0; i < n_chan; i++) {
		cf[i] = cf_mex[i];
	}

	/* 
	 * Handle COHC/CIHC/beta/offset/weight vectors by copying data from 
	 * MATLAB mxArray pointer to dynamically allocated array in C
	 */
	// Pull pointers for RHS arguments
	double *cohc_mex = mxGetPr(prhs[6]);
	double *cihc_mex = mxGetPr(prhs[7]);
	double *moc_beta_mex = mxGetPr(prhs[10]);
	double *moc_offset_mex = mxGetPr(prhs[11]);
	double *moc_weight_mex = mxGetPr(prhs[12]);
	// Allocate storage for equivalent vectors
	double *cohc = (double*) calloc(n_chan, sizeof(double));
	double *cihc = (double*) calloc(n_chan, sizeof(double));
	double *moc_beta = (double*) calloc(n_chan, sizeof(double));
	double *moc_offset = (double*) calloc(n_chan, sizeof(double));
	double *moc_weight = (double*) calloc(n_chan, sizeof(double));
	for (int i = 0; i < n_chan; i++) {
		cohc[i] = cohc_mex[i];
		cihc[i] = cihc_mex[i];
		moc_beta[i] = moc_beta_mex[i];
		moc_offset[i] = moc_offset_mex[i];
		moc_weight[i] = moc_weight_mex[i];
	}

	/* 
	 * Handle fGn matrices by copying data from MATLAB mxArray pointer to dynamically
	 * allocated matrix in C. Note that because MATLAB arrays are in column-major order 
	 * while C expects a row-major order, we must reorder the elements in memory before 
	 * passing to C.
	 */
	// Grab input pointers and allocate memory for row-major ordered noise matrices
	double *randNums_hsrtmp	= mxGetPr(prhs[1]);
	double *randNums_lsrtmp = mxGetPr(prhs[2]);
	double *randNums_hsr[n_chan];
	double *randNums_lsr[n_chan];
	for (int i = 0; i < n_chan; i++) {
        randNums_hsr[i] = (double*) calloc(totalstim, sizeof(double));
        randNums_lsr[i] = (double*) calloc(totalstim, sizeof(double));
    }

	// Loop through elements via linear indexing, store	elements in row-major order
	int idx_cf, idx_t;
	for (int i = 0; i < (totalstim * n_chan); i++) {
		// Determine Cartesian indices for linear index i given C-style row-major order
		idx_cf = (int) fmod(i, n_chan);  // index into channels
		idx_t = (int) (i/n_chan);        // index into time/samples

		// Store i-th element in corresponding location in randNums matrices
		randNums_hsr[idx_cf][idx_t] = randNums_hsrtmp[i];
		randNums_lsr[idx_cf][idx_t] = randNums_lsrtmp[i];
	}

	/* 
	 * Handle return matrices. All return matrices are size of (n_chan, totalstim). These 
	 * are dynamically allocated and passed to the C routine as pointers. Then, below,
	 * the data stored in these matrices will be extracted and reformatted into an mxArray
	 * for return as a pointer to MATLAB.
	 */	
	// Allocate output matrices
	double *controlout[n_chan];
	double *c1out[n_chan];
	double *c2out[n_chan];
	double *ihcout[n_chan];
	double *expout_hsr[n_chan];
	double *sout1_hsr[n_chan];
	double *sout2_hsr[n_chan];
	double *synout_hsr[n_chan];
	double *expout_lsr[n_chan];
	double *sout1_lsr[n_chan];
	double *sout2_lsr[n_chan];
	double *synout_lsr[n_chan];
	double *anrateout_hsr[n_chan];
	double *anrateout_lsr[n_chan];
	double *cnout[n_chan];
	double *icout[n_chan];
	double *mocwdr[n_chan];
	double *mocic[n_chan];
	double *gain[n_chan];
	double *gainpostmix[n_chan];
 	for (int i = 0; i < n_chan; i++) {
    	controlout[i]    = (double*) calloc(totalstim, sizeof(double));
    	c1out[i]         = (double*) calloc(totalstim, sizeof(double));
    	c2out[i]         = (double*) calloc(totalstim, sizeof(double));
    	ihcout[i]        = (double*) calloc(totalstim, sizeof(double));
    	expout_hsr[i]    = (double*) calloc(totalstim, sizeof(double));
    	sout1_hsr[i]     = (double*) calloc(totalstim, sizeof(double));
    	sout2_hsr[i]     = (double*) calloc(totalstim, sizeof(double));
    	synout_hsr[i]    = (double*) calloc(totalstim, sizeof(double));
    	expout_lsr[i]    = (double*) calloc(totalstim, sizeof(double));
    	sout1_lsr[i]     = (double*) calloc(totalstim, sizeof(double));
    	sout2_lsr[i]     = (double*) calloc(totalstim, sizeof(double));
    	synout_lsr[i]    = (double*) calloc(totalstim, sizeof(double));
		anrateout_hsr[i] = (double*) calloc(totalstim, sizeof(double));
		anrateout_lsr[i] = (double*) calloc(totalstim, sizeof(double));
		cnout[i]         = (double*) calloc(totalstim, sizeof(double));
		icout[i]         = (double*) calloc(totalstim, sizeof(double));
		mocwdr[i]        = (double*) calloc(totalstim, sizeof(double));
		mocic[i]         = (double*) calloc(totalstim, sizeof(double));
		gain[i]          = (double*) calloc(totalstim, sizeof(double));
		gainpostmix[i]   = (double*) calloc(totalstim, sizeof(double));
	}

	// Set gain and gain postmix to values of 1 by default
	for (int i = 0; i < n_chan; i++) {
		for (int j = 0; j < totalstim; j++) {
			gain[i][j] = 1.0;
			gainpostmix[i][j] = 1.0;
		}
	}

	// Run the efferent model via external call to C code
	model(
		px, 
		randNums_hsr,
		randNums_lsr,
		cf,
		n_chan,
		tdres,
		totalstim,
		cohc,
		cihc,
		species,
		100.0,                // spont
		powerlaw_mode,
		0.5e-3,               // default value for cn_tau_e (unused)
		2.0e-3,               // default value for cn_tau_i (unused)
		1.0e-3,               // default value for cn_delay (unused)
		1.5,                  // default value for cn_A (unused)
		0.6,                  // default value for cn_I (unused)
		1.0/(10.0*64.0),      // default value for ic_tau_e (unused)
		1.0/(10.0*64.0)*1.5,  // default value for ic_tau_i (unused)
		1.0/(10.0*64.0)*2.0,  // default value for ic_delay (unused)
		1.0,                  // default value for ic_A (unused)
		0.9,                  // default value for ic_I (unused)
		moc_cutoff,
		moc_beta,
		moc_offset,
		moc_minrate,
		moc_maxrate,
		moc_weight,
		moc_width,
		dur_settle,
		moc_delay,
		0,                    // moc_fix_gain = false
		controlout,
		c1out,
		c2out,
		ihcout,
		expout_hsr,
		sout1_hsr,
		sout2_hsr,
		synout_hsr,
		expout_lsr,
		sout1_lsr,
		sout2_lsr,
		synout_lsr,
		anrateout_hsr,
		anrateout_lsr,
		cnout,
		icout,
		mocwdr,
		mocic,
		gain,
		gainpostmix
	); 

	/*
	 * Handle returning results to MATLAB. First, we use mxCreateNumericArray to create
	 * appropriately sized output mxArrays. 
	 */
	// Create arrays for output
    mwSize size_output[2] = {n_chan, totalstim};
	plhs[0] = mxCreateNumericArray(2, size_output, mxDOUBLE_CLASS, mxREAL);    
    plhs[1] = mxCreateNumericArray(2, size_output, mxDOUBLE_CLASS, mxREAL); 
    plhs[2] = mxCreateNumericArray(2, size_output, mxDOUBLE_CLASS, mxREAL);  
    plhs[3] = mxCreateNumericArray(2, size_output, mxDOUBLE_CLASS, mxREAL); 

	// Obtain pointers to said arrays
	double *ihcouttmp = mxGetPr(plhs[0]);
	double *anrateout_hsrtmp = mxGetPr(plhs[1]);
	double *anrateout_lsrtmp = mxGetPr(plhs[2]);
	double *gainpostmixtmp = mxGetPr(plhs[3]);

	// Loop through elements via linear indexing, store	elements in column-major order
   	for (int i = 0; i < (n_chan * totalstim); i++) {
		// Determine Cartesian indices for linear index i given C-style row-major order
		idx_cf = (int) fmod(i, n_chan);  // index into channels
		idx_t = (int) (i/n_chan);        // index into time/samples

		// Store element at [idx_cf, idx_t] in i-th position in output matrix
		ihcouttmp[i] = ihcout[idx_cf][idx_t];
		anrateout_hsrtmp[i] = anrateout_hsr[idx_cf][idx_t];
		anrateout_lsrtmp[i] = anrateout_lsr[idx_cf][idx_t];
		gainpostmixtmp[i] = gainpostmix[idx_cf][idx_t];
	}

	// Free all dynamically allocated memory
	for (int i = 0; i < n_chan; i++) {
		free(randNums_hsr[i]);
		free(randNums_lsr[i]);
		free(controlout[i]);
		free(c1out[i]);
		free(c2out[i]);
		free(ihcout[i]);
		free(expout_hsr[i]);
		free(sout1_hsr[i]);
		free(sout2_hsr[i]);
		free(synout_hsr[i]);
		free(expout_lsr[i]);
		free(sout1_lsr[i]);
		free(sout2_lsr[i]);
		free(synout_lsr[i]);
		free(anrateout_hsr[i]);
		free(anrateout_lsr[i]);
		free(cnout[i]);
		free(icout[i]);
		free(mocwdr[i]);
		free(mocic[i]);
		free(gain[i]);
		free(gainpostmix[i]);
	}
	free(px);
	free(cf);
	free(cohc);
	free(cihc);
	free(moc_beta);
	free(moc_weight);
	free(moc_offset);
}