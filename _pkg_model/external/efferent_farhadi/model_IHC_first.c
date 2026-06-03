/* This is the BEZ2018 version of the code for auditory periphery model from the Carney, Bruce and Zilany labs.
 * This is for the first tw0o samples of efferent model 
 * This release implements the version of the model described in:
 *
 *   Bruce, I.C., Erfani, Y., and Zilany, M.S.A. (2018). "A Phenomenological
 *   model of the synapse between the inner hair cell and auditory nerve: 
 *   Implications of limited neurotransmitter release sites," to appear in
 *   Hearing Research. (Special Issue on "Computational Models in Hearing".)
 *
 * Please cite this paper if you publish any research
 * results obtained with this code or any modified versions of this code.
 *
 * See the file readme.txt for details of compiling and running the model.
 *
 * %%% Ian C. Bruce (ibruce@ieee.org), Yousof Erfani (erfani.yousof@gmail.com),
 *     Muhammad S. A. Zilany (msazilany@gmail.com) - December 2017 %%%
 *
 */

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

/* This function is the MEX "wrapper", to pass the input and output variables between the .mex* file and Matlab or Octave */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	
	double *px, cf, tdres, reptime, *cohc, cihc, *statein;
	int    nrep, pxbins, lp, cp, s,totalstim, species;
    mwSize outsize[2], outsize2[2];
    
	double *pxtmp, *cftmp, *nreptmp, *tdrestmp, *reptimetmp, *cohctmp, *cihctmp, *speciestmp, *stateintmp;
    double *ihcout, *state;
   
	void   IHCAN(double *, double, int, double, int, double *, double, int, double *, double *, double *);
	
	/* Check for proper number of arguments */
	
	if (nrhs != 9) 
	{
		mexErrMsgTxt("model_IHC requires 9 input arguments.");
	}; 

	if (nlhs !=2)  
	{
		mexErrMsgTxt("model_IHC requires 2 output argument.");
	};
	
	/* Assign pointers to the inputs */

	pxtmp		= mxGetPr(prhs[0]);
	cftmp		= mxGetPr(prhs[1]);
	nreptmp		= mxGetPr(prhs[2]);
	tdrestmp	= mxGetPr(prhs[3]);
	reptimetmp	= mxGetPr(prhs[4]);
    cohctmp		= mxGetPr(prhs[5]);
    cihctmp		= mxGetPr(prhs[6]);
    speciestmp	= mxGetPr(prhs[7]);
	stateintmp		= mxGetPr(prhs[8]);
	/* Check individual input arguments */

	pxbins = (int) mxGetN(prhs[0]);
	
	if (pxbins==1)
		mexErrMsgTxt("px must be a row vector\n");

    species = (int) speciestmp[0];
	if (speciestmp[0]!=species)
		mexErrMsgTxt("species must an integer.\n");
	if (species<1 || species>3)
		mexErrMsgTxt("Species must be 1 for cat, or 2 or 3 for human.\n");
 
    
	cf = cftmp[0];
    if (species==1)
    {
	if ((cf<124.9)|(cf>40.1e3))
	{
		mexPrintf("cf (= %1.1f Hz) must be between 125 Hz and 40 kHz for cat model\n",cf);
		mexErrMsgTxt("\n");
    }
    }
    if (species>1)
    {
  	if ((cf<124.9)|(cf>20.1e3))
	{
		mexPrintf("cf (= %1.1f Hz) must be between 125 Hz and 20 kHz for human model\n",cf);
		mexErrMsgTxt("\n");
    }
    }
    
	nrep = (int) nreptmp[0];
	if (nreptmp[0]!=nrep)
		mexErrMsgTxt("nrep must an integer.\n");
	if (nrep<1)
		mexErrMsgTxt("nrep must be greater that 0.\n");

    tdres = tdrestmp[0];
	
	reptime = reptimetmp[0];
	if (reptime<pxbins*tdres)  /* duration of stimulus = pxbins*tdres */
		mexErrMsgTxt("reptime should be equal to or longer than the stimulus duration.\n");

   
		
   
	

	cihc = cihctmp[0]; /* impairment in the IHC  */
	if ((cihc<0)|(cihc>1))
	{
		mexPrintf("cihc (= %1.1f) must be between 0 and 1\n",cihc);
		mexErrMsgTxt("\n");
	}
	
   
	/* Calculate number of samples for total repetition time */

	/*totalstim = (int)floor((reptime*1e3)/(tdres*1e3)); */ /*older definition*/
    totalstim = (int)floor(reptime/tdres+0.5);
    px = (double*)mxCalloc(totalstim,sizeof(double));
    cohc = (double*)mxCalloc(totalstim,sizeof(double)); 
	statein =(double*)mxCalloc(105,sizeof(double)); 
	/* Put stimulus waveform into pressure waveform */

	for (lp=0; lp<pxbins; lp++)
	{
		px[lp] = pxtmp[lp];
		    cohc[lp] = cohctmp[lp];
		if ((cohc[lp]<0)|(cohc[lp]>1))
		{
			mexPrintf("cohc (= %1.1f) must be between 0 and 1\n",cohc[lp]);
			mexErrMsgTxt("\n");
		}
	}
	/* Create an array for the return argument */
	for (s=0; s<105; s++)
	{
		statein[s]=stateintmp[s];
	}
	
    outsize[0] = 1;
	outsize[1] = totalstim*nrep;
    outsize2[0] = 1;
	outsize2[1] = 105;   // change this when adding states
	plhs[0] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericArray(2, outsize2, mxDOUBLE_CLASS, mxREAL);
	
	/* Assign pointers to the outputs */
	
	ihcout = mxGetPr(plhs[0]);
	state  = mxGetPr(plhs[1]);	
	
	/* run the model */

	IHCAN(px,cf,nrep,tdres,totalstim,cohc,cihc,species,ihcout,state,statein);

 mxFree(px);

}

void IHCAN(double *px, double cf, int nrep, double tdres, int totalstim,
                double *cohc, double cihc, int species, double *ihcout,  double *state, double *statein)
{	
    
    /*variables for middle-ear model */
	double megainmax;
    double *mey1, *mey2, *mey3, meout,c1filterouttmp,c2filterouttmp,c1vihctmp,c2vihctmp, phasein, *rin, *ohcin, *Ihcin;
    double fp,C,m11,m12,m13,m14,m15,m16,m21,m22,m23,m24,m25,m26,m31,m32,m33,m34,m35,m36;

	/*variables for the signal-path, control-path and onward */
	double *ihcouttmp,*tmpgain, *statetmp;
	int    grd;
	double C1output11[1],C1output21[1],C1output31[1],C1output41[1],C1output51[1],C1output12[1],C1output22[1],C1output32[1],C1output42[1],C1output52[1],C1input11[1], C1input12[1], C1input13[1], C1input21[1], C1input22[1], C1input23[1], C1input31[1], C1input32[1], C1input33[1], C1input41[1], C1input42[1], C1input43[1], C1input51[1], C1input52[1], C1input53[1], C1input61[1], C1input62[1], C1input63[1],c1gainout[1], c1phaseout[1];
    double bmplace,centerfreq,gain,taubm,ratiowb,bmTaubm,fcohc,TauWBMax,TauWBMin,tauwb,C2output11[1],C2output21[1],C2output31[1],C2output41[1],C2output51[1],C2output12[1],C2output22[1],C2output32[1],C2output42[1],C2output52[1];
    double Taumin[1],Taumax[1],bmTaumin[1],bmTaumax[1],ratiobm[1],lasttmpgain,wbgain,ohcasym,ihcasym,delay,wbphase[1],r0[1],r1[1],r2[1],r3[1],I0[1],I1[1],I2[1],I3[1],Ohc3[1],Ohc2[1],Ohc1[1],Ohc0[1],Ohcl3[1],Ohcl2[1],Ohcl1[1],Ohcl0[1], c2gainout[1], c2phaseout[1], C2input11[1], C2input12[1], C2input13[1], C2input21[1], C2input22[1], C2input23[1], C2input31[1], C2input32[1], C2input33[1], C2input41[1], C2input42[1], C2input43[1], C2input51[1], C2input52[1], C2input53[1], C2input61[1], C2input62[1], C2input63[1];
	int    i,n,delaypoint,grdelay[1],bmorder,wborder;
	double wbout1,wbout,ohcnonlinout,ohcout,tmptauc1,tauc1,rsigma,wb_gain;
	double Ihc3[1],Ihc2[1],Ihc1[1],Ihc0[1],Ihcl3[1],Ihcl2[1],Ihcl1[1],Ihcl0[1],Ihc7[1],Ihc6[1],Ihc5[1],Ihc4[1],Ihcl7[1],Ihcl6[1],Ihcl5[1],Ihcl4[1];
            
    /* Declarations of the functions used in the program */
	double C1ChirpFilt(double, double,double, int, double, double, double *, double *, double *, double *, double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *);
	double C2ChirpFilt(double, double,double, int, double, double, double *, double *, double *, double *, double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *);
    double WbGammaTone(double, double, double, int, double, double, int, double *, double, double *, double *, double *, double *, double *, double *,double *, double *,double *);

    double Get_tauwb(double, int, int, double *, double *);
	double Get_taubm(double, int, double, double *, double *, double *);
    double gain_groupdelay(double, double, double, double, int *);
    double delay_cat(double cf);
    double delay_human(double cf);

    double OhcLowPass(double, double, double, int, double, int,double *, double *, double *, double *, double *, double *,double *, double *,double *);
    double IhcLowPass(double, double, double, int, double, int,double *, double *, double *, double *, double *, double *,double *, double *,double *,double *, double *, double *, double *, double *,double *, double *,double *);
	double Boltzman(double, double, double, double, double);
    double NLafterohc(double, double, double, double);
	double ControlSignal(double, double, double, double, double);

    double NLogarithm(double, double, double, double);
    
    /* Allocate dynamic memory for the temporary variables */
	ihcouttmp  = (double*)mxCalloc(totalstim*nrep,sizeof(double));
	statetmp   = (double*)mxCalloc(105,sizeof(double));
	//statetmp2   = (double*)mxCalloc(totalstim,sizeof(double));
	mey1 = (double*)mxCalloc(totalstim,sizeof(double));
	mey2 = (double*)mxCalloc(totalstim,sizeof(double));
	mey3 = (double*)mxCalloc(totalstim,sizeof(double));
	rin  = (double*)mxCalloc(8,sizeof(double));
	ohcin =(double*)mxCalloc(8,sizeof(double));
	Ihcin =(double*)mxCalloc(16,sizeof(double));
	tmpgain = (double*)mxCalloc(totalstim,sizeof(double));
    
	/** Calculate the center frequency for the control-path wideband filter
	    from the location on basilar membrane, based on Greenwood (JASA 1990) */

	if (species==1) /* for cat */
    {
        /* Cat frequency shift corresponding to 1.2 mm */
        bmplace = 11.9 * log10(0.80 + cf / 456.0); /* Calculate the location on basilar membrane from CF */
        centerfreq = 456.0*(pow(10,(bmplace+1.2)/11.9)-0.80); /* shift the center freq */
    }

	if (species>1) /* for human */
    {
        /* Human frequency shift corresponding to 1.2 mm */
        bmplace = (35/2.1) * log10(1.0 + cf / 165.4); /* Calculate the location on basilar membrane from CF */
        centerfreq = 165.4*(pow(10,(bmplace+1.2)/(35/2.1))-1.0); /* shift the center freq */
    }
    
	/*==================================================================*/
	/*====== Parameters for the gain ===========*/
    
	if(species==1) gain = 52.0/2.0*(tanh(2.2*log10(cf/0.6e3)+0.15)+1.0); /* for cat */
    if(species>1) gain = 52.0/2.0*(tanh(2.2*log10(cf/0.6e3)+0.15)+1.0); /* for human */
    /*gain = 52/2*(tanh(2.2*log10(cf/1e3)+0.15)+1);*/
    if(gain>60.0) gain = 60.0;  
    if(gain<15.0) gain = 15.0;
    
	/*====== Parameters for the control-path wideband filter =======*/
	bmorder = 3;
	Get_tauwb(cf,species,bmorder,Taumax,Taumin);
	taubm   = cohc[0]*(Taumax[0]-Taumin[0])+Taumin[0];
	ratiowb = Taumin[0]/Taumax[0];
	/*====== Parameters for the signal-path C1 filter ======*/
	Get_taubm(cf,species,Taumax[0],bmTaumax,bmTaumin,ratiobm);
	bmTaubm  = cohc[0]*(bmTaumax[0]-bmTaumin[0])+bmTaumin[0];
	fcohc    = bmTaumax[0]/bmTaubm;
    /*====== Parameters for the control-path wideband filter =======*/
	wborder  = 3;
    TauWBMax = Taumin[0]+0.2*(Taumax[0]-Taumin[0]);
	TauWBMin = TauWBMax/Taumax[0]*Taumin[0];
    tauwb    = TauWBMax+(bmTaubm-bmTaumax[0])*(TauWBMax-TauWBMin)/(bmTaumax[0]-bmTaumin[0]);
	
	wbgain = gain_groupdelay(tdres,centerfreq,cf,tauwb,grdelay);
	tmpgain[0]   = wbgain; 
	lasttmpgain  = wbgain;
  	/*===============================================================*/
    /* Nonlinear asymmetry of OHC function and IHC C1 transduction function*/
	ohcasym  = 7.0;    
	ihcasym  = 3.0;
  	/*===============================================================*/
    /*===============================================================*/
    /* Prewarping and related constants for the middle ear */
     fp = 1e3;  /* prewarping frequency 1 kHz */
     C  = TWOPI*fp/tan(TWOPI/2*fp*tdres);
     if (species==1) /* for cat */
     {
         /* Cat middle-ear filter - simplified version from Bruce et al. (JASA 2003) */
         m11 = C/(C + 693.48);                    m12 = (693.48 - C)/C;            m13 = 0.0;
         m14 = 1.0;                               m15 = -1.0;                      m16 = 0.0;
         m21 = 1/(pow(C,2) + 11053*C + 1.163e8);  m22 = -2*pow(C,2) + 2.326e8;     m23 = pow(C,2) - 11053*C + 1.163e8; 
         m24 = pow(C,2) + 1356.3*C + 7.4417e8;    m25 = -2*pow(C,2) + 14.8834e8;   m26 = pow(C,2) - 1356.3*C + 7.4417e8;
         m31 = 1/(pow(C,2) + 4620*C + 909059944); m32 = -2*pow(C,2) + 2*909059944; m33 = pow(C,2) - 4620*C + 909059944;
         m34 = 5.7585e5*C + 7.1665e7;             m35 = 14.333e7;                  m36 = 7.1665e7 - 5.7585e5*C;
         megainmax=41.1405;
     };
     if (species>1) /* for human */
     {
         /* Human middle-ear filter - based on Pascal et al. (JASA 1998)  */
         m11=1/(pow(C,2)+5.9761e+003*C+2.5255e+007);m12=(-2*pow(C,2)+2*2.5255e+007);m13=(pow(C,2)-5.9761e+003*C+2.5255e+007);m14=(pow(C,2)+5.6665e+003*C);             m15=-2*pow(C,2);					m16=(pow(C,2)-5.6665e+003*C);
         m21=1/(pow(C,2)+6.4255e+003*C+1.3975e+008);m22=(-2*pow(C,2)+2*1.3975e+008);m23=(pow(C,2)-6.4255e+003*C+1.3975e+008);m24=(pow(C,2)+5.8934e+003*C+1.7926e+008); m25=(-2*pow(C,2)+2*1.7926e+008);	m26=(pow(C,2)-5.8934e+003*C+1.7926e+008);
         m31=1/(pow(C,2)+2.4891e+004*C+1.2700e+009);m32=(-2*pow(C,2)+2*1.2700e+009);m33=(pow(C,2)-2.4891e+004*C+1.2700e+009);m34=(3.1137e+003*C+6.9768e+008);     m35=2*6.9768e+008;				m36=(-3.1137e+003*C+6.9768e+008);
         megainmax=2;
     };
  	for (n=0;n<totalstim;n++) /* Start of the loop */
    {    
        
		
		// this part is for one sample first
	

	        if (n==0)  /* Start of the middle-ear filtering section  */
		{
	    	mey1[0]  = m11*px[0];
            if (species>1) mey1[0] = m11*m14*px[0];
            mey2[0]  = mey1[0]*m24*m21;
            mey3[0]  = mey2[0]*m34*m31;
            meout = mey3[0]/megainmax ;
        }
            
        else if (n==1)
		{
            mey1[1]  = m11*(-m12*mey1[0] + px[1]       - px[0]);
            if (species>1) mey1[1] = m11*(-m12*mey1[0]+m14*px[1]+m15*px[0]);
			mey2[1]  = m21*(-m22*mey2[0] + m24*mey1[1] + m25*mey1[0]);
            mey3[1]  = m31*(-m32*mey3[0] + m34*mey2[1] + m35*mey2[0]);
            meout = mey3[1]/megainmax;
		}
			
          
		
 

		/*
		// this part is for one sample 
	
		mey1[0]  = m11*(-m12*statein[1]  + px[0]         -statein[7]);
	    if (species>1) mey1[0]= m11*(-m12*statein[1]-m13*statein[0]+m14*px[0]+m15*statein[7]+m16*statein[6]);
            mey2[0]  = m21*(-m22*statein[3] - m23*statein[2] + m24*mey1[0] + m25*statein[1] + m26*statein[0]);
            mey3[0]  = m31*(-m32*statein[5] - m33*statein[4] + m34*mey2[0] + m35*statein[3] + m36*statein[2]);
            meout = mey3[0]/megainmax;
			
        
            
       
		
          mey1[1]  = m11*(-m12*mey1[0]  + px[1]         - px[0]);
            if (species>1) mey1[1]= m11*(-m12*mey1[0]-m13*statein[1]+m14*px[n]+m15*px[0]+m16*statein[7]);
            mey2[1]  = m21*(-m22*mey2[0] - m23*statein[3] + m24*mey1[1] + m25*mey1[0] + m26*statein[1]);
            mey3[1]  = m31*(-m32*mey3[0] - m33*statein[5] + m34*mey2[1] + m35*mey2[0] + m36*statein[3]);
            meout = mey3[1]/megainmax;

	*/
statetmp[0]=mey1[0];
statetmp[1]=mey1[1];
statetmp[2]=mey2[0];
statetmp[3]=mey2[1];
statetmp[4]=mey3[0];
statetmp[5]=mey3[1];
statetmp[6]=px[0];
statetmp[7]=px[1];
//totalsim==2

		/* Control-path filter */
        if (n==0) {
		
			phasein=statein[8];
			
			//wbgain=statein[18]; // this is for sample delete this for first
			//tauwb=statein[17];// this is for sample delete this for first
			rin[0]=statein[9];
			rin[1]=statein[10];
			rin[2]=statein[11];
			rin[3]=statein[12];
			rin[4]=statein[13];
			rin[5]=statein[14];
			rin[6]=statein[15];
			rin[7]=statein[16];
			
				  }
	 	
		
		wbout1= WbGammaTone(meout,tdres,centerfreq,n,tauwb,wbgain,wborder,wbphase,phasein,r0,r1,r2,r3,I0,I1,I2,I3,rin);

        

		statetmp[8]=wbphase[0];
        phasein=statetmp[8];
		
		statetmp[9]=r0[0];
		statetmp[10]=r1[0];
		statetmp[11]=r2[0];
		statetmp[12]=r3[0];
		statetmp[13]=I0[0];
		statetmp[14]=I1[0];
		statetmp[15]=I2[0];
		statetmp[16]=I3[0];
	
		for (int r=0; r<8;r++)
		{
			rin[r]=statetmp[r+9];	
		}	
		
	
		 
        wbout  = pow((tauwb/TauWBMax),wborder)*wbout1*10e3*__max(1,cf/5e3);
		
        ohcnonlinout = Boltzman(wbout,ohcasym,12.0,5.0,5.0); /* pass the control signal through OHC Nonlinear Function */

		     if (n==0) {  // this part should be commented for first time?
		
		ohcin[0]=statein[22];
		ohcin[1]=statein[23];
		ohcin[2]=statein[24];
		ohcin[3]=statein[25];
		ohcin[4]=statein[26];
		ohcin[5]=statein[27];
		ohcin[6]=statein[28];
		ohcin[7]=statein[29];
				  }
				  
		ohcout = OhcLowPass(ohcnonlinout,tdres,600,n,1.0,2,Ohc0,Ohc1,Ohc2,Ohc3,Ohcl0,Ohcl1,Ohcl2,Ohcl3,ohcin);/* lowpass filtering after the OHC nonlinearity */
     		
		
	    statetmp[22]=Ohc0[0];
		statetmp[23]=Ohc1[0];
		statetmp[24]=Ohc2[0];
		statetmp[25]=Ohc3[0];
		statetmp[26]=Ohcl0[0];
		statetmp[27]=Ohcl1[0];
		statetmp[28]=Ohcl2[0];
		statetmp[29]=Ohcl3[0];
		
		for (int o=0; o<8;o++)
		{
			ohcin[o]=statetmp[o+22];	
		}
		
		
		tmptauc1 = NLafterohc(ohcout,bmTaumin[0],bmTaumax[0],ohcasym); /* nonlinear function after OHC low-pass filter */
		
		tauc1    = cohc[n]*(tmptauc1-bmTaumin[0])+bmTaumin[0];  /* time -constant for the signal-path C1 filter */
	
		rsigma   = 1/tauc1-1/bmTaumax[0]; /* shift of the location of poles of the C1 filter from the initial positions */

		if (1/tauc1<0.0) mexErrMsgTxt("The poles are in the right-half plane; system is unstable.\n");
       
		tauwb = TauWBMax+(tauc1-bmTaumax[0])*(TauWBMax-TauWBMin)/(bmTaumax[0]-bmTaumin[0]);
		
statetmp[17]=tauwb;


	    wb_gain = gain_groupdelay(tdres,centerfreq,cf,tauwb,grdelay);
		
		  
		grd = grdelay[0]; 
if (n==0) statetmp[20]=grd; 
if (n==1) statetmp[21]=grd;
     if (n==0)
	 {
		 if (grd ==0)  {
			 wbgain=wb_gain;
		 }
		
		// if ((grd+n)<totalstim)
	      //   tmpgain[grd+n] = wb_gain;

       // if (tmpgain[n] == 0)
		//	tmpgain[n] = lasttmpgain;	
		
		if (grd=!0) {
			
			//wbgain = statein[18]; // for afagh>=3
			wbgain=lasttmpgain; // for afagh==1 and 2
		 }
			 lasttmpgain = wbgain;
	 }
	 
	 
	   
		
	if (n==0) 	statetmp[18]=wb_gain; 
	if (n==1)	 	statetmp[19]=wb_gain;


	
        /*====== Signal-path C1 filter ======*/
         

c1filterouttmp  = C1ChirpFilt(meout, tdres, cf, n, bmTaumax[0], rsigma, c1gainout, c1phaseout, C1input11, C1input12, C1input13, C1input21, C1input22, C1input23, C1input31, C1input32, C1input33, C1input41, C1input42, C1input43, C1input51, C1input52, C1input53, C1input61, C1input62, C1input63, C1output11, C1output21, C1output31, C1output41, C1output51, C1output12, C1output22, C1output32, C1output42, C1output52); /* parallel-filter output*/
 
statetmp[60]=c1gainout[0];
statetmp[61]=c1phaseout[0];
statetmp[62]=C1input11[0];
statetmp[63]=C1input12[0];
statetmp[64]=C1input13[0];
statetmp[65]=C1input21[0];
statetmp[66]=C1input22[0];
statetmp[67]=C1input23[0];
statetmp[68]=C1input31[0];
statetmp[69]=C1input32[0];
statetmp[70]=C1input33[0];
statetmp[71]=C1input41[0];
statetmp[72]=C1input42[0];
statetmp[73]=C1input43[0];
statetmp[74]=C1input51[0];
statetmp[75]=C1input52[0];
statetmp[76]=C1input53[0];
statetmp[77]=C1input61[0];
statetmp[78]=C1input62[0];
statetmp[79]=C1input63[0];

statetmp[80]=C1output11[0];
statetmp[81]=C1output21[0];
statetmp[82]=C1output31[0];
statetmp[83]=C1output41[0];
statetmp[84]=C1output51[0];
statetmp[85]=C1output12[0];
statetmp[86]=C1output22[0];
statetmp[87]=C1output32[0];
statetmp[88]=C1output42[0];
statetmp[89]=C1output52[0];
	 
        /*====== Parallel-path C2 filter ======*/
	
c2filterouttmp  = C2ChirpFilt(meout, tdres, cf, n, bmTaumax[0], 1/ratiobm[0], c2gainout, c2phaseout, C2input11, C2input12, C2input13, C2input21, C2input22, C2input23, C2input31, C2input32, C2input33, C2input41, C2input42, C2input43, C2input51, C2input52, C2input53, C2input61, C2input62, C2input63, C2output11, C2output21, C2output31, C2output41, C2output51, C2output12, C2output22, C2output32, C2output42, C2output52); /* parallel-filter output*/


statetmp[30]=c2gainout[0];
statetmp[31]=c2phaseout[0];
statetmp[32]=C2input11[0];
statetmp[33]=C2input12[0];
statetmp[34]=C2input13[0];
statetmp[35]=C2input21[0];
statetmp[36]=C2input22[0];
statetmp[37]=C2input23[0];
statetmp[38]=C2input31[0];
statetmp[39]=C2input32[0];
statetmp[40]=C2input33[0];
statetmp[41]=C2input41[0];
statetmp[42]=C2input42[0];
statetmp[43]=C2input43[0];
statetmp[44]=C2input51[0];
statetmp[45]=C2input52[0];
statetmp[46]=C2input53[0];
statetmp[47]=C2input61[0];
statetmp[48]=C2input62[0];
statetmp[49]=C2input63[0];

statetmp[50]=C2output11[0];
statetmp[51]=C2output21[0];
statetmp[52]=C2output31[0];
statetmp[53]=C2output41[0];
statetmp[54]=C2output51[0];
statetmp[55]=C2output12[0];
statetmp[56]=C2output22[0];
statetmp[57]=C2output32[0];
statetmp[58]=C2output42[0];
statetmp[59]=C2output52[0];

	    /*=== Run the inner hair cell (IHC) section: NL function and then lowpass filtering ===*/

        c1vihctmp  = NLogarithm(cihc*c1filterouttmp,0.1,ihcasym,cf);
	   
		c2vihctmp = -NLogarithm(c2filterouttmp*fabs(c2filterouttmp)*cf/10*cf/2e3,0.2,1.0,cf); /* C2 transduction output */
            
			
     		if (n==0) {  // this part should be commented for first time?
		
		Ihcin[0]=statein[90];
		Ihcin[1]=statein[91];
		Ihcin[2]=statein[92];
		Ihcin[3]=statein[93];
		Ihcin[4]=statein[94];
		Ihcin[5]=statein[95];
		Ihcin[6]=statein[96];
		Ihcin[7]=statein[97];
		Ihcin[8]=statein[98];
		Ihcin[9]=statein[99];
		Ihcin[10]=statein[100];
		Ihcin[11]=statein[101];
		Ihcin[12]=statein[102];
		Ihcin[13]=statein[103];
		Ihcin[14]=statein[104];
		Ihcin[15]=statein[105];
				  }	
			
			
			
			
				 
                ihcouttmp[n] = IhcLowPass(c1vihctmp+c2vihctmp,tdres,3000,n,1.0,7,Ihc0,Ihc1,Ihc2,Ihc3,Ihc4,Ihc5,Ihc6,Ihc7,Ihcl0,Ihcl1,Ihcl2,Ihcl3,Ihcl4,Ihcl5,Ihcl6,Ihcl7,Ihcin);
					ihcout[n]=ihcouttmp[n];
					
	    statetmp[90]=Ihc0[0];
		statetmp[91]=Ihc1[0];
		statetmp[92]=Ihc2[0];
		statetmp[93]=Ihc3[0];
	    statetmp[94]=Ihc4[0];
		statetmp[95]=Ihc5[0];
		statetmp[96]=Ihc6[0];
		statetmp[97]=Ihc7[0];
		
		statetmp[98]=Ihcl0[0];
		statetmp[99]=Ihcl1[0];
		statetmp[100]=Ihcl2[0];
		statetmp[101]=Ihcl3[0];
		statetmp[102]=Ihcl4[0];
		statetmp[103]=Ihcl5[0];
		statetmp[104]=Ihcl6[0];
		statetmp[105]=Ihcl7[0];
		
			for (int I=0; I<16;I++)
		{
			Ihcin[I]=statetmp[I+90];	
		}
		
		
		};  /* End of the loop */
		
	
   
    /* Stretched out the IHC output according to nrep (number of repetitions) */
   
    for(i=0;i<totalstim*nrep;i++)
	{
		ihcouttmp[i] = ihcouttmp[(int) (fmod(i,totalstim))];
  	};   
   	/* Adjust total path delay to IHC output signal */
    if (species==1)
        delay      = delay_cat(cf);
    if (species>1)
    {/*    delay      = delay_human(cf); */
        delay      = delay_cat(cf); /* signal delay changed back to cat function for version 5.2 */
    };
    delaypoint =__max(0,(int) ceil(delay/tdres));    
         
  /*for(i=delaypoint;i<totalstim*nrep;i++)
	{        
		ihcout[i] = ihcouttmp[i - delaypoint];
  	};   
	*/
 
	
state[0]=statetmp[0];
state[1]=statetmp[1];
state[2]=statetmp[2];
state[3]=statetmp[3];
state[4]=statetmp[4];
state[5]=statetmp[5];
state[6]=statetmp[6];
state[7]=statetmp[7];
state[8]=statetmp[8];
state[9]=statetmp[9];
state[10]=statetmp[10];
state[11]=statetmp[11];
state[12]=statetmp[12];
state[13]=statetmp[13];
state[14]=statetmp[14];
state[15]=statetmp[15];
state[16]=statetmp[16];
state[17]=statetmp[17];
state[18]=statetmp[18];
state[19]=statetmp[19];
state[20]=statetmp[20];
state[21]=statetmp[21];
state[22]=statetmp[22];
state[23]=statetmp[23];
state[24]=statetmp[24];
state[25]=statetmp[25];
state[26]=statetmp[26];
state[27]=statetmp[27];
state[28]=statetmp[28];
state[29]=statetmp[29];

state[30]=statetmp[30];
state[31]=statetmp[31];

state[32]=statetmp[32];
state[33]=statetmp[33];
state[34]=statetmp[34];
state[35]=statetmp[35];
state[36]=statetmp[36];
state[37]=statetmp[37];
state[38]=statetmp[38];
state[39]=statetmp[39];
state[40]=statetmp[40];
state[41]=statetmp[41];
state[42]=statetmp[42];
state[43]=statetmp[43];
state[44]=statetmp[44];
state[45]=statetmp[45];
state[46]=statetmp[46];
state[47]=statetmp[47];
state[48]=statetmp[48];
state[49]=statetmp[49];

state[50]=statetmp[50];
state[51]=statetmp[51];
state[52]=statetmp[52];
state[53]=statetmp[53];
state[54]=statetmp[54];
state[55]=statetmp[55];
state[56]=statetmp[56];
state[57]=statetmp[57];
state[58]=statetmp[58];
state[59]=statetmp[59];


state[60]=statetmp[60];
state[61]=statetmp[61];

state[62]=statetmp[62];
state[63]=statetmp[63];
state[64]=statetmp[64];
state[65]=statetmp[65];
state[66]=statetmp[66];
state[67]=statetmp[67];
state[68]=statetmp[68];
state[69]=statetmp[69];
state[70]=statetmp[70];
state[71]=statetmp[71];
state[72]=statetmp[72];
state[73]=statetmp[73];
state[74]=statetmp[74];
state[75]=statetmp[75];
state[76]=statetmp[76];
state[77]=statetmp[77];
state[78]=statetmp[78];
state[79]=statetmp[79];

state[80]=statetmp[80];
state[81]=statetmp[81];
state[82]=statetmp[82];
state[83]=statetmp[83];
state[84]=statetmp[84];
state[85]=statetmp[85];
state[86]=statetmp[86];
state[87]=statetmp[87];
state[88]=statetmp[88];
state[89]=statetmp[89];

state[90]=statetmp[90];
state[91]=statetmp[91];
state[92]=statetmp[92];
state[93]=statetmp[93];
state[94]=statetmp[94];
state[95]=statetmp[95];
state[96]=statetmp[96];
state[97]=statetmp[97];
state[98]=statetmp[98];
state[99]=statetmp[99];
state[100]=statetmp[100];
state[101]=statetmp[101];
state[102]=statetmp[102];
state[103]=statetmp[103];
state[104]=statetmp[104];
state[105]=statetmp[105];
    /* Freeing dynamic memory allocated earlier */

    mxFree(ihcouttmp); mxFree(statetmp); 
    mxFree(mey1); mxFree(mey2); mxFree(mey3);	
    mxFree(tmpgain); mxFree(rin);
	mxFree(ohcin);	mxFree(Ihcin);

} /* End of the SingleAN function */
/* -------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------- */
/** Get TauMax, TauMin for the tuning filter. The TauMax is determined by the bandwidth/Q10
    of the tuning filter at low level. The TauMin is determined by the gain change between high
    and low level */

double Get_tauwb(double cf, int species, int order, double *taumax,double *taumin)
{
  double Q10,bw,gain,ratio;
    
  if(species==1) gain = 52.0/2.0*(tanh(2.2*log10(cf/0.6e3)+0.15)+1.0); /* for cat */
  if(species>1) gain = 52.0/2.0*(tanh(2.2*log10(cf/0.6e3)+0.15)+1.0); /* for human */
  /*gain = 52/2*(tanh(2.2*log10(cf/1e3)+0.15)+1);*/ /* older values */

  if(gain>60.0) gain = 60.0;  
  if(gain<15.0) gain = 15.0;
   
  ratio = pow(10,(-gain/(20.0*order)));       /* ratio of TauMin/TauMax according to the gain, order */
  if (species==1) /* cat Q10 values */
  {
    Q10 = pow(10,0.4708*log10(cf/1e3)+0.4664);
  }
  if (species==2) /* human Q10 values from Shera et al. (PNAS 2002) */
  {
    Q10 = pow((cf/1000),0.3)*12.7*0.505+0.2085;
  }
  if (species==3) /* human Q10 values from Glasberg & Moore (Hear. Res. 1990) */
  {
    Q10 = cf/24.7/(4.37*(cf/1000)+1)*0.505+0.2085;
  }
  bw     = cf/Q10;
  taumax[0] = 2.0/(TWOPI*bw);
   
  taumin[0]   = taumax[0]*ratio;
  
  return 0;
}
/* -------------------------------------------------------------------------------------------- */
double Get_taubm(double cf, int species, double taumax,double *bmTaumax,double *bmTaumin, double *ratio)
{
  double gain,factor,bwfactor;
    
  if(species==1) gain = 52.0/2.0*(tanh(2.2*log10(cf/0.6e3)+0.15)+1.0); /* for cat */
  if(species>1) gain = 52.0/2.0*(tanh(2.2*log10(cf/0.6e3)+0.15)+1.0); /* for human */
  /*gain = 52/2*(tanh(2.2*log10(cf/1e3)+0.15)+1);*/ /* older values */

 
  if(gain>60.0) gain = 60.0;  
  if(gain<15.0) gain = 15.0;

  bwfactor = 0.7;
  factor   = 2.5;

  ratio[0]  = pow(10,(-gain/(20.0*factor))); 

  bmTaumax[0] = taumax/bwfactor;
  bmTaumin[0] = bmTaumax[0]*ratio[0];     
  return 0;
}
/* -------------------------------------------------------------------------------------------- */
/** Pass the signal through the signal-path C1 Tenth Order Nonlinear Chirp-Gammatone Filter */

double C1ChirpFilt(double x, double tdres,double cf, int n, double taumax, double rsigma, double*c1gainout, double*c1phaseout, double*C1input11, double*C1input12, double*C1input13, double*C1input21, double*C1input22, double*C1input23, double*C1input31, double*C1input32, double*C1input33, double*C1input41, double*C1input42, double*C1input43, double*C1input51, double*C1input52, double*C1input53, double*C1input61, double*C1input62, double*C1input63, double*C1output11, double*C1output21, double*C1output31, double*C1output41, double*C1output51, double*C1output12, double*C1output22, double*C1output32, double*C1output42, double*C1output52)
{
    static double C1gain_norm, C1initphase; 
    static double C1input[12][4], C1output[12][4];

    double ipw, ipb, rpa, pzero, rzero;
	double sigma0,fs_bilinear,CF,norm_gain,phase,c1filterout;
	int i,r,order_of_pole,half_order_pole,order_of_zero;
	double temp, dy, preal, pimg;

	COMPLEX p[11]; 
	
	/* Defining initial locations of the poles and zeros */
	/*======== setup the locations of poles and zeros =======*/
	  sigma0 = 1/taumax;
	  ipw    = 1.01*cf*TWOPI-50;
	  ipb    = 0.2343*TWOPI*cf-1104;
	  rpa    = pow(10, log10(cf)*0.9 + 0.55)+ 2000;
	  pzero  = pow(10,log10(cf)*0.7+1.6)+500;

	/*===============================================================*/     
         
     order_of_pole    = 10;             
     half_order_pole  = order_of_pole/2;
     order_of_zero    = half_order_pole;

	 fs_bilinear = TWOPI*cf/tan(TWOPI*cf*tdres/2);
     rzero       = -pzero;
	 CF          = TWOPI*cf;
   
   if (n==0)
   {		  
	p[1].x = -sigma0;     

    p[1].y = ipw;

	p[5].x = p[1].x - rpa; p[5].y = p[1].y - ipb;

    p[3].x = (p[1].x + p[5].x) * 0.5; p[3].y = (p[1].y + p[5].y) * 0.5;

    p[2]   = compconj(p[1]);    p[4] = compconj(p[3]); p[6] = compconj(p[5]);

    p[7]   = p[1]; p[8] = p[2]; p[9] = p[5]; p[10]= p[6];

	   C1initphase = 0.0;
       for (i=1;i<=half_order_pole;i++)          
	   {
           preal     = p[i*2-1].x;
		   pimg      = p[i*2-1].y;
	       C1initphase = C1initphase + atan(CF/(-rzero))-atan((CF-pimg)/(-preal))-atan((CF+pimg)/(-preal));
	   };

	/*===================== Initialize C1input & C1output =====================*/

      for (i=1;i<=(half_order_pole+1);i++)          
      {
		   C1input[i][3] = 0; 
		   C1input[i][2] = 0; 
		   C1input[i][1] = 0;
		   C1output[i][3] = 0; 
		   C1output[i][2] = 0; 
		   C1output[i][1] = 0;
      }

	/*===================== normalize the gain =====================*/
    
      C1gain_norm = 1.0;
      for (r=1; r<=order_of_pole; r++)
		   C1gain_norm = C1gain_norm*(pow((CF - p[r].y),2) + p[r].x*p[r].x);
      
   };
     
    norm_gain= sqrt(C1gain_norm)/pow(sqrt(CF*CF+rzero*rzero),order_of_zero);
	
	p[1].x = -sigma0 - rsigma;

	if (p[1].x>0.0) mexErrMsgTxt("The system becomes unstable.\n");
	
	p[1].y = ipw;

	p[5].x = p[1].x - rpa; p[5].y = p[1].y - ipb;

    p[3].x = (p[1].x + p[5].x) * 0.5; p[3].y = (p[1].y + p[5].y) * 0.5;

    p[2] = compconj(p[1]); p[4] = compconj(p[3]); p[6] = compconj(p[5]);

    p[7] = p[1]; p[8] = p[2]; p[9] = p[5]; p[10]= p[6];

    phase = 0.0;
    for (i=1;i<=half_order_pole;i++)          
    {
           preal = p[i*2-1].x;
		   pimg  = p[i*2-1].y;
	       phase = phase-atan((CF-pimg)/(-preal))-atan((CF+pimg)/(-preal));
	};

	rzero = -CF/tan((C1initphase-phase)/order_of_zero);

    if (rzero>0.0) mexErrMsgTxt("The zeros are in the right-half plane.\n");
	 
   /*%==================================================  */
	/*each loop below is for a pair of poles and one zero */
   /*%      time loop begins here                         */
   /*%==================================================  */
 
       C1input[1][3]=C1input[1][2]; 
	   C1input[1][2]=C1input[1][1]; 
	   C1input[1][1]= x;

       for (i=1;i<=half_order_pole;i++)          
       {
           preal = p[i*2-1].x;
		   pimg  = p[i*2-1].y;
		  	   
           temp  = pow((fs_bilinear-preal),2)+ pow(pimg,2);
		   

           /*dy = (input[i][1] + (1-(fs_bilinear+rzero)/(fs_bilinear-rzero))*input[i][2]
                                 - (fs_bilinear+rzero)/(fs_bilinear-rzero)*input[i][3] );
           dy = dy+2*output[i][1]*(fs_bilinear*fs_bilinear-preal*preal-pimg*pimg);

           dy = dy-output[i][2]*((fs_bilinear+preal)*(fs_bilinear+preal)+pimg*pimg);*/
		   
	       dy = C1input[i][1]*(fs_bilinear-rzero) - 2*rzero*C1input[i][2] - (fs_bilinear+rzero)*C1input[i][3]
                 +2*C1output[i][1]*(fs_bilinear*fs_bilinear-preal*preal-pimg*pimg)
			     -C1output[i][2]*((fs_bilinear+preal)*(fs_bilinear+preal)+pimg*pimg);

		   dy = dy/temp;

		   C1input[i+1][3] = C1output[i][2]; 
		   C1input[i+1][2] = C1output[i][1]; 
		   C1input[i+1][1] = dy;

		   C1output[i][2] = C1output[i][1]; 
		   C1output[i][1] = dy;
       }

	   dy = C1output[half_order_pole][1]*norm_gain;  /* don't forget the gain term */
	   c1filterout= dy/4.0;   /* signal path output is divided by 4 to give correct C1 filter gain */
	  c1gainout[0]=C1gain_norm;
	  c1phaseout[0]=C1initphase;
	  
	  C1input11[0]=C1input[1][1];
	  C1input12[0]=C1input[1][2];
	  C1input13[0]=C1input[1][3];
	  C1input21[0]=C1input[2][1];
	  C1input22[0]=C1input[2][2];
	  C1input23[0]=C1input[2][3];
	  C1input31[0]=C1input[3][1];
	  C1input32[0]=C1input[3][2];
	  C1input33[0]=C1input[3][3];
	  C1input41[0]=C1input[4][1];
	  C1input42[0]=C1input[4][2];
	  C1input43[0]=C1input[4][3];
	  C1input51[0]=C1input[5][1];
	  C1input52[0]=C1input[5][2];
	  C1input53[0]=C1input[5][3];
	  C1input61[0]=C1input[6][1];
	  C1input62[0]=C1input[6][2];
	  C1input63[0]=C1input[6][3];
	  C1output11[0]=C1output[1][1];
	  C1output21[0]=C1output[2][1];
	  C1output31[0]=C1output[3][1];
	  C1output41[0]=C1output[4][1];
	  C1output51[0]=C1output[5][1];
	  C1output12[0]=C1output[1][2];
	  C1output22[0]=C1output[2][2];
	  C1output32[0]=C1output[3][2];
	  C1output42[0]=C1output[4][2];
	  C1output52[0]=C1output[5][2];	                   
     return (c1filterout);
}  

/* -------------------------------------------------------------------------------------------- */
/** Parallelpath C2 filter: same as the signal-path C1 filter with the OHC completely impaired */

double C2ChirpFilt(double xx, double tdres,double cf, int n, double taumax, double fcohc, double*c2gainout, double*c2phaseout, double*C2input11, double*C2input12, double*C2input13, double*C2input21, double*C2input22, double*C2input23, double*C2input31, double*C2input32, double*C2input33, double*C2input41, double*C2input42, double*C2input43, double*C2input51, double*C2input52, double*C2input53, double*C2input61, double*C2input62, double*C2input63, double*C2output11, double*C2output21, double*C2output31, double*C2output41, double*C2output51, double*C2output12, double*C2output22, double*C2output32, double*C2output42, double*C2output52)
{
	
	static double C2gain_norm, C2initphase;
    static double C2input[12][4];  static double C2output[12][4];
   
	double ipw, ipb, rpa, pzero, rzero;

	double sigma0,fs_bilinear,CF,norm_gain,phase,c2filterout;
	int    i,r,order_of_pole,half_order_pole,order_of_zero;
	double temp, dy, preal, pimg;

	COMPLEX p[11]; 	
    
    /*================ setup the locations of poles and zeros =======*/

	  sigma0 = 1/taumax;
	  ipw    = 1.01*cf*TWOPI-50;
      ipb    = 0.2343*TWOPI*cf-1104;
	  rpa    = pow(10, log10(cf)*0.9 + 0.55)+ 2000;
	  pzero  = pow(10,log10(cf)*0.7+1.6)+500;
	/*===============================================================*/     
         
     order_of_pole    = 10;             
     half_order_pole  = order_of_pole/2;
     order_of_zero    = half_order_pole;

	 fs_bilinear = TWOPI*cf/tan(TWOPI*cf*tdres/2);
     rzero       = -pzero;
	 CF          = TWOPI*cf;
   	    
    if (n==0)
    {		  
	p[1].x = -sigma0;     

    p[1].y = ipw;

	p[5].x = p[1].x - rpa; p[5].y = p[1].y - ipb;

    p[3].x = (p[1].x + p[5].x) * 0.5; p[3].y = (p[1].y + p[5].y) * 0.5;

    p[2] = compconj(p[1]); p[4] = compconj(p[3]); p[6] = compconj(p[5]);

    p[7] = p[1]; p[8] = p[2]; p[9] = p[5]; p[10]= p[6];

	   C2initphase = 0.0;
       for (i=1;i<=half_order_pole;i++)         
	   {
           preal     = p[i*2-1].x;
		   pimg      = p[i*2-1].y;
	       C2initphase = C2initphase + atan(CF/(-rzero))-atan((CF-pimg)/(-preal))-atan((CF+pimg)/(-preal));
	   };

	/*===================== Initialize C2input & C2output =====================*/

      for (i=1;i<=(half_order_pole+1);i++)          
      {
		   C2input[i][3] = 0; 
		   C2input[i][2] = 0; 
		   C2input[i][1] = 0;
		   C2output[i][3] = 0; 
		   C2output[i][2] = 0; 
		   C2output[i][1] = 0;
      }
    
    /*===================== normalize the gain =====================*/
    
     C2gain_norm = 1.0;
     for (r=1; r<=order_of_pole; r++)
		   C2gain_norm = C2gain_norm*(pow((CF - p[r].y),2) + p[r].x*p[r].x);
    };
     
    norm_gain= sqrt(C2gain_norm)/pow(sqrt(CF*CF+rzero*rzero),order_of_zero);
    
	p[1].x = -sigma0*fcohc;

	if (p[1].x>0.0) mexErrMsgTxt("The system becomes unstable.\n");
	
	p[1].y = ipw;

	p[5].x = p[1].x - rpa; p[5].y = p[1].y - ipb;

    p[3].x = (p[1].x + p[5].x) * 0.5; p[3].y = (p[1].y + p[5].y) * 0.5;

    p[2] = compconj(p[1]); p[4] = compconj(p[3]); p[6] = compconj(p[5]);

    p[7] = p[1]; p[8] = p[2]; p[9] = p[5]; p[10]= p[6];

    phase = 0.0;
    for (i=1;i<=half_order_pole;i++)          
    {
           preal = p[i*2-1].x;
		   pimg  = p[i*2-1].y;
	       phase = phase-atan((CF-pimg)/(-preal))-atan((CF+pimg)/(-preal));
	};

	rzero = -CF/tan((C2initphase-phase)/order_of_zero);	
    if (rzero>0.0) mexErrMsgTxt("The zeros are in the right-hand plane.\n");
   /*%==================================================  */
   /*%      time loop begins here                         */
   /*%==================================================  */

       C2input[1][3]=C2input[1][2]; 
	   C2input[1][2]=C2input[1][1]; 
	   C2input[1][1]= xx;

      for (i=1;i<=half_order_pole;i++)          
      {
           preal = p[i*2-1].x;
		   pimg  = p[i*2-1].y;
		  	   
           temp  = pow((fs_bilinear-preal),2)+ pow(pimg,2);
		   
           /*dy = (input[i][1] + (1-(fs_bilinear+rzero)/(fs_bilinear-rzero))*input[i][2]
                                 - (fs_bilinear+rzero)/(fs_bilinear-rzero)*input[i][3] );
           dy = dy+2*output[i][1]*(fs_bilinear*fs_bilinear-preal*preal-pimg*pimg);

           dy = dy-output[i][2]*((fs_bilinear+preal)*(fs_bilinear+preal)+pimg*pimg);*/
		   
	      dy = C2input[i][1]*(fs_bilinear-rzero) - 2*rzero*C2input[i][2] - (fs_bilinear+rzero)*C2input[i][3]
                 +2*C2output[i][1]*(fs_bilinear*fs_bilinear-preal*preal-pimg*pimg)
			     -C2output[i][2]*((fs_bilinear+preal)*(fs_bilinear+preal)+pimg*pimg);

		   dy = dy/temp;

		   C2input[i+1][3] = C2output[i][2]; 
		   C2input[i+1][2] = C2output[i][1]; 
		   C2input[i+1][1] = dy;

		   C2output[i][2] = C2output[i][1]; 
		   C2output[i][1] = dy;

       };

	  dy = C2output[half_order_pole][1]*norm_gain;
	  c2filterout= dy/4.0;
	  
	  c2gainout[0]=C2gain_norm;
	  c2phaseout[0]=C2initphase;
	  
	  C2input11[0]=C2input[1][1];
	  C2input12[0]=C2input[1][2];
	  C2input13[0]=C2input[1][3];
	  C2input21[0]=C2input[2][1];
	  C2input22[0]=C2input[2][2];
	  C2input23[0]=C2input[2][3];
	  C2input31[0]=C2input[3][1];
	  C2input32[0]=C2input[3][2];
	  C2input33[0]=C2input[3][3];
	  C2input41[0]=C2input[4][1];
	  C2input42[0]=C2input[4][2];
	  C2input43[0]=C2input[4][3];
	  C2input51[0]=C2input[5][1];
	  C2input52[0]=C2input[5][2];
	  C2input53[0]=C2input[5][3];
	  C2input61[0]=C2input[6][1];
	  C2input62[0]=C2input[6][2];
	  C2input63[0]=C2input[6][3];
	  C2output11[0]=C2output[1][1];
	  C2output21[0]=C2output[2][1];
	  C2output31[0]=C2output[3][1];
	  C2output41[0]=C2output[4][1];
	  C2output51[0]=C2output[5][1];
	  C2output12[0]=C2output[1][2];
	  C2output22[0]=C2output[2][2];
	  C2output32[0]=C2output[3][2];
	  C2output42[0]=C2output[4][2];
	  C2output52[0]=C2output[5][2];
	  return (c2filterout); 
}   

/* -------------------------------------------------------------------------------------------- */
/** Pass the signal through the Control path Third Order Nonlinear Gammatone Filter */

double WbGammaTone(double x,double tdres,double centerfreq, int n, double tau,double gain,int order, double *wbphase, double phasein, double*r0,double*r1, double*r2, double*r3, double*I0 ,double*I1, double*I2, double*I3,double* rin )
{
   
  COMPLEX wbgtf[4], wbgtfl[4];

  double delta_phase,dtmp,c1LP,c2LP,out;
  int i,j;
  
  
      
            wbgtfl[0].x = rin[0];
			wbgtfl[0].y = rin[4];
			
			wbgtfl[1].x = rin[1];
			wbgtfl[1].y = rin[5];
			
			wbgtfl[2].x = rin[2];
			wbgtfl[2].y = rin[6];
			
			wbgtfl[3].x = rin[3];
			wbgtfl[3].y = rin[7];
			
		
		
		

		
		
		
  delta_phase = -TWOPI*centerfreq*tdres;
  wbphase[0] = phasein+delta_phase;
 
  
  dtmp = tau*2.0/tdres;
  c1LP = (dtmp-1)/(dtmp+1);
  c2LP = 1.0/(dtmp+1);
  wbgtf[0] = compmult(x,compexp(wbphase[0]));                 /* FREQUENCY SHIFT */
  
  for(j = 1; j <= order; j++)                              /* IIR Bilinear transformation LPF */
  wbgtf[j] = comp2sum(compmult(c2LP*gain,comp2sum(wbgtf[j-1],wbgtfl[j-1])),
      compmult(c1LP,wbgtfl[j]));
	  
  out = REAL(compprod(compexp(-wbphase[0]), wbgtf[order])); /* FREQ SHIFT BACK UP */
  
  for(i=0; i<=order;i++) wbgtfl[i] = wbgtf[i];
  r0[0]=REAL(wbgtfl[0]);
  r1[0]=REAL(wbgtfl[1]);
  r2[0]=REAL(wbgtfl[2]);
  r3[0]=REAL(wbgtfl[3]);
  I0[0]=IMAG(wbgtfl[0]);
  I1[0]=IMAG(wbgtfl[1]);
  I2[0]=IMAG(wbgtfl[2]);
  I3[0]=IMAG(wbgtfl[3]);
  return(out);
}

/* -------------------------------------------------------------------------------------------- */
/** Calculate the gain and group delay for the Control path Filter */

double gain_groupdelay(double tdres,double centerfreq, double cf, double tau,int *grdelay)
{ 
  double tmpcos,dtmp2,c1LP,c2LP,tmp1,tmp2,wb_gain;

  tmpcos = cos(TWOPI*(centerfreq-cf)*tdres);
  dtmp2 = tau*2.0/tdres;
  c1LP = (dtmp2-1)/(dtmp2+1);
  c2LP = 1.0/(dtmp2+1);
  tmp1 = 1+c1LP*c1LP-2*c1LP*tmpcos;
  tmp2 = 2*c2LP*c2LP*(1+tmpcos);
  
  wb_gain = pow(tmp1/tmp2, 1.0/2.0);
  
  grdelay[0] = (int)floor((0.5-(c1LP*c1LP-c1LP*tmpcos)/(1+c1LP*c1LP-2*c1LP*tmpcos)));

  return(wb_gain);
}
/* -------------------------------------------------------------------------------------------- */
/** Calculate the delay (basilar membrane, synapse, etc. for cat) */
double delay_cat(double cf)
{  
  double A0,A1,x,delay;

  A0    = 3.0;  
  A1    = 12.5;
  x     = 11.9 * log10(0.80 + cf / 456.0);      /* cat mapping */
  delay = A0 * exp( -x/A1 ) * 1e-3;
  
  return(delay);
}

/* Calculate the delay (basilar membrane, synapse, etc.) for human, based
        on Harte et al. (JASA 2009) */
double delay_human(double cf) 
{  
  double A,B,delay;

  A    = -0.37;  
  B    = 11.09/2;
  delay = B * pow(cf * 1e-3,A)*1e-3;
  
  return(delay);
}

/* -------------------------------------------------------------------------------------------- */
/* Get the output of the OHC Nonlinear Function (Boltzman Function) */

double Boltzman(double x, double asym, double s0, double s1, double x1)
  {
	double shift,x0,out1,out;

    shift = 1.0/(1.0+asym);  /* asym is the ratio of positive Max to negative Max*/
    x0    = s0*log((1.0/shift-1)/(1+exp(x1/s1)));
	    
    out1 = 1.0/(1.0+exp(-(x-x0)/s0)*(1.0+exp(-(x-x1)/s1)))-shift;
	out = out1/(1-shift);

    return(out);
  }  /* output of the nonlinear function, the output is normalized with maximum value of 1 */
  
/* -------------------------------------------------------------------------------------------- */
/* Get the output of the OHC Low Pass Filter in the Control path */

double OhcLowPass(double x,double tdres,double Fc, int n,double gain,int order, double* Ohc0,double* Ohc1, double* Ohc2, double* Ohc3, double* Ohcl0 ,double* Ohcl1, double* Ohcl2, double* Ohcl3, double* ohcin)
{ 
  double ohc[4],ohcl[4];

  double c,c1LP,c2LP;
  int i,j;
/*
  if (n==0)
  {
      for(i=0; i<(order+1);i++)
      {
          ohc[i] = 0;
          ohcl[i] = 0;
      }
  }    
  */
  
  ohc[0]=ohcin[0];
  ohc[1]=ohcin[1];
  ohc[2]=ohcin[2];
  ohc[3]=ohcin[3];
  ohcl[0]=ohcin[4];
  ohcl[1]=ohcin[5];
  ohcl[2]=ohcin[6];
  ohcl[3]=ohcin[7];  
  
  
  c = 2.0/tdres;
  c1LP = ( c - TWOPI*Fc ) / ( c + TWOPI*Fc );
  c2LP = TWOPI*Fc / (TWOPI*Fc + c);
  
  ohc[0] = x*gain;
  for(i=0; i<order;i++)
    ohc[i+1] = c1LP*ohcl[i+1] + c2LP*(ohc[i]+ohcl[i]);
  for(j=0; j<=order;j++) ohcl[j] = ohc[j];
  
  Ohc0[0]= ohc[0];
  Ohc1[0]= ohc[1];
  Ohc2[0]= ohc[2];
  Ohc3[0]= ohc[3];
  
  Ohcl0[0]=ohcl[0];
  Ohcl1[0]=ohcl[1];
  Ohcl2[0]=ohcl[2];
  Ohcl3[0]=ohcl[3];
  
  return(ohc[order]);
}
/* -------------------------------------------------------------------------------------------- */
/* Get the output of the IHC Low Pass Filter  */

double IhcLowPass(double x,double tdres,double Fc, int n,double gain,int order, double* Ihc0,double* Ihc1, double* Ihc2, double* Ihc3, double* Ihc4, double* Ihc5, double* Ihc6, double* Ihc7, double* Ihcl0 ,double* Ihcl1, double* Ihcl2, double* Ihcl3, double* Ihcl4 ,double* Ihcl5, double* Ihcl6, double* Ihcl7,double* Ihcin)
             
{
   double Ihc[8],Ihcl[8];
  
  double C,c1LP,c2LP;
  int i,j;
/*
  if (n==0)
  {
      for(i=0; i<(order+1);i++)
      {
          ihc[i] = 0;
          ihcl[i] = 0;
      }
  }     
  */
  Ihc[0]=Ihcin[0];
  Ihc[1]=Ihcin[1];
  Ihc[2]=Ihcin[2];
  Ihc[3]=Ihcin[3];
  Ihc[4]=Ihcin[4];
  Ihc[5]=Ihcin[5];
  Ihc[6]=Ihcin[6];
  Ihc[7]=Ihcin[7];
  
  Ihcl[0]=Ihcin[8];
  Ihcl[1]=Ihcin[9];
  Ihcl[2]=Ihcin[10];
  Ihcl[3]=Ihcin[11];
  Ihcl[4]=Ihcin[12];
  Ihcl[5]=Ihcin[13];
  Ihcl[6]=Ihcin[14];
  Ihcl[7]=Ihcin[15];
  
  
  C = 2.0/tdres;
  c1LP = ( C - TWOPI*Fc ) / ( C + TWOPI*Fc );
  c2LP = TWOPI*Fc / (TWOPI*Fc + C);
  
  Ihc[0] = x*gain;
  for(i=0; i<order;i++)
    Ihc[i+1] = c1LP*Ihcl[i+1] + c2LP*(Ihc[i]+Ihcl[i]);
  for(j=0; j<=order;j++) Ihcl[j] = Ihc[j];
  
  Ihc0[0]= Ihc[0];
  Ihc1[0]= Ihc[1];
  Ihc2[0]= Ihc[2];
  Ihc3[0]= Ihc[3];
  Ihc4[0]= Ihc[4];
  Ihc5[0]= Ihc[5];
  Ihc6[0]= Ihc[6];
  Ihc7[0]= Ihc[7];
  
  Ihcl0[0]=Ihcl[0];
  Ihcl1[0]=Ihcl[1];
  Ihcl2[0]=Ihcl[2];
  Ihcl3[0]=Ihcl[3];
  Ihcl4[0]=Ihcl[4];
  Ihcl5[0]=Ihcl[5];
  Ihcl6[0]=Ihcl[6];
  Ihcl7[0]=Ihcl[7];
  
  return(Ihc[order]);
}
/* -------------------------------------------------------------------------------------------- */
/* Get the output of the Control path using Nonlinear Function after OHC */

double NLafterohc(double x,double taumin, double taumax, double asym)
{    
	double R,dc,R1,s0,x1,out,minR;

	minR = 0.05;
    R  = taumin/taumax;
    
	if(R<minR) minR = 0.5*R;
    else       minR = minR;
    
    dc = (asym-1)/(asym+1.0)/2.0-minR;
    R1 = R-minR;

    /* This is for new nonlinearity */
    s0 = -dc/log(R1/(1-minR));
	
    x1  = fabs(x);
    out = taumax*(minR+(1.0-minR)*exp(-x1/s0));
	if (out<taumin) out = taumin; 
    if (out>taumax) out = taumax;
    return(out);
}
/* -------------------------------------------------------------------------------------------- */
/* Get the output of the IHC Nonlinear Function (Logarithmic Transduction Functions) */

double NLogarithm(double x, double slope, double asym, double cf)
{
	double corner,strength,xx,splx,asym_t;
	    
    corner    = 80; 
    strength  = 20.0e6/pow(10,corner/20);
            
    xx = log(1.0+strength*fabs(x))*slope;
    
    if(x<0)
	{
		splx   = 20*log10(-x/20e-6);
		asym_t = asym -(asym-1)/(1+exp(splx/5.0));
		xx = -1/asym_t*xx;
	};   
    return(xx);
}
/* -------------------------------------------------------------------------------------------- */