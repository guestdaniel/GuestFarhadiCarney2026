/* This is the 2019 version of the code for auditory periphery model from the Carney, Bruce and Zilany labs.
 * 
 * This release implement the efferents system and is updated  by Afagh Farhadi at University of Rochester with input from the Carney lab. 
 * This model input one sample of the stimuli in each run.
 *   
 *
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
	
	double *px, cf, tdres, reptime, cohc, cihc, *statein, noisein;
	int    nrep, pxbins, lp, cp, s,totalstim, species, lsr;
    mwSize outsize[2], outsize2[2];
    
	double *noiseintmp, *pxtmp, *cftmp, *nreptmp, *tdrestmp, *reptimetmp, *cohctmp, *cihctmp, *speciestmp, *lsrtmp, *stateintmp;
    double *ihcout, *state;
   
	void   IHCAN(double *, double, int, double, int, double , double, int, double *, double *, double *, double);
	
	/* Check for proper number of arguments */
	
	if (nrhs != 10) 
	{
		mexErrMsgTxt("model_IHC requires 10 input arguments.");
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
	stateintmp	= mxGetPr(prhs[8]);
	noiseintmp   = mxGetPr(prhs[9]);
//	lsrtmp	= mxGetPr(prhs[9]);
	/* Check individual input arguments */

	pxbins = (int) mxGetN(prhs[0]);
	//lsr = (int) mxGetN(prhs[9]);
	

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

	
    noisein=noiseintmp[0];
	cihc = cihctmp[0]; /* impairment in the IHC  */
	cohc = cohctmp[0]; /* impairment in the OHC  */
	
	if ((cihc<0)|(cihc>1))
	{
		mexPrintf("cihc (= %1.1f) must be between 0 and 1\n",cihc);
		mexErrMsgTxt("\n");
	}
	if ((cohc<0)|(cohc>1))
	{
		mexPrintf("cohc (= %1.1f) must be between 0 and 1\n",cohc);
		mexErrMsgTxt("\n");
	}
	
   
	/* Calculate number of samples for total repetition time */

	/*totalstim = (int)floor((reptime*1e3)/(tdres*1e3)); */ /*older definition*/
    totalstim = (int)floor(reptime/tdres+0.5);
    px = (double*)mxCalloc(totalstim,sizeof(double)); 
	statein =(double*)mxCalloc(112,sizeof(double)); 
	/* Put stimulus waveform into pressure waveform */

	for (lp=0; lp<pxbins; lp++)
	{
		px[lp] = pxtmp[lp];
		   
	}
	/* Create an array for the return argument */
	for (s=0; s<112; s++)
	{
		statein[s]=stateintmp[s];
	}
	
    outsize[0] = 1;
	outsize[1] = totalstim*nrep;
    outsize2[0] = 1;
	outsize2[1] = 112;   // change this when adding states
	plhs[0] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericArray(2, outsize2, mxDOUBLE_CLASS, mxREAL);
	
	/* Assign pointers to the outputs */
	
	ihcout = mxGetPr(plhs[0]);
	state  = mxGetPr(plhs[1]);	
	
	/* run the model */

	IHCAN(px,cf,nrep,tdres,totalstim,cohc,cihc,species,ihcout,state,statein,noisein);

 mxFree(px);

}

void IHCAN(double *px, double cf, int nrep, double tdres, int totalstim,
                double cohc, double cihc, int species, double *ihcout,  double *state, double *statein, double noisein)
{	
    
    /*variables for middle-ear model */
	double megainmax;
    double mey1, mey2, mey3, meout,c1filterouttmp,c2filterouttmp,c1vihctmp,c2vihctmp, phasein, *rin, *ohcin, *ohcstateout, *c2in, *c2out, *c1in, *c1out,*Ihcin,*Ihcstateout, *rout;
    double fp,C,m11,m12,m13,m14,m15,m16,m21,m22,m23,m24,m25,m26,m31,m32,m33,m34,m35,m36;

	/*variables for the signal-path, control-path and onward */
	double *ihcouttmp,*tmpgain, *statetmp;
	int    grd;

    double bmplace,centerfreq,gain,taubm,ratiowb,bmTaubm,fcohc,TauWBMax,TauWBMin,tauwb;
    double Taumin[1],Taumax[1],bmTaumin[1],bmTaumax[1],ratiobm[1],lasttmpgain,wbgain,ohcasym,ihcasym,delay,wbphase[1],c2gainin, c2phasein;
	double c1gainin, c1phasein;
	int    i,n,delaypoint,grdelay[1],bmorder,wborder;

	double wbout1,wbout,ohcnonlinout,ohcout,tmptauc1,tauc1,rsigma,wb_gain;
            
    /* Declarations of the functions used in the program */
    double C1ChirpFilt(double, double,double, int, double, double, double , double, double *,double *);
	double C2ChirpFilt(double, double,double, int, double, double, double , double, double *,double *);
    double WbGammaTone(double, double, double, int, double, double, int, double *, double, double *,double *);

    double Get_tauwb(double, int, int, double *, double *);
	double Get_taubm(double, int, double, double *, double *, double *);
    double gain_groupdelay(double, double, double, double, int *);
    double delay_cat(double cf);
    double delay_human(double cf);

    double OhcLowPass(double, double, double, int, double, int, double *,double *);
    double IhcLowPass(double, double, double, int, double, int, double *,double *);
	double Boltzman(double, double, double, double, double);
    double NLafterohc(double, double, double, double);
	double ControlSignal(double, double, double, double, double);

    double NLogarithm(double, double, double, double);
    
    /* Allocate dynamic memory for the temporary variables */
	ihcouttmp  = (double*)mxCalloc(totalstim*nrep,sizeof(double));
	statetmp   = (double*)mxCalloc(112,sizeof(double));
	//statetmp2   = (double*)mxCalloc(totalstim,sizeof(double));
	//mey1 = (double*)mxCalloc(totalstim,sizeof(double));
	//mey2 = (double*)mxCalloc(totalstim,sizeof(double));
	//mey3 = (double*)mxCalloc(totalstim,sizeof(double));
	rin  = (double*)mxCalloc(8,sizeof(double));
    rout = (double*)mxCalloc(8,sizeof(double));
	ohcin = (double*)mxCalloc(8,sizeof(double));
	ohcstateout= (double*)mxCalloc(8,sizeof(double));
	Ihcin =(double*)mxCalloc(16,sizeof(double));
	Ihcstateout=(double*)mxCalloc(16,sizeof(double));
	c2in =(double*)mxCalloc(28,sizeof(double));
	c2out =(double*)mxCalloc(28,sizeof(double));
	c1in =(double*)mxCalloc(28,sizeof(double));
	c1out =(double*)mxCalloc(28,sizeof(double));
	tmpgain = (double*)mxCalloc(totalstim,sizeof(double));
    
	/** Calculate the center frequency for the control-path wideband filter
	    from the location on basilar membrane, based on Greenwood (JASA 1990) */



	if (species>1) /* for human */
    {
        /* Human frequency shift corresponding to 1.2 mm */
        bmplace = (35/2.1) * log10(1.0 + cf / 165.4); /* Calculate the location on basilar membrane from CF */
        centerfreq = 165.4*(pow(10,(bmplace+1.2)/(35/2.1))-1.0); /* shift the center freq */
    }
    
	/*==================================================================*/
	/*====== Parameters for the gain ===========*/
    
	
    if(species>1) gain = 52.0/2.0*(tanh(2.2*log10(cf/0.6e3)+0.15)+1.0); /* for human */
    /*gain = 52/2*(tanh(2.2*log10(cf/1e3)+0.15)+1);*/
    if(gain>60.0) gain = 60.0;  
    if(gain<15.0) gain = 15.0;
    
	/*====== Parameters for the control-path wideband filter =======*/
	bmorder = 3;
	Get_tauwb(cf,species,bmorder,Taumax,Taumin);
	taubm   = cohc*(Taumax[0]-Taumin[0])+Taumin[0];
	ratiowb = Taumin[0]/Taumax[0];
	/*====== Parameters for the signal-path C1 filter ======*/
	Get_taubm(cf,species,Taumax[0],bmTaumax,bmTaumin,ratiobm);
	bmTaubm  = cohc*(bmTaumax[0]-bmTaumin[0])+bmTaumin[0];
	fcohc    = bmTaumax[0]/bmTaubm;
    /*====== Parameters for the control-path wideband filter =======*/
	wborder  = 3;
    TauWBMax = Taumin[0]+0.2*(Taumax[0]-Taumin[0]);
	TauWBMin = TauWBMax/Taumax[0]*Taumin[0];
  
  	/*===============================================================*/
    /* Nonlinear asymmetry of OHC function and IHC C1 transduction function*/
	ohcasym  = 7.0;    
	ihcasym  = 3.0;
  	/*===============================================================*/
    /*===============================================================*/
    /* Prewarping and related constants for the middle ear */
     fp = 1e3;  /* prewarping frequency 1 kHz */
     C  = TWOPI*fp/tan(TWOPI/2*fp*tdres);

     if (species>1) /* for human */
     {
         /* Human middle-ear filter - based on Pascal et al. (JASA 1998)  */
         m11=1/(pow(C,2)+5.9761e+003*C+2.5255e+007);m12=(-2*pow(C,2)+2*2.5255e+007);m13=(pow(C,2)-5.9761e+003*C+2.5255e+007);m14=(pow(C,2)+5.6665e+003*C);             m15=-2*pow(C,2);					m16=(pow(C,2)-5.6665e+003*C);
         m21=1/(pow(C,2)+6.4255e+003*C+1.3975e+008);m22=(-2*pow(C,2)+2*1.3975e+008);m23=(pow(C,2)-6.4255e+003*C+1.3975e+008);m24=(pow(C,2)+5.8934e+003*C+1.7926e+008); m25=(-2*pow(C,2)+2*1.7926e+008);	m26=(pow(C,2)-5.8934e+003*C+1.7926e+008);
         m31=1/(pow(C,2)+2.4891e+004*C+1.2700e+009);m32=(-2*pow(C,2)+2*1.2700e+009);m33=(pow(C,2)-2.4891e+004*C+1.2700e+009);m34=(3.1137e+003*C+6.9768e+008);     m35=2*6.9768e+008;				m36=(-3.1137e+003*C+6.9768e+008);
         megainmax=2;
     };
   
        
		
            mey1  = m11*(-m12*statein[0]  + px[0]         - statein[6]);
            if (species>1) mey1= m11*(-m12*statein[0]-m13*statein[1]+m14*px[0]+m15*statein[6]+m16*statein[7]);
            mey2 = m21*(-m22*statein[2] - m23*statein[3] + m24*mey1 + m25*statein[0] + m26*statein[1]);
            mey3  = m31*(-m32*statein[4]- m33*statein[5] + m34*mey2 + m35*statein[2] + m36*statein[3]);
            


	
statetmp[0]=mey1;
//statetmp[1]=mey1[1];
statetmp[2]=mey2;
//statetmp[3]=mey2[1];
statetmp[4]=mey3;
//statetmp[5]=mey3[1];
statetmp[6]=px[0];
//statetmp[7]=px[1];
//totalsim==2
meout=mey3/megainmax;


   
		
			phasein=statein[8];
			wbgain=statein[18]; 
			tauwb=statein[17];
			
			  for (int i=0; i<8 ; i++)
		  {
			  rin[i]=statein[9+i];
		  }
					
		
		wbout1= WbGammaTone(meout,tdres,centerfreq,0,tauwb,wbgain,wborder,wbphase,phasein,rout,rin);

 
		statetmp[8]=wbphase[0];
		
		 for (int i=0; i<8 ; i++)
		  {
			  statetmp[i+9]=rout[i];
		  }
		
   
        wbout  = pow((tauwb/TauWBMax),wborder)*wbout1*10e3*__max(1,cf/5e3);
		
        ohcnonlinout = Boltzman(wbout,ohcasym,12.0,5.0,5.0); /* pass the control signal through OHC Nonlinear Function */

		  for (int i=0; i<8 ; i++)
		  {
			  ohcin[i]=statein[22+i];
		  }
		
				  
		ohcout = OhcLowPass(ohcnonlinout,tdres,600,0,1.0,2,ohcstateout,ohcin);/* lowpass filtering after the OHC nonlinearity */
     		
		
		
		 for (int i=0; i<8 ; i++)
		  {
			  statetmp[i+22]=ohcstateout[i];
		  }
		

		
		tmptauc1 = NLafterohc(ohcout,bmTaumin[0],bmTaumax[0],ohcasym); /* nonlinear function after OHC low-pass filter */
		
		tauc1    = cohc*(tmptauc1-bmTaumin[0])+bmTaumin[0];  /* time -constant for the signal-path C1 filter */

		rsigma   = 1/tauc1-1/bmTaumax[0]; /* shift of the location of poles of the C1 filter from the initial positions */

		if (1/tauc1<0.0) mexErrMsgTxt("The poles are in the right-half plane; system is unstable.\n");
       
		tauwb = TauWBMax+(tauc1-bmTaumax[0])*(TauWBMax-TauWBMin)/(bmTaumax[0]-bmTaumin[0]);
		
        statetmp[17]=tauwb;

	    wb_gain = gain_groupdelay(tdres,centerfreq,cf,tauwb,grdelay);
		
		  
		grd = grdelay[0]; 
        statetmp[20]=grd; 
		statetmp[18]=wb_gain; 
	


	
        /*====== Signal-path C1 filter ======*/
    	c1gainin =statein[60];
		c1phasein=statein[61];
	
		for (int c1=0; c1<28;c1++)
		{
			c1in[c1]=statein[c1+62];	
		}     
		 
		 c1filterouttmp = C1ChirpFilt(meout, tdres, cf, 0, bmTaumax[0], rsigma, c1gainin, c1phasein, c1out, c1in); /* C1 filter output */

		 
statetmp[60]=c1gainin;
statetmp[61]=c1phasein;
		 for (int i=0; i<28 ; i++)
		  {
			  statetmp[i+62]=c1out[i];
		  }

		 
	 
        /*====== Parallel-path C2 filter ======*/
		c2gainin =statein[30];
		c2phasein=statein[31];
	
		for (int c2=0; c2<28;c2++)
		{
			c2in[c2]=statein[c2+32];	
		}
	
c2filterouttmp  = C2ChirpFilt(meout, tdres, cf, 0, bmTaumax[0], 1/ratiobm[0], c2gainin, c2phasein, c2out, c2in); /* parallel-filter output*/

	 	

statetmp[30]=c2gainin;
statetmp[31]=c2phasein;
for ( int i=0; i<18; i++) statetmp[32+i]=c2out[i];
for ( int i=0; i<10; i++) statetmp[50+i]=c2out[i+18];



	    /*=== Run the inner hair cell (IHC) section: NL function and then lowpass filtering ===*/

        c1vihctmp  = NLogarithm(cihc*c1filterouttmp,0.1,ihcasym,cf);
	    
		 
		c2vihctmp = -NLogarithm(c2filterouttmp*fabs(c2filterouttmp)*cf/10*cf/2e3,0.2,1.0,cf); /* C2 transduction output */
    
		   
     	
		for(int i=0; i<16; i++)
		{
			Ihcin[i]=statein[i+90];
		}
		

			
			
		
                ihcouttmp[0] = IhcLowPass(c1vihctmp+c2vihctmp,tdres,3000,0,1,7,Ihcstateout,Ihcin);
				           
				for (int i=0; i<16 ; i++)
				{
					statetmp[i+90]=Ihcstateout[i];
				}	
				
				
				/////////////////////////
				////////////////////////////////////
				///////////////////////////////////////////////////
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//
				//
				 
   
   
    
 
    double synstrength,synslope,CI,CL,PG,CG,VL,PL,VI,exponOut,expoafterfilter,powerlawafterfilter;
	double cf_factor,PImax,kslope,Ass,Asp,TauR,TauST,Ar_Ast,PTS,Aon,AR,AST,Prest,gamma1,gamma2,k1,k2;
	double VI0,VI1,alpha,beta,theta1,theta2,theta3,vsat,tmpst,tmp,PPI,CIlast,temp;
     int implnt=1;
     int spont=100; 	
//if (lsr==1) spont=0.1;
 spont=100; 

	 
	 
		if (spont==100) cf_factor = __min(800,pow(10,0.29*cf/1e3 + 0.7));
       if (spont==4)   cf_factor = __min(50,2.5e-4*cf*4+0.2);
       if (spont==0.1) cf_factor = __min(1.0,2.5e-4*cf*0.1+0.15);              
	         
	   PImax  = 0.6;                /* PI2 : Maximum of the PI(PI at steady state) */
       kslope = (1+50.0)/(5+50.0)*cf_factor*20.0*PImax;            
       /* Ass    = 300*TWOPI/2*(1+cf/100e3); */  /* Older value: Steady State Firing Rate eq.10 */
       Ass    = 800*(1+cf/100e3);    /* Steady State Firing Rate eq.10 */

       if (implnt==1) Asp = spont*3.0;   /* Spontaneous Firing Rate if actual implementation */
       if (implnt==0) Asp = spont*2.75; /* Spontaneous Firing Rate if approximate implementation */
       TauR   = 2e-3;               /* Rapid Time Constant eq.10 */
       TauST  = 60e-3;              /* Short Time Constant eq.10 */
       Ar_Ast = 6;                  /* Ratio of Ar/Ast */
       PTS    = 3;                  /* Peak to Steady State Ratio, characteristic of PSTH */
   
       /* now get the other parameters */
       Aon    = PTS*Ass;                          /* Onset rate = Ass+Ar+Ast eq.10 */
       AR     = (Aon-Ass)*Ar_Ast/(1+Ar_Ast);      /* Rapid component magnitude: eq.10 */
       AST    = Aon-Ass-AR;                       /* Short time component: eq.10 */
       Prest  = PImax/Aon*Asp;                    /* eq.A15 */
       CG  = (Asp*(Aon-Asp))/(Aon*Prest*(1-Asp/Ass));    /* eq.A16 */
       gamma1 = CG/Asp;                           /* eq.A19 */
       gamma2 = CG/Ass;                           /* eq.A20 */
       k1     = -1/TauR;                          /* eq.8 & eq.10 */
       k2     = -1/TauST;                         /* eq.8 & eq.10 */
               /* eq.A21 & eq.A22 */
       VI0    = (1-PImax/Prest)/(gamma1*(AR*(k1-k2)/CG/PImax+k2/Prest/gamma1-k2/PImax/gamma2));
       VI1    = (1-PImax/Prest)/(gamma1*(AST*(k2-k1)/CG/PImax+k1/Prest/gamma1-k1/PImax/gamma2));
       VI  = (VI0+VI1)/2;
       alpha  = gamma2/k1/k2;       /* eq.A23,eq.A24 or eq.7 */
       beta   = -(k1+k2)*alpha;     /* eq.A23 or eq.7 */
       theta1 = alpha*PImax/VI; 
       theta2 = VI/PImax;
       theta3 = gamma2-1/PImax;
  
       PL  = ((beta-theta2*theta3)/theta1-1)*PImax;  /* eq.4' */
       PG  = 1/(theta3-1/PL);                        /* eq.5' */
       VL  = theta1*PL*PG;                           /* eq.3' */
       CI  = Asp/Prest;                              /* CI at rest, from eq.A3,eq.A12 */
       CL  = CI*(Prest+PL)/PL;                       /* CL at rest, from eq.1 */
   	
       if(kslope>=0)  vsat = kslope+Prest;                
       tmpst  = log(2)*vsat/Prest;
       if(tmpst<400) synstrength = log(exp(tmpst)-1);
       else synstrength = tmpst;
       synslope = Prest/log(2)*synstrength;
       
       
            tmp = synstrength*(ihcouttmp[0]);   
            if(tmp<400) tmp = log(1+exp(tmp));
            PPI = synslope/synstrength*tmp;           
         
		 
		
				CI=statein[106];
				CL=statein[107];
			
            CIlast = CI; 
            CI = CI + (tdres/VI)*(-PPI*CI + PL*(CL-CI));
            CL = CL + (tdres/VL)*(-PL*(CL - CIlast) + PG*(CG - CL));
            if(CI<0)
            {
                temp = 1/PG+1/PL+1/PPI;
                CI = CG/(PPI*temp);
                CL = CI*(PPI+PL)/PL;
            };
			statetmp[106]=CI;
			statetmp[107]=CL;
            exponOut = CI*PPI; //afagh

				
/* 
		for(int i=0; i<16; i++)
		{
			Ihcin[i]=statein[i+132];
		}
			
            expoafterfilter = IhcLowPass(exponOut,tdres,10000,0,1.0,7,Ihcstateout,Ihcin);     
			
		for (int i=0; i<16 ; i++)
		{
			statetmp[i+132]=Ihcstateout[i];
		}
*/

//statetmp[148]=exponOut;
				   
				////////
				////////////////////////
				////////////////////////////////////////////
				////////////////////////////////////////////////////////////////
				double alpha1, I1, alpha2, I2, sout1, sout2,n1,n2,n3,m1,m2,m3,m4,m5;
				 alpha1 = 2.5e-6*100e3;
				 alpha2 = 1e-2*100e3;
				 
			     I1=statein[108];
			     I2=statein[109];
				 
				 sout1  = __max( 0, exponOut+ noisein - alpha1*I1); //*/   /* No fGn condition */
				 sout2  = __max( 0, exponOut - alpha2*I2);

			     statetmp[110]=sout1;
		         statetmp[111]=sout2;
				
				
              

         




				
			/*	
				if (statein[130]==1)
				{
                    n1 = 1.0e-3*sout2;
                    n2 = n1; n3= n2; I2 = n3;	
	                m1 = 0.2*sout1;
                    m2 = m1;	m3 = m2;			
                    m4 = m3;	m5 = m4;
					I1 = m5;
					
	
	
				}
				
				if (statein[131]==1)
				{
                    n1 = 1.992127932802320*statein[111]+ 1.0e-3*(sout2 - 0.994466986569624*statein[127]);
                    n2 = 1.999195329360981*statein[112]+ n1 - 1.997855276593802*statein[111];
                    n3 = -0.798261718183851*statein[113]+ n2 + 0.798261718184977*statein[112];
					
					I2=n3;
					
					m1 = 0.491115852967412*statein[119] + 0.2*(sout1 - 0.173492003319319*statein[125]);
                    m2 = 1.084520302502860*statein[120] + m1 - 0.803462163297112*statein[119];
                    m3 = 1.588427084535629*statein[121] + m2 - 1.416084732997016*statein[120];
                    m4 = 1.886287488516458*statein[122] + m3 - 1.830362725074550*statein[121];
                    m5 = 1.989549282714008*statein[123] + m4 - 1.983165053215032*statein[122];
					I1 = m5; 
	
				}				
			
				else if (statein[130]==0 && statein[131]==0)
				{
				
				
				 n1 = 1.992127932802320*statein[111] - 0.992140616993846*statein[108]+ 1.0e-3*(sout2 - 0.994466986569624*statein[127] + 0.000000000002347*statein[126]);
                 n2 = 1.999195329360981*statein[112] - 0.999195402928777*statein[109]+ n1  - 1.997855276593802*statein[111] + 0.997855827934345*statein[108];
                 n3 =-0.798261718183851*statein[113] - 0.199131619873480*statein[110]+ n2  + 0.798261718184977*statein[112] + 0.199131619874064*statein[109];
				 
				 
				 
				 
				 
				
				 I2 = n3;
				 
				 m1 = 0.491115852967412*statein[119] - 0.055050209956838*statein[114] + 0.2*(sout1- 0.173492003319319*statein[125]+ 0.000000172983796*statein[124]);
                 m2 = 1.084520302502860*statein[120] - 0.288760329320566*statein[115] + m1 - 0.803462163297112*statein[119] + 0.154962026341513*statein[114];
                 m3 = 1.588427084535629*statein[121] - 0.628138993662508*statein[116] + m2 - 1.416084732997016*statein[120] + 0.496615555008723*statein[115];
                 m4 = 1.886287488516458*statein[122] - 0.888972875389923*statein[117] + m3 - 1.830362725074550*statein[121] + 0.836399964176882*statein[116];
                 m5 = 1.989549282714008*statein[123] - 0.989558985673023*statein[118] + m4 - 1.983165053215032*statein[122] + 0.983193027347456*statein[117];				
				
				 I1 = m5; 
				}
			
			     statetmp[130]=0;
				 statetmp[131]=0;
				 statetmp[111]=n1;
				 statetmp[112]=n2;
				 statetmp[113]=n3;
				 
				 statetmp[119]=m1;
				 statetmp[120]=m2;
				 statetmp[121]=m3;
				 statetmp[122]=m4;
				 statetmp[123]=m5;
				 
				 statetmp[128]=I1;
				 statetmp[129]=I2;
				 
				//////////////////////////////////////////
				///////////////////////////////////////////////////
				///////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////////
				
				*/
				
	/*for(int i=0; i<16; i++)
		{
			Ihcin[i]=statein[i+150];
		}
  powerlawafterfilter = IhcLowPass(sout1,tdres,10000,0,1.0,7,Ihcstateout,Ihcin);     
			
		for (int i=0; i<16 ; i++)
		{
			statetmp[i+150]=Ihcstateout[i];
		}
*/

				ihcout[0]	 =sout1;//+sout2;
				ihcout[0]    = ihcout[0]/(1+0.75e-3*ihcout[0]);
			//	statetmp[149]=sout1;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			

 for(i=0;i<112;i++)
	{
state[i]=statetmp[i];
	}
  /* Freeing dynamic memory allocated earlier */

    mxFree(ihcouttmp); mxFree(statetmp); 
  
    mxFree(tmpgain); mxFree(rin); mxFree(ohcin);mxFree(rout);mxFree(ohcstateout);
	mxFree(Ihcin); mxFree(c2in); mxFree(c2out); mxFree(c1in);mxFree(c1out);mxFree(Ihcstateout);
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

double C1ChirpFilt(double x, double tdres,double cf, int n, double taumax, double rsigma,double c1gainin,double c1phasein, double *c1out, double *c1in)
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
	C1gain_norm=c1gainin;
	C1initphase=c1phasein;
	
	C1input[1][1]=c1in[0];
	C1input[1][2]=c1in[1];
	C1input[1][3]=c1in[2];
	C1input[2][1]=c1in[3];
	C1input[2][2]=c1in[4];
	C1input[2][3]=c1in[5];
	C1input[3][1]=c1in[6];
	C1input[3][2]=c1in[7];
	C1input[3][3]=c1in[8];
	C1input[4][1]=c1in[9];
	C1input[4][2]=c1in[10];
	C1input[4][3]=c1in[11];
	C1input[5][1]=c1in[12];
	C1input[5][2]=c1in[13];
	C1input[5][3]=c1in[14];
	C1input[6][1]=c1in[15];
	C1input[6][2]=c1in[16];
	C1input[6][3]=c1in[17];

	C1output[1][1]=c1in[18];
	C1output[2][1]=c1in[19];
	C1output[3][1]=c1in[20];
	C1output[4][1]=c1in[21];
	C1output[5][1]=c1in[22];
	C1output[1][2]=c1in[23];
	C1output[2][2]=c1in[24];
	C1output[3][2]=c1in[25];
	C1output[4][2]=c1in[26];
	C1output[5][2]=c1in[27];

	}

	C1gain_norm=c1gainin;
	C1initphase=c1phasein;
     
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
	      	  
	  c1out[0]=C1input[1][1];
	  c1out[1]=C1input[1][2];
	  c1out[2]=C1input[1][3];
	  c1out[3]=C1input[2][1];
	  c1out[4]=C1input[2][2];
	  c1out[5]=C1input[2][3];
	  c1out[6]=C1input[3][1];
	  c1out[7]=C1input[3][2];
	  c1out[8]=C1input[3][3];
	  c1out[9]=C1input[4][1];
	  c1out[10]=C1input[4][2];
	  c1out[11]=C1input[4][3];
	  c1out[12]=C1input[5][1];
	  c1out[13]=C1input[5][2];
	  c1out[14]=C1input[5][3];
	  c1out[15]=C1input[6][1];
	  c1out[16]=C1input[6][2];
	  c1out[17]=C1input[6][3];
	  
      c1out[18]=C1output[1][1];
	  c1out[19]=C1output[2][1];
	  c1out[20]=C1output[3][1];
	  c1out[21]=C1output[4][1];
	  c1out[22]=C1output[5][1];
	  c1out[23]=C1output[1][2];
	  c1out[24]=C1output[2][2];
	  c1out[25]=C1output[3][2];
	  c1out[26]=C1output[4][2];
	  c1out[27]=C1output[5][2];
	  
     return (c1filterout);
}  

/* -------------------------------------------------------------------------------------------- */
/** Parallelpath C2 filter: same as the signal-path C1 filter with the OHC completely impaired */

double C2ChirpFilt(double xx, double tdres,double cf, int n, double taumax, double fcohc,double c2gainin,double c2phasein,double *c2out, double *c2in)
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
	C2gain_norm=c2gainin;
	C2initphase=c2phasein;
	C2input[1][1]=c2in[0];
	C2input[1][2]=c2in[1];
	C2input[1][3]=c2in[2];
	C2input[2][1]=c2in[3];
	C2input[2][2]=c2in[4];
	C2input[2][3]=c2in[5];
	C2input[3][1]=c2in[6];
	C2input[3][2]=c2in[7];
	C2input[3][3]=c2in[8];
	C2input[4][1]=c2in[9];
	C2input[4][2]=c2in[10];
	C2input[4][3]=c2in[11];
	C2input[5][1]=c2in[12];
	C2input[5][2]=c2in[13];
	C2input[5][3]=c2in[14];
	C2input[6][1]=c2in[15];
	C2input[6][2]=c2in[16];
	C2input[6][3]=c2in[17];

	C2output[1][1]=c2in[18];
	C2output[2][1]=c2in[19];
	C2output[3][1]=c2in[20];
	C2output[4][1]=c2in[21];
	C2output[5][1]=c2in[22];
	C2output[1][2]=c2in[23];
	C2output[2][2]=c2in[24];
	C2output[3][2]=c2in[25];
	C2output[4][2]=c2in[26];
	C2output[5][2]=c2in[27];

	}

	C2gain_norm=c2gainin;
	C2initphase=c2phasein;



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
	  
	 
	  c2out[0]=C2input[1][1];
	  c2out[1]=C2input[1][2];
	  c2out[2]=C2input[1][3];
	  c2out[3]=C2input[2][1];
	  c2out[4]=C2input[2][2];
	  c2out[5]=C2input[2][3];
	  c2out[6]=C2input[3][1];
	  c2out[7]=C2input[3][2];
	  c2out[8]=C2input[3][3];
	  c2out[9]=C2input[4][1];
	  c2out[10]=C2input[4][2];
	  c2out[11]=C2input[4][3];
	  c2out[12]=C2input[5][1];
	  c2out[13]=C2input[5][2];
	  c2out[14]=C2input[5][3];
	  c2out[15]=C2input[6][1];
	  c2out[16]=C2input[6][2];
	  c2out[17]=C2input[6][3];
	  
      c2out[18]=C2output[1][1];
	  c2out[19]=C2output[2][1];
	  c2out[20]=C2output[3][1];
	  c2out[21]=C2output[4][1];
	  c2out[22]=C2output[5][1];
	  c2out[23]=C2output[1][2];
	  c2out[24]=C2output[2][2];
	  c2out[25]=C2output[3][2];
	  c2out[26]=C2output[4][2];
	  c2out[27]=C2output[5][2];
	  
	  return (c2filterout); 
}   

/* -------------------------------------------------------------------------------------------- */
/** Pass the signal through the Control path Third Order Nonlinear Gammatone Filter */

double WbGammaTone(double x,double tdres,double centerfreq, int n, double tau,double gain,int order, double *wbphase, double phasein, double*rout,double* rin )
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
  
  rout[0]=REAL(wbgtfl[0]);
  rout[1]=REAL(wbgtfl[1]);
  rout[2]=REAL(wbgtfl[2]);
  rout[3]=REAL(wbgtfl[3]);
  rout[4]=IMAG(wbgtfl[0]);
  rout[5]=IMAG(wbgtfl[1]);
  rout[6]=IMAG(wbgtfl[2]);
  rout[7]=IMAG(wbgtfl[3]);

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

double OhcLowPass(double x,double tdres,double Fc, int n,double gain,int order, double* ohcstateout, double* ohcin)
{ 
  double ohc[4],ohcl[4];

  double c,c1LP,c2LP;
  int i,j;

  
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
  
  ohcstateout[0]= ohc[0];
  ohcstateout[1]= ohc[1];
  ohcstateout[2]= ohc[2];
  ohcstateout[3]= ohc[3];
  
  ohcstateout[4]=ohcl[0];
  ohcstateout[5]=ohcl[1];
  ohcstateout[6]=ohcl[2];
  ohcstateout[7]=ohcl[3];
  
  return(ohc[order]);
}
/* -------------------------------------------------------------------------------------------- */
/* Get the output of the IHC Low Pass Filter  */

double IhcLowPass(double x,double tdres,double Fc, int n,double gain,int order, double* Ihcstateout, double* Ihcin)
             
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
  
  Ihcstateout[0]= Ihc[0];
  Ihcstateout[1]= Ihc[1];
  Ihcstateout[2]= Ihc[2];
  Ihcstateout[3]= Ihc[3];
  Ihcstateout[4]= Ihc[4];
  Ihcstateout[5]= Ihc[5];
  Ihcstateout[6]= Ihc[6];
  Ihcstateout[7]= Ihc[7];
  
  Ihcstateout[8]=Ihcl[0];
  Ihcstateout[9]=Ihcl[1];
  Ihcstateout[10]=Ihcl[2];
  Ihcstateout[11]=Ihcl[3];
  Ihcstateout[12]=Ihcl[4];
  Ihcstateout[13]=Ihcl[5];
  Ihcstateout[14]=Ihcl[6];
  Ihcstateout[15]=Ihcl[7];
  
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