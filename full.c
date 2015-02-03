#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/// Constante definition:
#define c 2.99792458e10
#define mp 1.67262178e-24
#define me 9.10938291e-28
#define sigmaT 6.6524574e-25
#define pi 3.14159265359
#define Bcrit 4.414e13
#define prefac_cool 3.247963784e-8

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/// Initial parameters of the shell
#define G0 300      // Initial Lorentz factor
#define E0 1e52     // Total energy of the blast wave: ergs
//#define M0 0.005

#define xi_B 1
#define xe_e 0.5
#define n 1
#define s -3
#define rzero 1e14



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/// Build-in parameters:
#define N 1000

//#define DT 0.0005
//#define DG 0.01




double abs_val(double x){
	if(x>0){
		return x;
	}
	else{
		return -x;
	}
}



double gammadot(double x){
	return -0.02*x*x;

}

double S(double gamma){
	return 0.;//pow(gamma,-2.2);
}

/*int explicitscheme(double *gamma, double *fgammat, double *fgammatp1){
	int i;
	double *G;
	G= (double *) calloc(N,sizeof(double));

	for(i=0;i<N;i++){
		G[i] = fgammat[i];// *gammadot(gamma[i]);
	}

	for(i=1;i<N-1;i++){
		fgammatp1[i] = DT/(2.*DG)*(G[i-1] - G[i+1])+fgammat[i];//+S(gamma[i]);
		//if(fgammatp1[i]<0){
		//	fgammatp1[i]=0.;
		//}
		fgammatp1[0]=0.;
		fgammatp1[N]=0.;
	}
	free(G);
	return 0;
}*/


//Source(gamma[j],gmin,gmax)
double Source(double x, double gmin, double gmax){
	if((x<gmin) || (x>gmax)){
		return 0.;
	}
	else{
		return pow(x,s);
	}
}

//int main(void){
int evol_distrib_function(double *G, double *gamma, double *fgamma, double *fgammatp1, double ub, double normS, double *coolingG, double gmin, double gmax, double DT, int itnum, double dV){


	int j,k;
	double deltaGamma;
    double epsilonb=sqrt(8.*pi*ub)/Bcrit;
//	double *gamma;
//	gamma= (double *) calloc(N,sizeof(double));
//	double *fgamma;
//	fgamma= (double *) calloc(N,sizeof(double));

	double *fgammatemp;
	fgammatemp= (double *) calloc(N,sizeof(double));
//	double *fgammatp1;
//	fgammatp1= (double *) calloc(N,sizeof(double));

//	double *G;
//	G= (double *) calloc(N+1,sizeof(double));

	double gdotp,gdotm,a;
	//double step;
	double integralp, integralm;
	//step=exp(1./N*log(1e5));
	double V2,V3;
/*	for(i=0;i<N;i++){
		//gamma[i]=1.+DG*i;
		gamma[i]=pow(step,i);
		if(i<N-1){
			G[i]=0.5*(gamma[i]+gamma[i]*step);
		}
		else{
			G[N-1]=0.5*(gamma[i]+gamma[i]*step);
		}
		fgamma[i]=pow(gamma[i],-2.2);
		if((gamma[i]<3000.) || (gamma[i]>8.e4)){
			fgamma[i]=0.;
		}

	} */

	FILE *output;
	//FILE *nbscat;
	output=fopen("distrib.txt","a");
	//nbscat=fopen("partnum1.txt","a");
	/*		for(i=0;i<N;i++){
				fprintf(output,"%g\t%g\n",gamma[i],fgamma[i]);
			}
*/
	double Tot=0;
	//for(i=0;i<2500;i++){

		/////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////
		// Initializing round:
		//printf("i=%d\n",i);
		deltaGamma = G[N-1]-G[N-2];
        fgammatp1[N-1]=fgamma[N-1]/(1. + (DT*gamma[N-1]*gamma[N-1] *prefac_cool*ub)/ deltaGamma);
	     	//fgammatp1[0]=fgamma[N-1]/(1. + (DT*0.005*gamma[N-1]*gamma[N-1] )/ deltaGamma);
        for(j = N-2; j>=1; j--){

			// Gamma is not computed in the middle of the cell. Try with middle.
//			deltaGamma = 0.5*(gamma[j+1]-gamma[j-1]); //Half steps are at j+.5 and j-.5
			deltaGamma = 0.5*(G[j]-G[j-1]); //Half steps are at j+.5 and j-.5
			//deltaGamma = pow(10.,0.5*(log(gamma[j+1])-log(gamma[j])))-pow(10.,0.5*(log(gamma[j])-log(gamma[j-1])));

			//Set the coeffs. At some point this should be virtualized
			gdotp=prefac_cool*ub*gamma[j+1]*gamma[j+1];		// Initialisation: no SSC
			gdotm=prefac_cool*ub*gamma[j]*gamma[j];
			V3 = (DT*gdotp )/ deltaGamma;
			V2 = 1. + (DT*gdotm )/ deltaGamma;
		  	fgammatp1[j] = (fgamma[j]+Source(gamma[j],gmin,gmax)*DT + V3*fgammatp1[j+1])/V2;
	//  spectra[i+1][j] = (source[j] - V3[j]*spectra[i+1][j+1])/V2[j];
		}
		deltaGamma = G[1]-G[0];
		fgammatp1[0]=fgamma[0]+(DT*prefac_cool*ub*gamma[1]*gamma[1]* fgammatp1[1])/ deltaGamma;
      		//printf("%g\t%g\n",deltaGamma,fgammatp1[0],fgamma[0]);
		//getchar();



		/////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////
		// Loop round:
		// Criteria !!
		double norm=1.;
		int counter=0.;
		while(norm>1e-9){
			norm=0.;
			fgammatemp[N-1]=fgammatp1[N-1];
			counter+=1;
			//printf("counter%d\n",counter);
			for(j = N-2; j>=1; j--){

				deltaGamma = 0.5*(G[j]-G[j-1]); //Half steps are at j+.5 and j-.5

				// COmpute the integral
				integralp=0.;
				if(epsilonb<1./G[j+1]){
					for(k=j+1;k<N-2;k++){
						if(1./(epsilonb*G[j+1])>G[k]){
							integralp += G[k]*G[k]*fgammatp1[k]*(G[k+1]-G[k]);
						}
					}
				}
				integralm=0.;
				if(epsilonb<1./G[j]){
					for(k=j;k<N-2;k++){
						if(1./(epsilonb*G[j])>G[k]){
							integralm += G[k]*G[k]*fgammatp1[k]*(G[k+1]-G[k]);
						}
					}
				}
				//integral=0.;
				//printf("integral=%g\n",integral);
				gdotp=prefac_cool*gamma[j+1]*gamma[j+1]*ub*(1.+me*c*c*c*sigmaT*2.*integralp/(3.*dV));		// Initialisation: no SSC
				gdotm=prefac_cool*gamma[j]*gamma[j]*ub*(1.+me*c*c*c*sigmaT*2.*integralm/(3.*dV));
                /*printf("Cooling: %g\tdV%g\n",me*c*c*c*sigmaT*2.*integralm/(3.*dV),dV);
                if(j==500){
                    getchar();
                }*/
				coolingG[j]=gdotm;
				V3 = (DT*gdotp )/ deltaGamma;
				V2 = 1. + (DT*gdotm )/ deltaGamma;
			  	fgammatemp[j] = (fgamma[j] +normS*Source(gamma[j],gmin,gmax)*DT+ V3*fgammatemp[j+1])/V2;

				// Norm computation:
				if(fgammatp1[j]!=0){
					a=abs_val((fgammatemp[j]-fgammatp1[j])/fgammatp1[j]);
					//printf("loop a=%.13g\tfgammatemp[j]=%g\tfgammatp1[j]=%g\n",a,fgammatemp[j],fgammatp1[j]);
					//getchar();
				}
				else{
					a=fgammatemp[j];
				}
				if(norm<a){
					norm=a;
				}
			}
			//printf("norm=%g\n",norm);
			//getchar();

			for(j = N-2; j>=1; j--){		// As to be done after because the values are needed to compute the integral

				fgammatp1[j]=fgammatemp[j];
				if(fgammatp1[j]<1e-20){
					fgammatp1[j]=0.;
				}
				//if(counter%10==0){
				//	fprintf(output,"%.10g\t%.10g\n",gamma[j],fgammatp1[j]);
				//	getchar();
				//}
			}
				Tot=0;
				for(k=0;k<N-1;k++){
					Tot+=fgamma[k]*(G[k+1]-G[k]);
				}
			//	fprintf(nbscat,"%.13g\n",Tot);
			//	if(counter%10==0){
					//fprintf(output,"%.10g\t%.10g\n",gamma[j],fgammatp1[j]);
					//getchar();
			//	}
				//printf("norm=%g\n",norm);
		}




		for(j=0;j<N;j++){
			//fgamma[j]=fgammatp1[j];
			if(fgammatp1[j]<1e-30){
				fgammatp1[j]=0.;
			}
			/*if((gamma[j]<1000) && (fgamma[j]>1)){
                printf("Problem in cooling\n");
                getchar();
			}*/
		}
			for(j = N-2; j>=1; j--){
				if(itnum%1000==0){
					fprintf(output,"%g\t%g\n",gamma[j],fgammatp1[j]);
				}
            }

        fclose(output);



	//	if((j%50==0) && (j!=0)){
	/*	if((i%50==0)){
			for(k=0;k<N;k++){
				fprintf(output,"%.13g\t%.13g\n",gamma[k],fgamma[k]);
			}
		}
		Tot=0;
		for(j=0;j<N-1;j++){
			Tot+=fgamma[j]*(G[j+1]-G[j]);
		}
		fprintf(nbscat,"%g\n",Tot);
	//}
	//fclose(nbscat);
	//fclose(output);
*/
	//free(G);
	free(fgammatemp);
	//free(fgammatp1);
	//free(fgamma);
	//free(gamma);
    if(itnum%1000==0){
	//	getchar();
	}
    //printf("End of function\n");
    //getchar();

	return 0;
}







int dynamical_evol(double *DQa, double *DQb, double r){

    // Lorentz factor
	DQa[0]=-(DQb[0]*DQb[0]-1.)*r*r/(DQb[2]);
	//printf("dgamma=%g\n",DQa[0]);

	// Internal energy
	DQa[1]=(DQb[0]-1.)*r*r;
	//printf("dE=%g\n",DQa[1]);

	// Mass
	DQa[2]=(DQb[0]-1.)*r*r+r*r;    // Purely radiative, no factor epsilon.
	//printf("dM=%g\n\n",DQa[2]);

	return 0;
}

//				dynamical_evol_cool( k[i],  summ,  t+d[i]*step , hydron,epsilon);
//dynamical_evol_cool( k[i],  summ,  t+d[i]*step , hydronp1);
int dynamical_evol_cool(double *DQa, double *DQb, double r, double epsilon){


	DQa[0]=-(DQb[0]*DQb[0]-1.)/(DQb[2]) * r*r;
	DQa[1]=(DQb[0]-1)*r*r;
	DQa[2]=(1.-epsilon)*(DQb[0]-1)*r*r+r*r;

	return 0;
}


double radiative_efficiency(double *fgamma, double *fgammatp1, double *G, double *coolingG){
	int i;
	double epsilon=0.;
	for(i=1;i<N;i++){
		epsilon+=me*c*c*0.5*(fgamma[i]+fgammatp1[i])*coolingG[i]*(G[i+1]-G[i]);
	}



	return epsilon;
}

int main(void){



	// Dynamics
	double a[4][4],b[4],d[4],k[4][3];
	a[0][0]=0;
	a[1][0]=0.5;
	a[2][0]=0;
	a[3][0]=0;
	a[0][1]=0;
	a[1][1]=0;
	a[2][1]=0.5;
	a[3][1]=0;
	a[0][2]=0;
	a[1][2]=0;
	a[2][2]=0;
	a[3][2]=1;
	a[0][3]=0;
	a[1][3]=0;
	a[2][3]=0;
	a[3][3]=0;


	b[0]=1./6.;
	b[1]=1./3.;
	b[2]=1./3.;
	b[3]=1./6.;

	d[0]=0;
	d[1]=1./2.;
	d[2]=1./2.;
	d[3]=1.;

	int i,ii,itnum;
	int j;
	double t,r,ta;
    r=1.;
	double step=5.;
	double rmax=200;
	t=1;
	ta=0.;
    double dtcom;
    double dV;
	//double rho,M,m,gam,dr,dt;
	 // rho=1.;
	double hydron[3];
	double hydrontemp[3];
	double hydronp1[3];
    double fullsol[4];
	double summ[3],sum1[3];

	hydron[0]= G0;
	hydron[1]= 0;
	hydron[2]= E0/(4.*pi*mp*n*rzero*rzero*rzero*c*c*G0);

	hydrontemp[0]= 10.*G0;
	hydrontemp[1]= 1;
	hydrontemp[2]= 1e5;

	hydronp1[0]= 10.*G0;
	hydronp1[1]= 5;
	hydrontemp[2]= 1e5;

	FILE *output;
    FILE *output1;
    //FILE *output2;
	output=fopen("dynamics.txt","a");
    output1=fopen("dynamics_unit.txt","a");
    //output2=fopen("dynamics_adia.txt","a");

    printf("G=%g\tM=%g\n",hydron[0],hydron[2]);
    getchar();

	////// Distribution function at r should be also known:
	//double deltaGamma;
	double *gamma;
	gamma= (double *) calloc(N,sizeof(double));
	double *fgamma;
	fgamma= (double *) calloc(N,sizeof(double));

	double *fgammatemp;
	fgammatemp= (double *) calloc(N,sizeof(double));
	double *fgammatp1;
	fgammatp1= (double *) calloc(N,sizeof(double));

	double *G;
	G= (double *) calloc(N+1,sizeof(double));

	double *coolingG;
	coolingG= (double *) calloc(N+1,sizeof(double));

	//double gdotp,gdotm,a;
	double stepG;
	stepG=exp(1./N*log(1e9));			// Increase energy max !
	for(i=0;i<N;i++){
		//gamma[i]=1.+DG*i;
		gamma[i]=pow(stepG,i);
		if(i<N-1){
			G[i]=0.5*(gamma[i]+gamma[i]*stepG);
		}
		else{
			G[N-1]=0.5*(gamma[i]+gamma[i]*stepG);
		}
		fgamma[i]=0.;//pow(gamma[i],-2.2);
		fgammatemp[i]=0.;
		//if((gamma[i]<3000.) || (gamma[i]>8.e4)){
		//	fgamma[i]=0.;
		//}

	}

	double norm=0.;
	double dE_intdt=0.;
	double Nzero=0.;
	double bet=0.;
	double B2=0.;
	double ub=0.;
	double gmin=0.;
	double gmax=0.;
	double epsilon=0.;

    itnum=0;

	while(r<rmax){
        itnum+=1;
        printf("t=%g\tr=%g\tepsilon=%g\n",t,r,epsilon);
        //getchar();
		// Initialization to find the radiative efficiency: adiabatic approximation
		for(i=0;i<4;i++){
			for(ii=0;ii<3;ii++){
				summ[ii]=0.;
                //printf("here, summ=%g\n",summ[ii]);
				for(j=0;j<i;j++){
					summ[ii]+=a[i][j]*k[j][ii];
					//printf("here, summ=%g\n",summ[ii]);
                    //printf("i=%d\tj=%d\tsum=%g\tk=%g\ta=%g\n",i,j,summ[ii],k[j][ii],a[i][j]);
				//getchar();
				}
				summ[ii]= hydron[ii]+step*summ[ii];
                //printf("summ=%g\thydro=%g\tstep=%g\n",summ[ii],hydron[ii], step);
			}
			//dynamical_evol(double *DQa, double *DQb, double r){
			dynamical_evol( k[i],  summ,  r+d[i]*step);
            //printf("ki0=%g\tki1=%g\tki2=%g\n",k[i][0],k[i][1],k[i][2]);
			//	k[ii][i] = function(t+d[i]*step,hydron[ii]+step*summ);

		}

		for(ii=0;ii<3;ii++){
			sum1[ii]=0.;
			for(i=0;i<4;i++){
				sum1[ii]+=step*b[i]*k[i][ii];
			}
			hydronp1[ii] = hydron[ii]+sum1[ii];
		}

        fullsol[0]=hydronp1[0];
        fullsol[1]=hydronp1[1]*4.*pi*mp*n*rzero*rzero*rzero*c*c;
        fullsol[2]=hydronp1[2]*4.*pi*mp*n*rzero*rzero*rzero;
        fullsol[3]=r*rzero;

        //printf("Gamma=%.13g\tEint=%.13g\tM=%.13g\n",hydronp1[0],hydronp1[1],hydronp1[2]);
        //printf("r=%g\tGamma=%.13g\tEint=%.13g\tM=%.13g\n",fullsol[3],fullsol[0],fullsol[1],fullsol[2]);
        //getchar();
		/// Loop for convergence with cooling
		norm=1.;


		while(norm>1e-9){

			// Compute new distribution function:
			// U_b is needed, as well as the normalisation of the source function:
			B2=32.*pi*xi_B*n*mp*fullsol[0]*fullsol[0]*c*c;
			ub = B2/(8.*pi);		// B2 is B^2
			gmin = mp*fullsol[0]/me;
			gmax=4e7/sqrt(sqrt(B2));

            //printf("B=%g\tgmin=%g\tgmax=%g\n",sqrt(B2),gmin,gmax);
			// Normalisation of the rate function:
			bet=sqrt(hydronp1[0]*hydronp1[0]-1.)/hydronp1[0];

			dE_intdt=hydronp1[0]*(hydronp1[0]-1.)*4.*pi*r*r*mp*bet*c*c*c*n*rzero*rzero;

			Nzero=xe_e*(-s-2.)*pow(gmin,-(s+2.))*dE_intdt/(me*c*c);
			dtcom=step*rzero/(bet*c*hydronp1[0]);
			dV=4.*pi*r*r*r*rzero*rzero*rzero/hydronp1[0];
            //printf("bet=%g\tdE_intdt=%g\tNzero=%g\tdt=%g\n",bet,dE_intdt,Nzero,dtcom);
            //getchar();

			// int evol_distrib_function(double *G, double *gamma, double *fgamma, double *fgammatp1, double ub, double normS, double DT){
			evol_distrib_function(G,gamma, fgamma, fgammatp1, ub, Nzero, coolingG , gmin,gmax,dtcom,itnum,dV);       // fgamma distribution function at t, fgammatp1 at d+dt with the parameters computed from RK4, ub magnetic energy density, normS normalisation of source function     To be computed before !!

			// Compute epsilon:
			epsilon=radiative_efficiency(fgammatemp,fgammatp1,G,coolingG)/dE_intdt;
            //printf("epsilon=%g\n",epsilon);
            //getchar();

			for(i=0;i<4;i++){
				for(ii=0;ii<3;ii++){
					summ[ii]=0.;
					for(j=0;j<i;j++){
						summ[ii]+=a[i][j]*k[j][ii];
					//printf("i=%d\tj=%d\tsum=%g\tk=%g\ta=%g\n",i,j,summ,k[j-1],a[i][j]);
					//getchar();
					}
					summ[ii]= hydron[ii]+step*summ[ii];
				}
				//dynamical_evol(double *DQa, double *DQb, double r){
				dynamical_evol_cool( k[i],  summ,  r+d[i]*step ,epsilon);    /// Use previous hydrop1 to compute cooling
				//	k[ii][i] = function(t+d[i]*step,hydron[ii]+step*summ);

			}

			for(ii=0;ii<3;ii++){
				sum1[ii]=0.;
				for(i=0;i<4;i++){
					sum1[ii]+=step*b[i]*k[i][ii];
				}
				hydrontemp[ii] = hydron[ii]+sum1[ii];
			}


			norm=abs_val((hydrontemp[0]-hydronp1[0])/hydronp1[0])+abs_val((hydrontemp[1]-hydronp1[1])/hydronp1[1]);
			//printf("norm=%g\t%g\t%g\t%g\n",norm,hydrontemp[0],hydronp1[0],hydronp1[2]);
			//fprintf(output,"%g\t%g\t%g\t%g\n",t,hydronp1[0],hydronp1[1],hydronp1[2]);
			//printf("%g\n",t);

			for(ii=0;ii<3;ii++){
				hydronp1[ii]=hydrontemp[ii];
			}

		}

        for(i=0;i<N;i++){
            fgammatemp[i]=fgamma[i];
            fgamma[i]=fgammatp1[i];
        }

		for(ii=0;ii<3;ii++){
			hydron[ii]=hydronp1[ii];
		}
		r=r+step;
		ta+=step*rzero/(2.*hydron[0]*hydron[0]*c*bet);
		fprintf(output,"%.13g\t%.13g\t%.13g\t%.13g\t%.13g\n",r,ta,hydron[0],hydron[1],hydron[2]);
		fprintf(output1,"%.13g\t%.13g\t%.13g\t%.13g\t%.13g\n",fullsol[3],ta,fullsol[0],fullsol[1],fullsol[2]);
        //r=r+step;
        //printf("r=%g\tGamma=%.13g\tEint=%.13g\tM=%.13g\n",r,hydron[0],hydron[1],hydron[2]);
        //printf("End of loop\n\n");
	}

	fclose(output1);
	fclose(output);
	return 0;
}




