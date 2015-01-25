#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#define c 3e10
#define N 500
#define mp 1.67e-24
#define G0 200
#define M0 0.005
#define epsilon 1
#define epsilon_b 1.e-8


//#define DT 0.00005
//#define DG 0.03

#define DT 0.0005
#define DG 0.01

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

int explicitscheme(double *gamma, double *fgammat, double *fgammatp1){
	int i;
	double *G;
	G= (double *) calloc(N,sizeof(double));

	for(i=0;i<N;i++){
		G[i] = fgammat[i];//*gammadot(gamma[i]);
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
}

double Source(double x){
	if((x<3000.) || (x>8.e4)){
		return 0.;
	}
	else{
		return pow(x,-2.2);
	}
}

int main(void){
	
	int i,j,k;
	double deltaGamma;
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

	double gdotp,gdotm,a;
	double step;
	double integralp, integralm;
	step=exp(1./N*log(1e5));
	double V2,V3;
	for(i=0;i<N;i++){
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
		
	}

	FILE *output;
	FILE *nbscat;
	output=fopen("test2.txt","a");
	nbscat=fopen("partnum1.txt","a");
			for(i=0;i<N;i++){
				fprintf(output,"%g\t%g\n",gamma[i],fgamma[i]);
			}

	double Tot=0;
	for(i=0;i<2500;i++){

		/////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////
		// Initializing round:
		printf("i=%d\n",i);
		deltaGamma = G[N-1]-G[N-2];
 	     	fgammatp1[N-1]=fgamma[N-1]/(1. + (DT*0.005*gamma[N-1]*gamma[N-1] )/ deltaGamma);
	     	//fgammatp1[0]=fgamma[N-1]/(1. + (DT*0.005*gamma[N-1]*gamma[N-1] )/ deltaGamma);
	      	for(j = N-2; j>=1; j--){

			// Gamma is not computed in the middle of the cell. Try with middle.
//			deltaGamma = 0.5*(gamma[j+1]-gamma[j-1]); //Half steps are at j+.5 and j-.5
			deltaGamma = 0.5*(G[j]-G[j-1]); //Half steps are at j+.5 and j-.5
			//deltaGamma = pow(10.,0.5*(log(gamma[j+1])-log(gamma[j])))-pow(10.,0.5*(log(gamma[j])-log(gamma[j-1])));

			//Set the coeffs. At some point this should be virtualized
			gdotp=0.005*gamma[j+1]*gamma[j+1];		// Initialisation: no SSC
			gdotm=0.005*gamma[j]*gamma[j];
			V3 = (DT*gdotp )/ deltaGamma;
			V2 = 1. + (DT*gdotm )/ deltaGamma;
		  	fgammatp1[j] = (fgamma[j]+Source(gamma[j])*DT + V3*fgammatp1[j+1])/V2;
	//  spectra[i+1][j] = (source[j] - V3[j]*spectra[i+1][j+1])/V2[j];	 	
		}
		deltaGamma = G[1]-G[0];
		fgammatp1[0]=fgamma[0]+(DT*0.005*gamma[1]*gamma[1]* fgammatp1[1])/ deltaGamma;
      		//printf("%g\t%g\n",deltaGamma,fgammatp1[0],fgamma[0]);
		//getchar();

/*

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
			printf("%d\n",counter);
			for(j = N-2; j>=1; j--){

				deltaGamma = 0.5*(G[j]-G[j-1]); //Half steps are at j+.5 and j-.5

				// COmpute the integral
				integralp=0.;
				if(epsilon_b<1./G[j+1]){
					for(k=j+1;k<N-2;k++){
						if(1./(epsilon_b*G[j+1])>G[k]){
							integralp += G[k]*G[k]*fgammatp1[k]*(G[k+1]-G[k]);
						}
					}
				}
				integralm=0.;
				if(epsilon_b<1./G[j]){
					for(k=j;k<N-2;k++){
						if(1./(epsilon_b*G[j])>G[k]){
							integralm += G[k]*G[k]*fgammatp1[k]*(G[k+1]-G[k]);
						}
					}
				}	
				//integral=0.;			
				//printf("integral=%g\n",integral);
				gdotp=0.005*gamma[j+1]*gamma[j+1]*(1.+0.01*integralp);		// Initialisation: no SSC
				gdotm=0.005*gamma[j]*gamma[j]*(1.+0.01*integralm);		
				V3 = (DT*gdotp )/ deltaGamma;
				V2 = 1. + (DT*gdotm )/ deltaGamma;
			  	fgammatemp[j] = (fgamma[j] + V3*fgammatemp[j+1])/V2;

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
				if(i%50==0){
					fprintf(output,"%.13g\t%.13g\n",gamma[j],fgammatp1[j]);
				}
			}
				Tot=0;
				for(k=0;k<N-1;k++){
					Tot+=fgamma[k]*(G[k+1]-G[k]);
				}
				fprintf(nbscat,"%.13g\n",Tot);	
		}


*/

		for(j=0;j<N;j++){
			fgamma[j]=fgammatp1[j];
			if(fgamma[j]<1e-30){
				fgamma[j]=0.;
			}
		}







	//	if((j%50==0) && (j!=0)){
		if((i%50==0)){
			for(k=0;k<N;k++){
				fprintf(output,"%.13g\t%.13g\n",gamma[k],fgamma[k]);
			}
		}
		Tot=0;
		for(j=0;j<N-1;j++){
			Tot+=fgamma[j]*(G[j+1]-G[j]);
		}
		fprintf(nbscat,"%g\n",Tot);
	}
	fclose(nbscat);
	fclose(output);

	free(G);
	free(fgammatemp);
	free(fgammatp1);
	free(fgamma);
	free(gamma);

	return 0;
}


