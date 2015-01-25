#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#define c 3e10
#define mp 1.67e-24
#define G0 200
#define M0 0.005
#define epsilon 1

double function(double a, double b){
	return -sin(a);

}


int dynamical_evol(double *DQa, double *DQb, double r){

	DQa[0]=-(DQb[0]*DQb[0]-1.)/(DQb[2]) * r*r;
	DQa[1]=(DQb[0]-1)*r*r;
	DQa[2]=(1.-epsilon)*(DQb[0]-1)*r*r+r*r;
	
	return 0;
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

	int i,ii;
	int j;
	double t,r;
	double step=0.01;
	double tmax=10000;
	t=1;
	r=0;
	double rho,M,m,gam,dr,dt;
	rho=1.;
	double hydron[3];
	double hydronp1[3];
	double summ[3],sum1[3];

	hydron[0]= 100;
	hydron[1]= 0;
	hydron[2]= 10000;

	FILE *output;
	output=fopen("test1.txt","a");

	while(t<tmax){
			
		for(i=0;i<4;i++){
			for(ii=0;ii<3;ii++){
				summ[ii]=0.;
				for(j=0;j<i;j++){
					summ[ii]+=a[i][j]*k[j-1][ii];
				//printf("i=%d\tj=%d\tsum=%g\tk=%g\ta=%g\n",i,j,summ,k[j-1],a[i][j]);
				//getchar();
				}
				summ[ii]= hydron[ii]+step*summ[ii];
			}
			//dynamical_evol(double *DQa, double *DQb, double r){
			dynamical_evol( k[i],  summ,  t+d[i]*step);
			//	k[ii][i] = function(t+d[i]*step,hydron[ii]+step*summ);

		}
		
		for(ii=0;ii<3;ii++){
			sum1[ii]=0.;
			for(i=0;i<4;i++){
				sum1[ii]+=step*b[i]*k[i][ii];
			}
			hydronp1[ii] = hydron[ii]+sum1[ii];
		}
		t=t+step;
		fprintf(output,"%g\t%g\t%g\t%g\t%g\n",t,hydronp1[0],hydronp1[1],hydronp1[2]);
		//printf("%g\n",t);
		for(ii=0;ii<3;ii++){
			hydron[ii]=hydronp1[ii];
		}

	}


	fclose(output);
	return 0;
}


