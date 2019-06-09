#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
double const PI=3.14159265359;//definition of PI
double const TAU_0=0.25;
double const TAU_F=20.0;
double const XI_0=0.000000000000000;  // Momentum-space anisotrpy
double const T_0=1.5; // Mev%c
double const ETA=(double)1/((double)4*3.14159265359);
double const TSTEP=0.005;

double ex(double tau)
{
	return ((1+XI_0)*(tau*tau/(TAU_0*TAU_0)))-1.0;
}

double Temp(double E)
{
        /*Returns the TEMPERATURE*/
        return sqrt(sqrt(E*3.0)/PI);
}

double tau_eq (double tempr)
{
	double gamma_eq=5.0/ETA;
	return gamma_eq/tempr;
}

double Ht(double x)
{
	if(x==0.000000000000000){
		//This if ensures the 0/0 case!
		return (double)1;
	}
	double a=atan(sqrt(x))/sqrt(x);
	double b=1.0/(1.0+x);
	return 0.5*(b+a);
}

double Ht_0()
{
	if(XI_0==0.000000000000000){
		return (double)1;
	}
	else{
		return 0.5*(atan(sqrt(XI_0))/sqrt(XI_0)+(1/(1+XI_0)));
	}
}

void D_matrix(int currPos,double *Dv,double **D,double *Tempr)
{

	int j;
	D[currPos][currPos]=(double)1;
	for (j=currPos-1;j>=0;j--)
	{
		D[currPos][j]=Dv[currPos]/Dv[j];
	}
}

void D_evolve(int currPos, double *Dv, double *Tempr)
{
	if(currPos==0){
		//In the first iteration, exp(0)=1
		Dv[0]=(double)1;
	}
	else{
		double dT=TSTEP;
		Dv[currPos]=Dv[currPos-1]*exp(-dT/tau_eq(Tempr[currPos]));//This is indeed Tempr[currPos] because it is defined that way!
	}
}

void debug(double a,char* name){
	printf("\n\n       DEBUGG:\n          %s=%f\n\n",name,a);
	char c=getchar();
}

void print_to_file(FILE* f,double x, double y)
{
	fprintf(f,"%.15f %.15f\n",x,y);
}

int main (int argc, char* argv[])
{
	double Energy_0=(PI*PI*T_0*T_0*T_0*T_0)/(double)3;
	int i,j;
	double **D,*Tempr,*deriT, *Energy,*deriE,*Dv, *Pressure_L, *Pressure_T;
	double currTau,currTauPrime,dT=TSTEP;
	int niter=(int)(TAU_F-TAU_0)/TSTEP;
	double S_t,F_t,sn,Py,etaEff;
	 //-----------------------FILE INIT-----------------------------
        if(argc<5)
        {
                printf("\n[!]ERROR. Please type the names of the output");
                printf(".dat files!\n\n[info]run %s ",argv[0]);
                printf("<1>.dat <2>.dat <3>.dat <4>.dat\n\n\n");
                return 1;
        }
        FILE *fp,*fp1,*fp2,*fp3,*fp4;

        if(strcmp(argv[1],argv[2])==0 || strcmp(argv[2],argv[3])==0 || strcmp(argv[1],argv[3])==0 || strcmp(argv[1],argv[4])==0 || strcmp(argv[2],argv[4])==0 || strcmp(argv[3],argv[4])==0 || strcmp(argv[1],argv[5])==0||strcmp(argv[2],argv[5])==0||strcmp(argv[3],argv[5])==0||strcmp(argv[4],argv[5])==0)
        {
                printf("\n[!]ERROR. Filenames must be DIFFERENT!\n\n");
                return 2;
        }
        fp=fopen(argv[1],"w");
        fp1=fopen(argv[2],"w");
        fp2=fopen(argv[3],"w");
        fp3=fopen(argv[4],"w");
	fp4=fopen(argv[5],"w");
        //-------------------------------------------------------------

	/*%%%%%%%%%%%%%%%%%%%%%%% D matrix initialization%%%%%%%%%%%%%%%*/
	D = (double **)malloc((niter+1)*sizeof(double));
        if(D==NULL)
        {
                printf("NOT Enough memory!\n Exitting...\n");
                return 4;
        }
        for(i=0;i<niter;i++)
        {
                D[i]=(double *)malloc((niter+1)*sizeof(double));
                if(D[i]==NULL)
                {
                        printf("NOT Enough memory!\n Exitting...\n");
                        return 5;
                }
        }
	Dv=(double *)malloc((niter+1)*sizeof(double));
	Tempr=(double *)malloc((niter+1)*sizeof(double));
	deriT=(double *)malloc((niter+1)*sizeof(double));
	Energy=(double *)malloc((niter+1)*sizeof(double));
	deriE=(double *)malloc((niter+1)*sizeof(double));
	Pressure_L=(double *)malloc((niter+1)*sizeof(double));
	Pressure_T=(double *)malloc((niter+1)*sizeof(double));
	if(Dv==NULL||Tempr==NULL||deriT==NULL||Energy==NULL||deriE==NULL||Pressure_L==NULL||Pressure_T==NULL){
		printf("NOT Enough memory!\n Exitting...\n");
		return 6;
	}
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	Energy[0]=Energy_0;
	Tempr[0]=Temp(Energy[0]);
	/*%%%%%%%%%%%%%%%%%%%%%%% TIME INTEGRATION %%%%%%%%%%%%%%%%%%%%*/
	for(i=0;i<niter;i++)
	{
		currTau=TAU_0+(i+0.5)*dT;
		D_evolve(i, Dv, Tempr);
		D_matrix(i,Dv,D,Tempr);
		F_t=D[i][0]*Ht(ex(currTau))/Ht_0();
		S_t=0;
		for(j=0;j<=i;j++)
		{
			currTauPrime=TAU_0+(j+0.5)*dT;
			S_t+=D[i][j]*Energy[j]*Ht((currTau/currTauPrime)*(currTau/currTauPrime)-1)*dT/tau_eq(Tempr[i]);
		}

		Energy[i+1]=S_t+F_t;
		Tempr[i+1]=Temp(Energy[i+1]);
		deriE[i]=(Energy[i+1]-Energy[i])/dT;
		deriT[i]=(Tempr[i+1]-Tempr[i])/dT;

		Py=currTau*deriE[i]+4*Energy[i]/3;
		Pressure_L[i]=Energy[i+1]/(double)3 - Py;
		Pressure_T[i]=Energy[i+1]/(double)3 + Py/(double)2;

		//This has to be a very small number
		sn=fabs(currTau*deriE[i]/(Energy[i]+Pressure_L[i]))-(double)1;
		etaEff=9.0*currTau*currTau/4.0*(deriT[i]+Tempr[i+1]/(3.0*currTau));

		printf("i=%d ;Energy=%.15f  ;deriE=%.15f  ;Temperature=%.15f   ;deriT=%.15f  ;SN = %.15f\n",i,Energy[i],deriE[i],Tempr[i],deriT[i], sn);
		//%%%%%%%%%%%%%%%% FILE PRINTING %%%%%%%%%%%%%%%%%%%%%%%%%%%
		print_to_file(fp,currTau,etaEff);
		print_to_file(fp1,currTau,sn);
		print_to_file(fp2,currTau,Pressure_L[i]/Pressure_T[i]);
		print_to_file(fp3,currTau,Energy[i+1]*pow(currTau/TAU_0,3/(double)2)*TAU_0*TAU_0*TAU_0*TAU_0*TAU_0);
		//print_to_file(fp4,currTau,Pressure_L[i]/Pressure_T[i]);
		//add time functionality as well.
	}
	fclose(fp);
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	free(Tempr);
	free(Dv);
	free(D);
	free(deriT);
	free(Energy);
	free(deriE);
	free(Pressure_T);
	free(Pressure_L);
	return 0;
}
