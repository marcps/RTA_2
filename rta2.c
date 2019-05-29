#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
double const PI=3.14159265359;//definition of PI
double const TAU_0=0.1;
double const TAU_F=100;
double const EX_0=1;
double const ENERGY_0=1;
int const NITER=10000;

double ex(double tau)
{
	return (1+EX_0)*(tau*tau/TAU_0/TAU_0);
}

double Temp(double E)
{
        /*Returns the TEMPERATURE*/
        return sqrt(PI*sqrt(E/3));
}

double tau_eq (double tempr)
{
	double gamma=1.0;
	return gamma/tempr;
}

double Ht(double exx)
{
	double a=(PI/180)*atan(sqrt(exx))/sqrt(exx);
	double b=1/(1+exx);
	return 0.5*(b+a);
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
		Dv[0]=(double)1;
	}
	int i,j,sum;
	double dT;
	dT=(TAU_F-TAU_0)/NITER;
	Dv[currPos]=Dv[currPos-1]*exp(-dT/tau_eq(Tempr[currPos]));//This is indeed Tempr[currPos] because it is defined that way!
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
	int i,j;
	double **D,*Tempr,*deriT, *Energy,*deriE,*Dv;
	double currTau,currTauPrime,dT=(TAU_F-TAU_0)/NITER;
	double S_t,F_t,sn;
	 //-----------------------FILE INIT-----------------------------
        if(argc<5)
        {
                printf("\n[!]ERROR. Please type the names of the output");
                printf(".dat files!\n\n[info]run %s ",argv[0]);
                printf("<1>.dat <2>.dat <3>.dat <4>.dat\n\n\n");
                return 1;
        }
        FILE *fp,*fp1,*fp2,*fp3;

        if(strcmp(argv[1],argv[2])==0 || strcmp(argv[2],argv[3])==0 || strcmp(argv[1],argv[3])==0)
        {
                printf("\n[!]ERROR. Filenames must be DIFFERENT!\n\n");
                return 2;
        }
        fp=fopen(argv[1],"w");
        fp1=fopen(argv[2],"w");
        fp2=fopen(argv[3],"w");
        fp3=fopen(argv[4],"w");
        //-------------------------------------------------------------

	/*%%%%%%%%%%%%%%%%%%%%%%% D matrix initialization%%%%%%%%%%%%%%%*/
	D = (double **)malloc((NITER+1)*sizeof(double));
        if(D==NULL)
        {
                printf("NOT Enough memory!\n Exitting...\n");
                return 4;
        }
        for(i=0;i<NITER;i++)
        {
                D[i]=(double *)malloc((NITER+1)*sizeof(double));
                if(D[i]==NULL)
                {
                        printf("NOT Enough memory!\n Exitting...\n");
                        return 5;
                }
        }
	Dv=(double *)malloc((NITER+1)*sizeof(double));
	Tempr=(double *)malloc((NITER+1)*sizeof(double));
	deriT=(double *)malloc((NITER+1)*sizeof(double));
	Energy=(double *)malloc((NITER+1)*sizeof(double));
	deriE=(double *)malloc((NITER+1)*sizeof(double));
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	Energy[0]=ENERGY_0;
	Tempr[0]=Temp(Energy[0]);
	/*%%%%%%%%%%%%%%%%%%%%%%% TIME INTEGRATION %%%%%%%%%%%%%%%%%%%%*/
	for(i=0;i<NITER;i++)
	{
		currTau=TAU_0+(i+0.5)*dT;
		D_evolve(i, Dv, Tempr);
		D_matrix(i,Dv,D,Tempr);
		S_t=D[i][i]*Ht(ex(currTau));
		F_t=0;
		for(j=0;j<=i;j++)
		{
			currTauPrime=currTau-(j+0.5)*dT;
			F_t+=D[i][j]*Energy[j] * Ht((currTau/currTauPrime)*(currTau/currTauPrime)-1) * dT/tau_eq(Tempr[i]);
		}
		printf("Energy=%.15f  Temperature=%.15f\n",Energy[i],Tempr[i]);

		Energy[i+1]=S_t+F_t;
		Tempr[i+1]=Temp(Energy[i]);
		deriE[i]=(Energy[i+1]-Energy[i])/dT;
		deriT[i]=(Tempr[i+1]-Tempr[i])/dT;
		//sn=fabs(currTau*deriE[i]/(Energy[i]+ /\/\/\/\/\/\/\/\/\/\/\/\/\ WE LACK A PRESSURE EXPRESSION!!!! /\/\/\/\/\/\/\/\/\/\/\
		//%%%%%%%%%%%%%%%% FILE PRINTING %%%%%%%%%%%%%%%%%%%%%%%%%%%
		print_to_file(fp,currTau*Tempr[i],Tempr[i]*pow(currTau,1/3));
		print_to_file(fp1,currTau*Tempr[i],1); //HERE GOES SN!!! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
		print_to_file(fp2,currTau*Temprr[i],currTau*deriT[i]/Tempr[i]);
		print_to_file(fp3,currTau*Tempr[i],deriE[i]/Energy[i]);
		//add time functionality as well.
	}
	return 0;
}

