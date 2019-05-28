#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double const PI=3.14159265359;//definition of PI
double const TAU_0=;
double const TAU_F=;
double const EX_0;
int const NITER=;

double ex(double tau)
{
	return (1+EX_0)*(tau*tau/TAU_0/TAU_0);
}

double Tempr(double E)
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
	double a=atan(sqrt(exx))/sqrt(exx);
	double b=1/(1+exx);
	return 0.5*(b+a);
}

void D_write_column(double **D, int j, double val_d)
{
	//Writes the value val_d in the whole j column
	int i;
	for (i=j;i<NITER;i++)
	{
		D[i][j]=val_d;
	}
}
double D_calculate(int posf,int posi, double *Tempr)
{
	//posf is the current position, posi is the position to where we wanna
	//integrate.
	int i;
	double sum=0;
	//we cannot sum from the current position because there's no defined
	//Temperature!!!
	for (i=posf-1 ; i>=posi ; i--)
	{
		sum+=dT/tau_eq(Tempr[i]);
	}
	return exp(-sum);
}

int main (int argc, char* argv[])
{
	int i,j;
	double **D,*Tempr, *Energy;
	double currTau,dT=(TAU_F-TAU_0)/NITER;
	double eny_1,eny_2;
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
	Tempr=(double *)malloc((NITER+1)*sizeof(double));
	Energy=(double *)malloc((NITER+1)*sizeof(double));
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	Energy[0]=ENERGY_0;
	Tempr[0]=Tempr(Energy[0]);
	D_write_column(D,0,ENERGY_0);
	/*%%%%%%%%%%%%%%%%%%%%%%% TIME INTEGRATION %%%%%%%%%%%%%%%%%%%%*/
	for(i=1;i<NITER;i++)
	{
		currTau=TAU_0+(i+0.5)dT;

		eny_1=D_calculate(i,0,Tempr)*Ht(ex(currTau));
