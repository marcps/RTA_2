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

void D_matrix(int currPos, ,double *Dv,double **D,double *Tempr)
{
	//posf is the current position, posi is the position to where we wanna
	//integrate.
	int i;
	double sum=0;
	double dT=(TAU_F-TAU_0)/NITER;
	//we cannot sum from the current position because there's no defined
	//Temperature!!!
	for (j=currPos;j>=0;j--)
	{
		D[currPos][j]=Dv[j]/Dv[currPos];
	}
}

void D_evolve(int currPos, double *Dv, double *Tempr)
{
	int i,j,sum;
	dT=(TAU_F-TAU_0)/NITER;
	Dv=Dv*exp(-dT/tau_eq(Tempr[currPos];
}

int main (int argc, char* argv[])
{
	int i,j;
	double **D,*Tempr, *Energy,*Dv;
	double currTau,currTauPrime,dT=(TAU_F-TAU_0)/NITER;
	double S_t,F_t;
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
	Energy=(double *)malloc((NITER+1)*sizeof(double));
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	Energy[0]=ENERGY_0;
	Tempr[0]=Tempr(Energy[0]);
	/*%%%%%%%%%%%%%%%%%%%%%%% TIME INTEGRATION %%%%%%%%%%%%%%%%%%%%*/
	for(i=0;i<NITER;i++)
	{
		currTau=TAU_0+(i+0.5)dT;
		D_evolve(i, Dv, Tempr);
		S_t=D[i][i]*Ht(ex(currTau));
		F_t=0;
		for(j=0;j<=i;j++)
		{
			currTauPrime=currTau+(j+0.5)*dT;
			F_t+=D[i][j] * Energy[i] * Ht((currTau/currTauPrime)*(currTau/currTauPrime)-1) * dT/tau_eq(Tempr[i]);
		}
		Energy[i]=S_t+F_t;
		Tempr[i]=Tempr(Energy[i]);
	}
	return 0;
}

		
