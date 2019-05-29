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

double DD(int posf,int posi, double *Tempr)
{
	//posf is the current position, posi is the position to where we wanna
	//integrate.
	int i;
	double sum=0;
	double dT=(TAU_F-TAU_0)/NITER;
	//we cannot sum from the current position because there's no defined
	//Temperature!!!
	for (i=posi ; i<posf ; i++)
	{
		sum+=1.0/tau_eq(Tempr[i]);
	}
	return exp(-dT*sum);
}

void D_matrix(int currPos, double **D, double *Tempr)
{
	int i,j,sum;
	//This calculates all the necessary columns for the row currPos
	//starting from D[currPos][currPos] up to D[currPos][1]????
	for(i=0;i<currPos;i++)
	{
		D[currPos][currPos-i]=DD(currPos,i,Tempr[i]);//NO ESTA BE
	}
}
int main (int argc, char* argv[])
{
	int i,j;
	double **D,*Tempr, *Energy;
	double currTau,currTauPrime,dT=(TAU_F-TAU_0)/NITER;
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
	/*%%%%%%%%%%%%%%%%%%%%%%% TIME INTEGRATION %%%%%%%%%%%%%%%%%%%%*/
	for(i=1;i<NITER;i++)
	{
		currTau=TAU_0+(i+0.5)dT;
		D_matrix(i-1,D,Tempr);
		eny_1=D[i-1][i-1]*Ht(ex(currTau));
		eny_2=0;
		for(j=0;j<i-1;j++)
		{
			//This is clever but kind of convoluted
			currTauPrime=currTau-(j+0.5)*dT;
			eny_2+=D[i-1][j]*Energy[i-1-j]/tau_eq(Tempr[i-1-j])*Ht((currTau*currTau/currTauPrime/currTauPrime)-1);
		}
		Energy[i]=eny_1+dT*eny_2;
		Tempr[i]=Tempr(Energy[i]);
	}
	return 0;
}

		
