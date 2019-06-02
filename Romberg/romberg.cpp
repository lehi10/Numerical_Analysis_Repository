#include <iostream> 
#include <cmath>
using namespace std;

double trapecioComp(double (*fun)(double), double a, double b, int n)
{
  double h=(b-a)/(n*1.0);
  double sol=0;
  for(int i=1;i<n;i++)
    sol+=fun((i*h)+a);
  
  return h/2.0*( fun(a)+2.0* sol + fun(b) ) ;
}


double R( double (*fun)(double) ,int i, int j, double a, double b)
{
	if(j==1)	
		return trapecioComp(fun,a,b,pow(2,i-1));	
	double res =( pow(4,j-1)*R(fun,i,j-1,a,b)-R(fun,i-1,j-1,a,b) ) / ( pow(4, j-1) - 1) ;
	return res;
}


void romberg(double (*fun)(double) ,double  a,double b, int n)
{
	for(int i=1;i<=n;i++)
	{
		for(int j=1;j<=i;j++)
		{	
			cout<< R(fun,i,j,a,b) <<"\t";
		}
		cout<<endl;
	}
}

double fun(double n)
{
  return sin(n);
}

int main()
{
	
romberg(fun,0,M_PI/2,7);

return 0;
}
