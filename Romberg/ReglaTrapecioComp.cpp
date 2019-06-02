#include <iostream>
#include <cmath>


using namespace std;

double simpson38(double (*f)(double),double a,double b,int n)
{
    if(n%3==0)
    {
        double h=(b-a)/n;
        double xi=0;
        double s1=0,s2=0,s3=0,rs;
        int limite2=n-3;

        if(limite2!=0)
        {
            int i=1,j=2,k=3;
            while(k<=n-3)
            {
                s3=s3+f(a+h*k);
                k=k+3;
                while(i<=n-2)
                {
                    s1=s1+f(a+h*i);
                    i=i+3;
                    while(j<=n-1)
                    {
                        s2=s2+f(a+h*j);
                        j=j+3;
                    }
                }
            }
        }
        else
        {
            s3=0;
            int x=1,y=2;
            while(y<=n-1)
            {
                s2=s2+f(a+h*y);
                y=y+3;
                while(x<=n-2)
                {
                    s1=s1+f(a+h*x);
                    x=x+3;
                }
            }
        }
        rs=((3*h)/8)*(f(a)+3*s1+3*s2+2*s3+f(b));
        return rs;
    }
}


double trapecioComp(double (*fun)(double), double a, double b, int n)
{
  double h=(b-a)/(n*1.0);
  double sol=0;
  for(int i=1;i<n;i++)
  {
    sol+=fun((i*h)+a);
  }

  return h/2.0*( fun(a)+2.0* sol + fun(b) ) ;
}



double simpson(double (*fun)(double), double a, double b, int n)
{

  double h=(b-a)/n;


  double temp1=0.0;
  for(int j=1;j<=(n/2-1);j++)
    temp1+=fun(a+((2*j)*h));

  double temp2=0.0;
  for(int j=1;j<=n/2;j++)
    temp2+=fun(a+((2*j-1.0)*h));

  double res=(h/3.0)*(fun(a) + 2.0 * temp1 + 4.0*temp2 + fun(b) );

  return res;
}


double simpson3oct(double (*fun)(double), double a, double b, int n)
{
  double res=fun(a)+fun(b);
  double h=(b-a)/n;

  for(int i=1;i<n-2;i=i+3)
  {
    res+=3.0*fun(a+i*h);
  }

  for(int i=2;i<n-1;i=i+3)
  {
    res+=3.0*fun(a+i*h);
  }

  for(int i=3;i<n-3;i=i+3)
  {
    res+=2.0*fun(a+i*h);
  }

  return (3.0/8.0)*h *res;

}


double pow2(double n)
{
  return sin(pow(n,2));
}





int main()
{
  cout<<simpson(pow2,0,1,12)<<endl;
  cout<<trapecioComp(pow2,0,1,12)<<endl;
  cout<<simpson38(pow2,0,1,12)<<endl;
  return 0;
}
