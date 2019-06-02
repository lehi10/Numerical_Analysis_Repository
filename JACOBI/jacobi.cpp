#include <iostream>
#include <stdlib.h>
#include <math.h>
using namespace std;


double **crear(int n, int m)
{
    double **A=new double*[n];
    for(int i=0;i<n;i++)
        A[i]=new double[m];
    return A;
}

double *crear(int n)
{
    double *vec=new double[n];
    return vec;
}

void llenarManual(double **&A, int n, int m)
{
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            cin>>A[i][j];
}

void llenarManual(double *&A, int n)
{
    for(int i=0;i<n;i++)
            cin>>A[i];
}

void llenar(double **&A, int n, int m)
{
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            A[i][j]=rand()%10+1;
}

void llenar(double *&b, int n)
{
    for(int i=0;i<n;i++)
        b[i]=rand()%30+1;
}

void copiar(double *&x, double *&xant,int n)
{
    for(int i=0;i<n;i++)
        xant[i]=x[i];
}

void imprimir(double **&A, int n, int m)
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
            cout<<A[i][j]<<" ";
        cout<<endl;
    }
}


void imprimir(double *&b,int n)
{
    for(int i=0;i<n;i++)    
        cout<<b[i]<<" ";
    cout<<endl;
}

double norma(double *&x, double *&xant,int m)
{
    double sum=0;
    for(int i=0;i<m;i++)
        sum+=fabs(x[i]-xant[i]);
    return sum;
}


double *jacobi(double **&A,double *&b, int n,int m,double *&x, double tol)
{
    double *xant=crear(m);
    int cont=0;
    do
    {  
        copiar(x,xant,n);
        for(int i=0;i<n;i++)
        {
            double denom=1;
            double sum=b[i];
            for(int j=0;j<m;j++)
            {
                if(i==j)
                    denom*=A[i][j];
                else
                    sum-=A[i][j]*x[j];
            }
            x[i]=sum/denom;
        }
    }
    while(norma(x,xant,m)>tol);
    
    return x;
}

int main()
{
    
    int     n=4,m=4;
    double  tol=0.000001;

    double **A=crear(n,m);
    double *b=crear(n);
    double *x0=crear(m);

    cout<<"Ingresar matriz A de tamaño "<<n<<" x "<<m<<" : "<<endl;
    llenarManual(A,n,m);
    cout<<"Ingresar vector de terminos independientes de tamaño "<<n<<" :"<<endl;
    llenarManual(b,n);
    cout<<"Vector inicial"<<endl;
    llenar(x0,m);    
    imprimir(x0,m);
    
    double *res=jacobi(A,b,n,m,x0,tol);
    cout<<"---------------------------------------------------"<<endl;
    cout<<"Solución aproximada aplicando el metodo de Jacobi :"<<endl;
    imprimir(res,m);

    return 0;
}
