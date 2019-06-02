
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


void identidad(double **&A, int n, int m)
{
    for(int i=0;i<n;i++)
        A[i][i]=1;
}

void copiar(double **&A, double **&B,int n, int m)
{
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            B[i][j]=A[i][j];
}

double ** productoMat(double **&A, int n,int m,double**&B,int h,int w)
{
    double **result=crear(n,w);

    for(int f=0;f<n;f++)
        for(int c=0;c<w;c++)
            for(int i=0;i<m;i++)
                result[f][c]+=A[f][i]*B[i][c];
    
    return result;
        
}

double * productoMatVec(double **&A, int n,int m,double*&b,int w)
{
    double *result=crear(n);    
    for(int f=0;f<n;f++)
        for(int i=0;i<m;i++)
            result[f]+=A[f][i]*b[i];
    
    return result;
        
}

double **sumaMat(double **&A, double**&B ,int n,int m)
{
    double **result=crear(n,m);
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            result[i][j]=A[i][j]+B[i][j];
    return result;
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

void pivoteo(double **&A,int n,int m,int irowMax, int colStart,int rowtoChange)
{
	if(irowMax==rowtoChange)
		return;
	double temp;

	for(int i=colStart; i<m ; i++)
	{
		temp=A[irowMax][i];
		A[irowMax][i]=A[rowtoChange][i];
		A[rowtoChange][i]=temp;
	}
}

void maxIndex(double **&A,int n, int m,int irow ,int icol  ,int &indexMax)
{
	indexMax=irow;
	for(int i=irow;i<n;i++)
		if(fabs(A[i][icol]) > fabs(A[indexMax][icol]))		
			indexMax=i;		
	
}

void escalonaPiv(double **&U,double **&P, double **&L,int n,int m)
{
	for(int i=0;i<m;i++)
	{
		int indexMax;
        maxIndex(U,n,m,i,i,indexMax);    

        pivoteo(U,n,m,indexMax,0,i);
        pivoteo(P,n,n,indexMax,0,i);
        pivoteo(L,n,n,indexMax,0,i);
    
		for(int j=i+1;j<n;j++)
		{
			double num=U[j][i];
			double denom=U[i][i];

            L[j][i]=num/denom;

			for(int k=i;k<m;k++)			
				U[j][k]=U[j][k]-((U[i][k]/denom)*num );
		}	

	}
}


void escalona(double **&U,double **&L,int n,int m)
{
	for(int i=0;i<m;i++)
	{   
		for(int j=i+1;j<n;j++)
		{
			double num=U[j][i];
			double denom=U[i][i];

            L[j][i]=num/denom;
            if(L[j][i] == 0)
                cout<<"[ERROR] No Admite descomposicón en LU. "<<endl;
			for(int k=i;k<m;k++)			
				U[j][k]=U[j][k]-((U[i][k]/denom)*num );
		}	
	}
}


void sustRegre(double **&A, int n, int m, double *&res)
{
	double suma=0;
	res[n-1]=A[n-1][m-1]/A[n-1][n-1];	
	for(int i=n-2;i>=0;i--)
	{
		suma=0;
		for(int j=i+1;j<n;j++)
			suma=suma+A[i][j]*res[j];
		res[i]=(A[i][m-1]-suma)/A[i][i];
	}
}

void sustProg(double **&A, int n, int m, double *&res)
{
	double suma=0;
	res[0]=A[0][m-1]/A[0][0];
    
	for(int i=1;i<n;i++)
	{
		suma=0;
		for(int j=0;j<i;j++)
			suma=suma+A[i][j]*res[j];
		res[i]=(A[i][m-1]-suma)/A[i][i];
	}    
    
}

double** contact(double **A,int n, int m, double *b)
{
    double **res=crear(n,m+1);
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            res[i][j]=A[i][j];
    for(int i=0;i<n;i++)
        res[i][m]=b[i];
    return res;
}

double * ResolverSistConLU(double **A, double *b, int n, int m)
{

    //HALLANDO LA DESPOMPOSICÓN EN LU
    //Creamos las matrices L y U
    double ** L=crear(n,n);
    double ** U=crear(n,m);
    copiar(A,U,n,m);
    escalona(U,L,n,m);
    //Sumar la matriz L con la identidad
    double **I=crear(n,n);
    identidad(I,n,n);
    L=sumaMat(L,I,n,n);
 
    
    //Sustitución Progresiva
    // Ly=Pb
    double **Ly=contact(L,n,n,b);
    double *y=crear(n);
    sustProg(Ly,n,m+1,y);

    //Sustitución Regresiva
    // Ux =y
    double **Ux=contact(U,n,m,y);
    double *x=crear(n);

    sustRegre(Ux,n,m+1,x);
    
    return x;
}

double * ResolverSistConPLU(double **A, double *b, int n, int m)
{
    double ** P=crear(n,n);
    double ** L=crear(n,n);
    double ** U=crear(n,m);
    
    copiar(A,U,n,m);
    identidad(P,n,m);

    //Escalonamiento con Pivoteo
    escalonaPiv(U,P,L,n,m);

    //Sumar la matriz L con una identidad
    double **I=crear(n,n);
    identidad(I,n,n);
    L=sumaMat(L,I,n,n);

    //Verificación
        //double **PA=productoMat(P,n,n,A,n,m);
        //double **LU=productoMat(L,n,n,U,n,m);
        //imprimir(PA,n,m);
        //imprimir(LU,n,m);


    //Sustitución Progresiva
    // Ly=Pb
    double *Pb=productoMatVec(P,n,n,b,n);
    double **Ly=contact(L,n,n,Pb);
    double *y=crear(n);
    sustProg(Ly,n,m+1,y);
    
    //Sustitución Regresiva
    //  Ux=y
    double **Ux=contact(U,n,m,y);
    double *x=crear(n);
    sustRegre(Ux,n,n+1,x);
    return x;
}



int main()
{
    int n=4, m=4;
    double ** A=crear(n,m);
    double   *b=crear(n);
    
    //Llena Manualmente
    llenarManual(A,n,m);
    llenarManual(b,n);

    //LLena con datos aleatorios
    //llenar(A,n,m);
    //llenar(b,n);
    

    double *res2=ResolverSistConPLU(A,b,n,m);
    double *res1=ResolverSistConLU(A,b,n,m); 
    
    imprimir(res1,n);
    imprimir(res2,n);

    return 0;
}