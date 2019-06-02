#include <iostream>
#include <stdlib.h>
using namespace std;

void crear(double **&A,int n, int m)
{
    A=new double*[n];
    for(int i=0;i<n;i++)
        A[i]=new double[m];
}

void llenar(double **&A, int n, int m)
{
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            A[i][j]=rand()%10+1;
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


void escalona(double **&A,int n,int m)
{
	for(int i=0;i<m;i++)
		for(int j=i+1;j<n;j++)
		{
			double num=A[j][i];
			double denom=A[i][i];
			for(int k=i;k<m;k++)			
				A[j][k]=A[j][k]-((A[i][k]/denom)*num );
		}
}

void maxIndex(double **&A,int n, int m,int irow ,int icol  ,int &indexMax)
{
	indexMax=irow;
	for(int i=irow;i<n;i++)
		if(A[i][icol] > A[indexMax][icol])		
			indexMax=i;		
	
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


void escalonaPiv(double **&A,int n,int m)
{
	for(int i=0;i<m;i++)
	{
		int indexMax;
		maxIndex(A,n,m,i,i,indexMax);
		pivoteo(A,n,m,indexMax,i,i);
		for(int j=i+1;j<n;j++)
		{
			double num=A[j][i];
			double denom=A[i][i];
			for(int k=i;k<m;k++)			
				A[j][k]=A[j][k]-((A[i][k]/denom)*num );
		}	
	}
}


void sustRegre(double **&A, int n, int m, double *&res)
{

	double suma=0;
	res[n-1]=A[n-1][m-1]/A[n-1][n-1];
	
	for(int z=n-2;z>=0;z--)
	{
		suma=0;
		for(int y=z+1;y<n;y++)
		{
			suma=suma+A[z][y]*res[y];
		}
		res[z]=(A[z][m-1]-suma)/A[z][z];
	}
}


void metGauss(double **&A,int n, int m,double *res)
{
    cout<<endl<<"Matriz escalonada"<<endl;
	escalona(A,n,m);
	imprimir(A,n,m);   
	sustRegre(A,n,m,res);
}

void metGaussPiv(double **&A,int n, int m,double *res)
{
    cout<<endl<<"Matriz escalonada con pivoteo"<<endl;
	escalonaPiv(A,n,m);
	imprimir(A,n,m);   
	sustRegre(A,n,m,res);
}

int main()
{
    double **A;
    int n=4,m=5;
 	double *res = new double[m-1];
	crear(A,n,m);
    llenar(A,n,m);
    cout<<"Matriz sin escalonar"<<endl;
    imprimir(A,n,m);   

    metGauss(A,n,m,res);
 /*
 	metGaussPiv(A,n,m,res);
 */	   
 
	cout<<endl<<"La respuesta del sistema de ecuaciones es : "<<endl;
   	for(int i=0;i<m-1;i++)
		cout<<res[i]<<" ";
	cout<<endl;
	
    
    return 0;
}
