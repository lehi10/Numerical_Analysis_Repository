#include <iostream>
#include <stdlib.h>
#include <math.h>

using namespace std;



/**
 * Reserva espacio de memoria para una matriz de n x m 
 */
double **crear(int n, int m){
    double **A=new double*[n];
    for(int i=0;i<n;i++)
        A[i]=new double[m];
    return A;
}

/**
 * Reserva espacio de Memoria para un array de tamaño n
 */
double *crear(int n){
    double *vec=new double[n];
    return vec;
}

/**
 * Llenar Manualmente(desde consola) una matriz
 */
void llenarManual(double **&A, int n, int m){
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            cin>>A[i][j];
}

/**
 * Llenar Manualmente(desde consola) un Array
 */
void llenarManual(double *&A, int n){
    for(int i=0;i<n;i++)
            cin>>A[i];
}

/**
 * Llenar Matriz
 */
void llenar(double **&A, int n, int m){
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            A[i][j]=rand()%10+1;
}

/**
 * Llenar Array
 */
void llenar(double *&b, int n){
    for(int i=0;i<n;i++)
        b[i]=rand()%30+1;
}

/**
 * Imprimir contenido de una matriz
 **/
void imprimir(double **&A, int n, int m){
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++)
            cout<<A[i][j]<<" ";
        cout<<endl;
    }
}

/**
 * Imprimir contenido de un array
 **/
void imprimir(double *&b,int n)
{
    for(int i=0;i<n;i++)    
        cout<<b[i]<<" ";
    cout<<endl;
}

// RESOLUCION DE SISTEMA DE ECUACIONES POR MEDIO DE PLU

void identidad(double **&A, int n, int m){
    for(int i=0;i<n;i++)
        A[i][i]=1;
}

void copiar(double **&A, double **&B,int n, int m){
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            B[i][j]=A[i][j];
}

double ** productoMat(double **&A, int n,int m,double**&B,int h,int w){
    double **result=crear(n,w);
    for(int f=0;f<n;f++)
        for(int c=0;c<w;c++)
            for(int i=0;i<m;i++)
                result[f][c]+=A[f][i]*B[i][c];
    return result;
}

double * productoMatVec(double **&A, int n,int m,double*&b,int w){
    double *result=crear(n);    
    for(int f=0;f<n;f++)
        for(int i=0;i<m;i++)
            result[f]+=A[f][i]*b[i];
    
    return result;
        
}

double **sumaMat(double **&A, double**&B ,int n,int m){
    double **result=crear(n,m);
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            result[i][j]=A[i][j]+B[i][j];
    return result;
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


void restaVec(double *&res ,double *& v1,double *&v2 ,int &sz)
{
    for(int i=0;i<sz;i++)
        res[i]=v1[i]-v2[i];
}

void jacobiano2(double **&res ,double (*JF[][2])(double*) , int JF_f,int JF_c,double *X ,int X_sz)
{
    for(int i=0;i<JF_f;i++)
        for(int j=0;j<JF_c;j++)
            res[i][j]=(*JF[i][j])(X);
}

void jacobiano3(double **&res ,double (*JF[][3])(double*) , int JF_f,int JF_c,double *X ,int X_sz)
{
    for(int i=0;i<JF_f;i++)
        for(int j=0;j<JF_c;j++)
            res[i][j]=(*JF[i][j])(X);
}



void evalFunciones(double *res, double (*F[])(double*) , int F_sz, double *X, int X_sz )
{
    for(int i=0;i<F_sz;i++)
        res[i]=(*F[i])(X);
}

void copiarV(double *&v1,double *&v2,int size)
{
    for(int i=0;i<size;i++)
        v1[i]=v2[i];
}

double distancia(double *&v1,double *&v2,int size)
{
    double res=0;
    for(int i=0;i<size;i++)
        res+=fabs(v1[i]-v2[i]);
    return res;
}



double * metNewtonSist2( double (*F[])(double*) ,int F_sz,double (*JF[][2])(double *),int JF_f,int JF_c ,double *X0 ,int X0_sz, double tol,int iter ){

    double *Xi =X0;
    double **A  =crear(JF_f,JF_c);                    
    double *Fx=crear(F_sz);
    double *Xant=crear(X0_sz);

    //do
    for(int i=0;i<iter;i++)
    {
        copiarV(Xant,Xi,X0_sz);
        jacobiano2(A,JF,JF_f,JF_c,Xant,X0_sz);
        evalFunciones(Fx,F,F_sz,Xant,X0_sz);
        double *_X =ResolverSistConPLU(A,Fx,JF_f,JF_c);
        restaVec(Xi,Xant,_X,X0_sz);
    }
    //while(distancia(Xi,Xant,X0_sz) > tol );

    return Xi;
}


double * metNewtonSist3( double (*F[])(double*) ,int F_sz,double (*JF[][3])(double *),int JF_f,int JF_c ,double *X0 ,int X0_sz, double tol ){

    double *Xi =X0;
    double **A  =crear(JF_f,JF_c);                    
    double *Fx=crear(F_sz);
    double *Xant=crear(X0_sz);

    do
    {
        copiarV(Xant,Xi,X0_sz);
        jacobiano3(A,JF,JF_f,JF_c,Xant,X0_sz);
        evalFunciones(Fx,F,F_sz,Xant,X0_sz);
        double *_X =ResolverSistConPLU(A,Fx,JF_f,JF_c);
        restaVec(Xi,Xant,_X,X0_sz);
    }
    while(distancia(Xi,Xant,X0_sz) > tol );

    return Xi;
}




double * metNewtonSistJaprox( double (*F[])(double*) ,int F_sz,double *X0 ,int X0_sz, double tol ,int iter){

    double *Xi =X0;
    int JF_f=F_sz,JF_c=X0_sz;
    double **A  =crear(JF_f,JF_c);                    
    double *Fx=crear(F_sz);
    double *Xant=crear(X0_sz);

    //do
    for(int i=0;i<iter;i++)
    {
        copiarV(Xant,Xi,X0_sz);
        //jacobianoAprox(A,JF_f,JF_c,Xant,X0_sz);
        
        evalFunciones(Fx,F,F_sz,Xant,X0_sz);
        double *_X =ResolverSistConPLU(A,Fx,JF_f,JF_c);
        restaVec(Xi,Xant,_X,X0_sz);
    }
    //while(distancia(Xi,Xant,X0_sz) > tol );
    
    return Xi;
}


//Funciones de F
double f1(double *x){return pow(x[0],2)+pow(x[1],2)-4 ;}
double f2(double *x){return x[0]-pow(x[1],2)-1 ;}

//Derivadas parciales de las funciones de F
double df1(double *x){return 2*x[0];}
double df2(double *x){return 2*x[1];}
double df3(double *x){return 1;}
double df4(double *x){return -2*x[1];}


/*
//Funciones de F
double f1(double *x){return pow(x[0],2)+pow(x[1],2)+pow(x[2],2)-2;}
double f2(double *x){return x[0]+x[1]+x[2]-1;}
double f3(double *x){return pow(x[0],2)+pow(x[1],2)-x[2];}

//Derivadas parciales de las funciones de F
double df1(double *x){return 2*x[0];}
double df2(double *x){return 2*x[1];}
double df3(double *x){return 2*x[2];}
double df4(double *x){return 1;}
double df5(double *x){return 1;}
double df6(double *x){return 1;}
double df7(double *x){return 2*x[0];}
double df8(double *x){return 2*x[1];}
double df9(double *x){return -1;}
*/


int main()
{
//Con 2 variables

    int F_sz=2;
    int JF_f=2 ,JF_c=2;
    int X0_sz=2;
    double tol=0.00001;
    
    double *X0=crear(X0_sz);
    X0[0]=1;
    X0[1]=1;
    double (*F[2])(double*)={f1,f2};
    double (*JF[2][2])(double*)={{df1,df2},{df3,df4}};
    double *X=metNewtonSist2(F,F_sz,JF,JF_f,JF_c,X0,X0_sz,tol,5);
    imprimir(X,X0_sz); 


/*
//Con 3 variables

    int F_sz=3;
    int JF_f=3 ,JF_c=3;
    int X0_sz=3;
    double tol=0.00001;
    
    double *X0=crear(X0_sz);
    X0[0]=1;
    X0[1]=1;
    X0[1]=0;

    //Funciones F
    double (*F[F_sz])(double*)={f1,f2,f3};
    //Matriz Jacobiana de las funciones de F
    double (*JF[3][3])(double*)=   {{df1,df2,df3},
                                    {df4,df5,df6},
                                    {df7,df8,df9}};

    double *X=metNewtonSist3(F,F_sz,JF,JF_f,JF_c,X0,X0_sz,tol);
    
    
    imprimir(X,X0_sz);

   */
   return 0;
}
