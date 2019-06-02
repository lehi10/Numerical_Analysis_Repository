#include <iostream>
#include <math.h>
#include <fstream>
#define _USE_MATH_DEFINES

using namespace std; 


double funcion1(double x)
{
    return 3*pow(x,5)+2*x*exp(-1*x)+5*cos(x)-4;
}

double funcion2(double x)
{
    return  6* atan((4*x)-5)+exp(-2*x+1)*sin(x+M_PI)-pow(x,2);
}



float biseccion(double (*ptr_fun)(double),double begin , double end, int  it_num, string fileName)
{
    
    ofstream  ofs;
    ofs.open (fileName+".csv", std::ofstream::out);
    
    double a=begin, b=end;
    double m= (a+b)/2;
    int i=0;

    ofs<<i<<" , "<< a << ", "<<b<<" , "<< m << "," <<ptr_fun(m) <<endl;
    
    for( i=1;i<it_num;i++ )
	{
        if((ptr_fun(a)*ptr_fun(m)) < 0  )
            b=m;
        else
            a=m;
        
        ofs<<i<<" , "<< a << ", "<<b<<" , "<< m << "," <<ptr_fun(m) <<endl;
        m=(a+b)/2;
	}
}

float biseccion2(double (*ptr_fun)(double) ,double begin, double end, double  tol, string fileName)
{
    
    ofstream  ofs;
    ofs.open (fileName+".csv", std::ofstream::out);
    
    int k=0;
    double a=begin, b=end;
    double m= (a+b)/2;


    ofs<<k<<" , "<< a << ", "<<b<<" , "<< m << "," <<ptr_fun(m) <<endl;
    k++;
	while( abs(ptr_fun(m)) >= tol )
	{
        if((ptr_fun(a)*ptr_fun(m)) < 0  )
            b=m;
        else
            a=m;
        
        ofs<<k<<" , "<< a << ", "<<b<<" , "<< m << "," <<ptr_fun(m) <<endl;
        m=(a+b)/2;
        k++;
	}
}

float biseccion3(double (*ptr_fun)(double),double begin, double end, double  tol, string fileName)
{
    
    ofstream  ofs;
    ofs.open (fileName+".csv", std::ofstream::out);
    
    double a=begin, b=end;
    int k=0;
    double m= (a+b)/2;

    ofs<<k<<" , "<< a << ", "<<b<<" , "<< m << "," <<abs(a-b) << "," <<ptr_fun(m) <<endl;
    k++;
	while( abs(a-b)  >= tol )
	{
        if((ptr_fun(a)*ptr_fun(m)) < 0  )
            b=m;
        else
            a=m;
        
        ofs<<k<<" , "<< a << ", "<<b<<" , "<< m << "," <<abs(a-b) << "," <<ptr_fun(m) <<endl;
        m=(a+b)/2;
        k++;
	}
}


int main()
{
    biseccion(funcion1,-2,5,20,"resultados biseccion 1 ecuacion 1");
    biseccion2(funcion1,-2,5,0.001,"resultados biseccion 2 ecuacion 1");
    biseccion3(funcion1,-2,5,0.001,"resultados biseccion 3 ecuacion 1");

    biseccion(funcion2,-2,5,20,"resultados biseccion 1 ecuacion 2");
    biseccion2(funcion2,-2,5,0.001,"resultados biseccion 2 ecuacion 2");
    biseccion3(funcion2,-2,5,0.001,"resultados biseccion 3 ecuacion 2");
	return 0;
}