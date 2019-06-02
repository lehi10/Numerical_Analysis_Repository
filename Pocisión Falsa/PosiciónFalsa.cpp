#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#define _USE_MATH_DEFINES

using namespace std; 

double funcion1(double x)
{
    return 10*pow(x,3)+x*atan(x+2*M_PI)+5;
}

double funcion2(double x)
{
    return  pow(x,2)*cos(2*x-1)+exp(x-2);
}



float posicionFalsa(double (*ptr_fun)(double),double begin , double end, int  it_num, string fileName)
{
    
    ofstream  ofs;
    ofs.open (fileName+".csv", std::ofstream::out);
    
    double a=begin, b=end;

    //double m= (ptr_fun(b)*a - ptr_fun(a)*b)/(ptr_fun(b) - ptr_fun(a));
    double m = b-ptr_fun(b) *( (b-a) / ( ptr_fun(b) - ptr_fun(a) ) );
    int i=0;

    ofs<<i<<" , "<< a << ", "<<b<<" , "<< m << "," <<ptr_fun(m) <<endl;
    
    for( i=1;i<it_num;i++ )
	{
        if((ptr_fun(a)*ptr_fun(m)) < 0  )
            b=m;
        else
            a=m;
        
        ofs<<i<<" , "<< a << ", "<<b<<" , "<< m << "," <<ptr_fun(m) <<endl;
        m= (ptr_fun(b)*a - ptr_fun(a)*b)/(ptr_fun(b) - ptr_fun(a));
	}
}

float posicionFalsa2(double (*ptr_fun)(double) ,double begin, double end, double  tol, string fileName)
{
    
    ofstream  ofs;
    ofs.open (fileName+".csv", std::ofstream::out);
    
    int k=0;
    double a=begin, b=end;
    double m= (ptr_fun(b)*a - ptr_fun(a)*b)/(ptr_fun(b) - ptr_fun(a));


    ofs<<k<<" , "<< a << ", "<<b<<" , "<< m << "," <<ptr_fun(m) <<endl;
    k++;
	while( fabs(ptr_fun(m)) >= tol )
	{
        if((ptr_fun(a)*ptr_fun(m)) < 0  )
            b=m;
        else
            a=m;
        
        ofs<<k<<" , "<< a << ", "<<b<<" , "<< m << "," <<ptr_fun(m) <<endl;
        m= (ptr_fun(b)*a - ptr_fun(a)*b)/(ptr_fun(b) - ptr_fun(a));
        k++;
	}
}

float posicionFalsa3(double (*ptr_fun)(double),double begin, double end, double  tol, string fileName)
{
    
    ofstream  ofs;
    ofs.open (fileName+".csv", std::ofstream::out);
    
    double a=begin, b=end;
    int k=0;
    double m= (ptr_fun(b)*a - ptr_fun(a)*b)/(ptr_fun(b) - ptr_fun(a));

    ofs<<k<<" , "<< a << ", "<<b<<" , "<< m << "," <<fabs(a-b) << "," <<ptr_fun(m) <<endl;
    k++;
	while( fabs(a-b)  >= tol )
	{
        if((ptr_fun(a)*ptr_fun(m)) < 0  )
            b=m;
        else
            a=m;
        
        ofs<<k<<" , "<< a << ", "<<b<<" , "<< m << "," <<fabs(a-b) << "," <<ptr_fun(m) <<endl;
        m= (ptr_fun(b)*a - ptr_fun(a)*b)/(ptr_fun(b) - ptr_fun(a));
        k++;
	}
}


int main()
{
    //posicionFalsa(funcion1,-2,5,100,"resultados posicionFalsa 1 ecuacion 1");
    cout<<"terminó"<<endl;
    //posicionFalsa2(funcion1,-2,5,0.01,"resultados posicionFalsa 2 ecuacion 1");
    cout<<"terminó2"<<endl;
    //posicionFalsa3(funcion1,-2,5,0.1,"resultados posicionFalsa 3 ecuacion 1");
    cout<<"terminó3"<<endl;
    posicionFalsa(funcion2,-2,5,50,"resultados posicionFalsa 1 ecuacion 2");
    cout<<"terminó4"<<endl;
    posicionFalsa2(funcion2,-2,5,0.0001,"resultados posicionFalsa 2 ecuacion 2");
    cout<<"terminó5"<<endl;
    posicionFalsa3(funcion2,-2,5,0.0001,"resultados posicionFalsa 3 ecuacion 2");
    cout<<"terminó6"<<endl;
    
	return 0;
}