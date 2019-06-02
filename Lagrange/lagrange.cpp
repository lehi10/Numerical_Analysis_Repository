#include <iostream>
#include <vector>

using namespace std;

vector<double> convolution(vector<double> v1, vector<double> v2)
{
  int v1Sz=v1.size();
  int v2Sz=v2.size();
  int resGrado=v1Sz-1+v2Sz-1;
  vector<double> res(resGrado+1);
  int itTemp=resGrado;
  for(int i=0;i<v2Sz;i++)
  {
    itTemp=i;
    for(int j=0;j<v1Sz;j++)
    {
        res[itTemp]+=v2[j]*v1[i];
        itTemp++;
    }
  }
  return res;
}




double langrange(vector<vector<double> > puntos , int nPuntos,int dim, double X)
{
  double res=0;
  for(int i=0;i<nPuntos;i++)
  {

    double numerador=1;
    double denominador=1;

    for(int j=0;j<nPuntos;j++)
    {
      if(i!=j)
      {
        numerador*=(X- puntos[j][0] );
        denominador*=(puntos[i][0]-puntos[j][0]);
      }
    }

    double multiplicatoria=numerador/denominador;

    res+=puntos[i][1]*multiplicatoria;
  }
return res;
}

int main()
{
int numPuntos=4;
int dim=2;
vector<vector<double> >  puntos(numPuntos);

for(int i=0;i<numPuntos;i++)
  puntos[i]=vector<double>(2);

puntos[0][0]=1;
puntos[0][1]=3;

puntos[1][0]=2;
puntos[1][1]=14;

puntos[2][0]=0;
puntos[2][1]=2;

puntos[3][0]=-1;
puntos[3][1]=-1;

/*
vector<double> v1(2);
v1[0]=1;
v1[1]=-2;
vector<double> v2(2);
v2[0]=1;
v2[1]=2;
vector<double> res=convolution(v1,v2);
for(int i=0;i<res.size();i++)
  cout<<res[i]<<" ";
*/

cout<<langrange(puntos,numPuntos,2, 3)<<endl;


return 0;
}
