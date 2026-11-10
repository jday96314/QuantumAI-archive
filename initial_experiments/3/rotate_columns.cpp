#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
int main(int argc,char**argv){
 if(argc!=6)return 2;int C=stoi(argv[3]),B=stoi(argv[4]),steps=stoi(argv[5]);
 vector<double>A(C*B),norm(C);ifstream in(argv[1],ios::binary);in.read((char*)A.data(),A.size()*8);
 auto l1=[&](int c){double s=0;for(int k=0;k<B;k++)s+=abs(A[c*B+k]);return s;};
 for(int c=0;c<C;c++)norm[c]=l1(c);
 uint64_t seed=19260817;auto rnd=[&](){seed^=seed<<13;seed^=seed>>7;seed^=seed<<17;return seed;};
 for(int t=0;t<steps;t++){
  int i=rnd()%C,j=rnd()%(C-1);if(j>=i)j++;
  double*a=&A[i*B],*b=&A[j*B];double old=norm[i]*norm[i]+norm[j]*norm[j];
  double best=old,theta=0,na=norm[i],nb=norm[j];
  auto check=[&](double angle){double c=cos(angle),s=sin(angle),x=0,y=0;
   for(int k=0;k<B;k++){x+=abs(c*a[k]+s*b[k]);y+=abs(c*b[k]-s*a[k]);}
   double score=x*x+y*y;if(score<best){best=score;theta=angle;na=x;nb=y;}
  };
  const double pi=acos(-1.);for(int k=1;k<16;k++)check(k*pi/32);
  double center=theta;for(int k=-7;k<=7;k++)check(center+k*pi/256);
  center=theta;for(int k=-7;k<=7;k++)check(center+k*pi/2048);
  if(best<old-1e-12){double c=cos(theta),s=sin(theta);for(int k=0;k<B;k++){double x=a[k],y=b[k];a[k]=c*x+s*y;b[k]=c*y-s*x;}norm[i]=na;norm[j]=nb;}
  if(t%50000==0){double cost=0;for(double x:norm)cost+=x*x;cerr<<t<<" lambda_h2 "<<.25*cost<<'\n';}
 }
 ofstream out(argv[2],ios::binary);out.write((char*)A.data(),A.size()*8);
 double cost=0;for(double x:norm)cost+=x*x;cout<<.25*cost<<'\n';
}
