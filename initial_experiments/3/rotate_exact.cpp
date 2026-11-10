// Exact pair minimization of sum of squared column L1 norms.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
struct Event {double x;int row,kind;bool operator<(const Event&o)const{return x<o.x;}};
int main(int argc,char**argv){
 if(argc!=6)return 2;int C=stoi(argv[3]),B=stoi(argv[4]),steps=stoi(argv[5]);
 vector<double>A(C*B),norm(C);ifstream in(argv[1],ios::binary);in.read((char*)A.data(),A.size()*8);
 auto l1=[&](int c){double s=0;for(int k=0;k<B;k++)s+=abs(A[c*B+k]);return s;};
 for(int c=0;c<C;c++)norm[c]=l1(c);
 uint64_t seed=19260817;auto rnd=[&](){seed^=seed<<13;seed^=seed>>7;seed^=seed<<17;return seed;};
 const double pi=acos(-1.),lo=-pi/4,hi=pi/4;vector<Event> ev;ev.reserve(B+2);
 for(int t=0;t<steps;t++){
  int i=rnd()%C,j=rnd()%(C-1);if(j>=i)j++;
  double*a=&A[i*B],*b=&A[j*B];double old=norm[i]*norm[i]+norm[j]*norm[j];
  ev.clear();double U=0,V=0,P=0,Q=0,cs=cos(lo),sn=sin(lo);
  for(int k=0;k<B;k++){
   double sa=cs*a[k]+sn*b[k]>=0?1:-1,sb=cs*b[k]-sn*a[k]>=0?1:-1;
   U+=sa*a[k];V+=sa*b[k];P+=sb*b[k];Q-=sb*a[k];
   double x=atan2(-a[k],b[k]);while(x<lo)x+=pi;while(x>hi)x-=pi;
   if(x>lo&&x<hi)ev.push_back({x,k,0});
   x=atan2(b[k],a[k]);while(x<lo)x+=pi;while(x>hi)x-=pi;
   if(x>lo&&x<hi)ev.push_back({x,k,1});
  }
  sort(ev.begin(),ev.end());ev.push_back({hi,-1,0});
  double best=old,theta=0,left=lo;
  auto check=[&](double x){double c=cos(x),s=sin(x),p=U*c+V*s,q=P*c+Q*s,val=p*p+q*q;
   if(val<best){best=val;theta=x;}};
  for(auto e:ev){
   check(left);check(e.x);
   double aa=U*U+P*P,bb=V*V+Q*Q,cc=U*V+P*Q;
   double opt=.5*atan2(2*cc,aa-bb)+pi/2;while(opt>hi)opt-=pi;
   if(opt>=left&&opt<=e.x)check(opt);
   if(e.row>=0){int k=e.row;double c=cos(e.x),s=sin(e.x);
    if(e.kind==0){double ds=(-a[k]*s+b[k]*c)>=0?2:-2;U+=ds*a[k];V+=ds*b[k];}
    else{double ds=(-b[k]*s-a[k]*c)>=0?2:-2;P+=ds*b[k];Q-=ds*a[k];}
   }left=e.x;
  }
  if(best<old-1e-12){double c=cos(theta),s=sin(theta),na=0,nb=0;
   for(int k=0;k<B;k++){na+=abs(c*a[k]+s*b[k]);nb+=abs(c*b[k]-s*a[k]);}
   if(na*na+nb*nb<old-1e-12){for(int k=0;k<B;k++){double x=a[k],y=b[k];a[k]=c*x+s*y;b[k]=c*y-s*x;}norm[i]=na;norm[j]=nb;}
  }
  if(t%50000==0){double cost=0;for(double x:norm)cost+=x*x;cerr<<t<<" lambda_h2 "<<.25*cost<<'\n';}
 }
 ofstream out(argv[2],ios::binary);out.write((char*)A.data(),A.size()*8);
 double cost=0;for(double x:norm)cost+=x*x;cout.precision(16);cout<<.25*cost<<'\n';
}
