// Independent exact signed radix-four row-rounded integer Horner fitter.
// No circuit builder or previous table is used. Floating values only seed and
// rank fits; certify.py subsequently checks integer phase against isqrt targets.
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <vector>
using namespace std;
using I=long long;
struct Fit {double err=1e99; I a=0,b=0,c=0;};
struct Leaf {int p,lo,depth,origin,degree; Fit fit;};
I rounded(I a,int e) {return e>=0?a<<e:(a>>-e)+((a>>(-e-1))&1);}
I product(I a,int z,int L,int e) {
 I out=0; auto bit=[&](int k){return k<0?0:(z>>min(k,L-1))&1;};
 for(int g=0;g<(L+1)/2;g++) {
  int d=bit(2*g-1)+bit(2*g)-2*bit(2*g+1);
  if(d)out+=(d>0?1:-1)*rounded(a,e+2*g+(abs(d)==2));
 } return out;
}
int main(int argc,char**argv) {
 if(argc!=10){cerr<<"fit mode base wc w1 w2 tolerance maxdepth output radius\n";return 2;}
 string mode=argv[1];int base=stoi(argv[2]),wc=stoi(argv[3]),w1=stoi(argv[4]),w2=stoi(argv[5]);
 double tol=stod(argv[6]);int maxdepth=stoi(argv[7]);ofstream out(argv[8]);
 int radius=stoi(argv[9]);vector<Leaf> leaves;double worst=0;bool ok=true;
 int mixed_cut=1<<20;auto at=mode.find('@');
 if(at!=string::npos)mixed_cut=llround(stod(mode.substr(at+1))*(1<<21));
 auto fit=[&](int parity,int lo,int depth,int origin,int degree) {
  int size=1<<(21-depth),L=21-origin;
  double center=1.+(double((lo>>(21-origin))<<(21-origin))+(1<<(20-origin)))/(1<<21);
  double low=1.+double(lo)/(1<<21),high=low+double(size)/(1<<21),scale=parity?1.:sqrt(.5);
  // Fit around the shared origin, including when the leaf is asymmetric.
  double mid=(low+high)/2, h=(high-low)/2;
  double c0=scale/sqrt(mid),c1=.5*scale/pow(mid,1.5),c2=degree==2?.375*scale/pow(mid,2.5):0.;
  // Chebyshev correction for interpolation, then shift polynomial origin.
  if(degree==1)c0+=.5*c2*h*h;
  else c0+=.5*c2*h*h;
  double delta=center-mid;
  c0-=c1*delta;c0+=c2*delta*delta;c1-=2*c2*delta;
  I a=llround(ldexp(c0,wc)),b=llround(ldexp(c1,w1)),c=llround(ldexp(c2,w2));
  auto evaluate=[&](I bb,I cc,bool full) {
   int stride=full?1:max(1,size/192);vector<pair<I,double>> values;
   double mn=1e99,mx=-1e99;double factor=ldexp(scale,wc);
   for(int off=0;off<size;off+=stride) {
    int n=(1<<21)+lo+off;
    if(parity ? n%2!=0 : (n>3133452&&n%4!=0))continue;
    int z=(n&((1<<L)-1))^(1<<(L-1));
    I slope=bb-(degree==2?product(cc,z,L,w1-w2-21):0);
    if(slope<0)return Fit{};
    I y=a-product(slope,z,L,wc-w1-21);
    double ideal=factor/sqrt(double(n)/(1<<21));
    mn=min(mn,ideal-y);mx=max(mx,ideal-y);values.push_back({y,ideal});
   }
   Fit best;
   if(values.empty())return Fit{0,a,bb,cc};
   I bias=llround((mn+mx)/2);
   for(I shift=bias-2;shift<=bias+2;shift++) {
    double err=0;for(auto [y,ideal]:values) {
     if(y+shift<0||y+shift>=(1LL<<wc)){err=1e99;break;}
     err=max(err,abs((y+shift)/ideal-1));
    }
    if(err<best.err)best={err,a+shift,bb,cc};
   }return best;
  };
  vector<Fit> fits;
  for(int db=-radius;db<=radius;db++)for(int dc=degree==2?-3:0;dc<=(degree==2?3:0);dc++)
   if(b+db>=0&&c+dc>=0)fits.push_back(evaluate(b+db,c+dc,false));
  sort(fits.begin(),fits.end(),[](Fit a,Fit b){return a.err<b.err;});Fit best;
  for(int j=0;j<min(8,int(fits.size()));j++) {
   auto f=evaluate(fits[j].b,fits[j].c,true);if(f.err<best.err)best=f;
  }return best;
 };
 for(int p=0;p<2;p++) {
  function<void(int,int)> visit=[&](int lo,int depth) {
   int origin=mode.find("shared")==0?base:depth;
   bool linear=mode.find("linear")!=string::npos || (mode.find("mixed")!=string::npos&&lo>=mixed_cut);
   int degree=linear?1:2;auto f=fit(p,lo,depth,origin,degree);
   if(f.err>tol&&depth<maxdepth) {
    visit(lo,depth+1);visit(lo+(1<<(20-depth)),depth+1);return;
   }
   ok&=f.err<tol;worst=max(worst,f.err);leaves.push_back({p,lo,depth,origin,degree,f});
  };
  for(int lo=0;lo<(1<<21);lo+=1<<(21-base))visit(lo,base);
 }
 out.precision(17);out<<"{\"mode\":\""<<mode<<"\",\"base\":"<<base<<",\"widths\":["<<wc<<","<<w1<<","<<w2<<"],\"tolerance\":"<<tol<<",\"max_depth\":"<<maxdepth<<",\"radius\":"<<radius<<",\"sample_points\":192,\"full_candidates\":8,\"mixed_cut\":"<<mixed_cut<<",\"leaves\":[\n";
 for(int j=0;j<int(leaves.size());j++) {
  auto l=leaves[j];out<<"["<<l.p<<","<<l.lo<<","<<l.depth<<","<<l.origin<<","<<l.degree<<","<<l.fit.a<<","<<l.fit.b<<","<<l.fit.c<<","<<l.fit.err<<"]"<<(j+1==int(leaves.size())?"\n":",\n");
 }
 out<<"],\"pass\":"<<(ok?"true":"false")<<",\"max_relative_error\":"<<worst<<"}\n";
 cout.precision(17);cout<<mode<<" "<<wc<<","<<w1<<","<<w2<<" leaves "<<leaves.size()<<" error "<<worst<<" pass "<<ok<<endl;
}
