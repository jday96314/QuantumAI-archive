// Fit a dyadic coefficient table while keeping one common local coordinate.
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <functional>
using namespace std;
struct Fit{double err;long long a,b,c;};
int main(int argc,char**argv){
 if(argc!=9)return 2;int f=stoi(argv[1]),fm=stoi(argv[2]),wc=stoi(argv[3]),w1=stoi(argv[4]),w2=stoi(argv[5]);
 string option=argv[8];bool booth=option.find("booth")==0;double tolerance=ldexp(1.,-24);
 auto colon=option.find(':');if(colon!=string::npos)tolerance=stod(option.substr(colon+1));
 ifstream in(argv[6]);int count;in>>count;vector<array<long long,3>> init(count);for(auto&r:init)in>>r[0]>>r[1]>>r[2];
 int span=1<<(21-f),unit=1<<(21-fm),rep=1<<(fm-f),leaves=0;double worst=0;bool pass=true;
 vector<array<long long,3>> output(1<<(fm+1));
 for(int a=0;a<count;a++){
  auto old=init[a];int base=((1<<f)+(a% (1<<f)))*span;double scale=a>=(1<<f)?sqrt(2.):1.;
  auto eval=[&](int offset,int end,long long slope,long long curve,bool full){
   int stride=full?1:max(1,(end-offset)/256),len=(end-offset+stride-1)/stride;
   vector<long long>s(len,slope),y(len,old[0]);vector<double> ideal(len);vector<bool> valid(len);
   if(booth){
    int L=21-f;
    for(int g=0;g<(L+1)/2;g++)for(int t=0;t<len;t++){
     int z=(offset+t*stride)^(1<<(L-1));auto bit=[&](int k){return k<0?0:(z>>min(k,L-1))&1;};int d=bit(2*g-1)+bit(2*g)-2*bit(2*g+1);if(!d)continue;
     int e=w1-w2-21+2*g+(abs(d)==2);long long v=e>=0?curve<<e:(curve>>-e)+((curve>>(-e-1))&1);s[t]+=d<0?v:-v;
    }
    for(int g=0;g<(L+1)/2;g++)for(int t=0;t<len;t++){
     int z=(offset+t*stride)^(1<<(L-1));auto bit=[&](int k){return k<0?0:(z>>min(k,L-1))&1;};int d=bit(2*g-1)+bit(2*g)-2*bit(2*g+1);if(!d)continue;
     int e=wc-w1-21+2*g+(abs(d)==2);long long v=e>=0?s[t]<<e:(s[t]>>-e)+((s[t]>>(-e-1))&1);y[t]+=d<0?v:-v;
    }
   }else{
   for(int m=f+1;m<=21;m++){
    int e=w1-w2-m;long long v=e>=0?curve<<e:(curve>>-e)+((curve>>(-e-1))&1);
    for(int t=0;t<len;t++){int bit=((offset+t*stride)>>(21-m))&1;s[t]+=(m==f+1?(1-bit):-bit)*v;}
   }
   for(int m=f+1;m<=21;m++){
    int e=wc-w1-m;
    for(int t=0;t<len;t++){long long v=e>=0?s[t]<<e:(s[t]>>-e)+((s[t]>>(-e-1))&1);int bit=((offset+t*stride)>>(21-m))&1;y[t]+=(m==f+1?(1-bit):-bit)*v;}
   }
   }
   double lo=1e100,hi=-1e100;
   for(int t=0;t<len;t++){int nn=base+offset+t*stride;valid[t]=a>=(1<<f)?nn%2==0:(nn<=3133452||nn%4==0);ideal[t]=ldexp(scale,wc+10)/sqrt(nn);if(valid[t]){lo=min(lo,ideal[t]-y[t]);hi=max(hi,ideal[t]-y[t]);}}
   long long mid=llround((lo+hi)/2),bias=mid;double best=1e100;
   for(long long b=mid-2;b<=mid+2;b++){double err=0;for(int t=0;t<len;t++)if(valid[t]){if(y[t]+b<0||y[t]+b>=(1LL<<wc)){err=1e100;break;}err=max(err,abs((y[t]+b)/ideal[t]-1));}if(err<best){best=err;bias=b;}}
   return Fit{best,old[0]+bias,slope,curve};
  };
  function<void(int,int)> visit=[&](int begin,int end){
   vector<Fit> opts;for(int db=-18;db<=18;db++)for(int dc=-3;dc<=3;dc++)if(old[1]+db>=0&&old[2]+dc>=0)opts.push_back(eval(begin,end,old[1]+db,old[2]+dc,false));
   sort(opts.begin(),opts.end(),[](Fit x,Fit y){return x.err<y.err;});Fit best{1e100,0,0,0};
   for(int k=0;k<min(12,(int)opts.size());k++){auto x=eval(begin,end,opts[k].b,opts[k].c,true);if(x.err<best.err)best=x;}
   if(option.find("deep")!=string::npos&&best.err>=tolerance&&end-begin<=unit){
    for(int k=12;k<(int)opts.size();k++){auto x=eval(begin,end,opts[k].b,opts[k].c,true);if(x.err<best.err)best=x;if(best.err<tolerance)break;}
   }
   if(best.err>=tolerance&&end-begin>unit){int mid=(begin+end)/2;visit(begin,mid);visit(mid,end);return;}
   leaves++;worst=max(worst,best.err);if(best.err>=tolerance)pass=false;
   for(int k=begin/unit;k<end/unit;k++)output[a*rep+k]={best.a,best.b,best.c};
  };
  visit(0,span);
  if(!pass)break;
 }
 ofstream out(argv[7]);out<<output.size()<<'\n';for(auto r:output)out<<r[0]<<' '<<r[1]<<' '<<r[2]<<'\n';
 cout.precision(17);cout<<worst<<' '<<leaves<<' '<<pass<<'\n';
}
