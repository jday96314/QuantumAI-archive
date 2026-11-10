#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
#include <array>
#include <atomic>
#include <omp.h>
using namespace std;
struct Candidate{double error;int s,c,b;};
int main(int argc,char**argv){
 if(argc!=9)return 2;
 string mode=argv[1];bool binary=mode.find("binary")==0,centered=mode.find("centered")!=string::npos;
 bool reachable=mode.find("reachable")!=string::npos;
 bool deep=mode.find("deep")!=string::npos;
 bool rounded=mode.find("rounded")!=string::npos;
 bool fixedcurve=mode.find("fixedcurve")!=string::npos;
 bool allerrors=mode.find("allerrors")!=string::npos;
 bool booth=mode.find("booth")!=string::npos;
 auto tp=mode.find("ties");int ties=tp==string::npos?0:stoi(mode.substr(tp+4));bool oddties=mode.find("oddties")!=string::npos;
 bool masked=mode.find("roundmask")!=string::npos;unsigned mask0=~0u,mask1=~0u;
 if(masked){string masks=argv[8];auto comma=masks.find(',');mask0=stoul(masks.substr(0,comma));mask1=stoul(masks.substr(comma+1));}
 auto fp=mode.find("fraction");int extra=fp==string::npos?0:stoi(mode.substr(fp+8));int fine=1<<extra;
 bool global=mode.find("globalphase")!=string::npos;double phase=global?stod(argv[8]):0.;
 int f=stoi(argv[2]),wc=stoi(argv[3]),w1=stoi(argv[4]),w2=stoi(argv[5]);
 ifstream in(argv[6]);int count;in>>count;vector<array<int64_t,3>> tab(count);for(auto&r:tab)in>>r[0]>>r[1]>>r[2];
 vector<double> errors(count);int top=binary?21:20;int step=1<<(top-f);
 int first=binary?0:1<<f;
 int bits0=0;for(auto&r:tab)while((1LL<<bits0)<=r[0])bits0++;
 int64_t limit=1LL<<bits0;
 atomic<bool> failed(false);
 #pragma omp parallel for schedule(dynamic)
 for(int a=first;a<count;a++){
  if(!allerrors&&failed.load())continue;
  auto old=tab[a];vector<Candidate> opts;
  auto eval=[&](int slope,int curve,bool full){
   auto increment=[&](int64_t v,int shift){
    int out=(v>>(shift-1))&1;
    if(shift<=ties&&(v&((1LL<<shift)-1))==(1LL<<(shift-1)))out=((v>>shift)&1)^oddties;
    return out;
   };
   int stride=full?1:max(1,step/256);int len=(step+stride-1)/stride;
   vector<int64_t> ss(len,slope),vv(len,old[0]);
   if(booth){
    int L=top-f;
    if(mode.find("boothlast")==string::npos){
    for(int g=0;g<(L+1)/2;g++)for(int t=0;t<len;t++){
     int z=(t*stride)^(1<<(L-1));auto bit=[&](int k){return k<0?0:(z>>min(k,L-1))&1;};
     int d=bit(2*g-1)+bit(2*g)-2*bit(2*g+1);if(!d)continue;
     int e=w1-w2-top+2*g+(abs(d)==2);int64_t v=e>=0?((int64_t)curve<<e):curve>>-e;
     if(rounded&&e<0)v+=increment(curve,-e);ss[t]+=d<0?v:-v;
    }
    }else{
     for(int m=f+1;m<=top;m++){
      int e=w1-w2-m;int64_t v=e>=0?((int64_t)curve<<e):curve>>-e;
      if(rounded&&e<0)v+=increment(curve,-e);
      for(int t=0;t<len;t++){int bit=(t*stride>>(top-m))&1;ss[t]+=(m==f+1?(1-bit):-bit)*v;}
     }
    }
    for(int g=0;g<(L+1)/2;g++)for(int t=0;t<len;t++){
     int z=(t*stride)^(1<<(L-1));auto bit=[&](int k){return k<0?0:(z>>min(k,L-1))&1;};
     int d=bit(2*g-1)+bit(2*g)-2*bit(2*g+1);if(!d)continue;
     int e=wc-w1-top+2*g+(abs(d)==2);int64_t v=e>=0?ss[t]<<e:ss[t]>>-e;
     if(rounded&&e<0)v+=increment(ss[t],-e);vv[t]+=d<0?v:-v;
    }
   }else{
   for(int m=f+1;m<=top;m++){
    int e=w1-w2-m;int64_t v=e>=0?((int64_t)curve<<e):curve>>-e;
    if(rounded&&e<0&&((mask1>>(m-f-1))&1))v+=increment(curve,-e);
    for(int t=0;t<len;t++){
     int bit=(t*stride>>(top-m))&1;
     if(centered&&m==f+1)ss[t]+=(1-bit)*v;else ss[t]-=bit*v;
    }
   }
   for(int m=f+1;m<=top;m++){
    int e=wc-w1-m;
    for(int t=0;t<len;t++){
     int64_t v=e>=0?ss[t]<<e:ss[t]>>-e;int bit=(t*stride>>(top-m))&1;
     if(rounded&&e<0&&((mask0>>(m-f-1))&1))v+=increment(ss[t],-e);
     if(centered&&m==f+1)vv[t]+=(1-bit)*v;else vv[t]-=bit*v;
    }
   }
   }
   double lo=1e100,hi=-1e100;
   vector<double> ideals(len),corrlo(len,0),corrhi(len,0);vector<bool> valid(len,true);
   for(int t=0;t<len;t++){
    double norm=binary?((1<<f)+(a% (1<<f)))*((double)step)+t*stride:((double)a*step)+t*stride;
    if(reachable){int nn=(int)norm;valid[t]=(a>=(1<<f))?(nn%2==0):(nn<=3133452||nn%4==0);}
    if(global&&valid[t]){
     int nn=(int)norm,parity=a>=(1<<f),kmin=parity?1:(nn<=3133452?0:2);
     int kmax=min(21,__builtin_ctz((unsigned)nn));if((kmax&1)!=parity)kmax--;
     corrlo[t]=ldexp(phase,wc+10-kmin/2);corrhi[t]=ldexp(phase,wc+10-kmax/2);
    }
    double scale=binary&&a>=(1<<f)?sqrt(2.):1.;
    double ideal=ldexp(scale,wc+10)/sqrt(norm);ideals[t]=ideal;
    if(valid[t]){lo=min(lo,min(ideal-vv[t]-corrlo[t],ideal-vv[t]-corrhi[t]));hi=max(hi,max(ideal-vv[t]-corrlo[t],ideal-vv[t]-corrhi[t]));}
   }
   int center=llround((lo+hi)*.5*fine),bestb=center;double best=1e100;
   for(int bias=center-2;bias<=center+2;bias++){
    double err=0.;bool overflow=false;
    for(int t=0;t<len;t++){
     if(!valid[t])continue;
     double value=vv[t]+(double)bias/fine;
     if(value>=limit||value<0){overflow=true;break;}
     err=max(err,max(abs((value+corrlo[t])/ideals[t]-1),abs((value+corrhi[t])/ideals[t]-1)));
    }
    if(!overflow&&err<best){best=err;bestb=bias;}
   }
   return Candidate{best,slope,curve,bestb};
  };
  for(int ds=-18;ds<=18;ds++)for(int dc=(fixedcurve?0:-3);dc<=(fixedcurve?0:3);dc++)if(old[1]+ds>=0&&old[2]+dc>=0)
   opts.push_back(eval(old[1]+ds,old[2]+dc,false));
  sort(opts.begin(),opts.end(),[](auto&a,auto&b){return a.error<b.error;});
  Candidate best{1e100,0,0,0};
  for(int i=0;i<(deep?(int)opts.size():min(12,(int)opts.size()));i++){
   if(deep&&i>=12&&opts[i].error>=best.error)break;
   auto x=eval(opts[i].s,opts[i].c,true);if(x.error<best.error)best=x;
  }
  tab[a]={old[0]*fine+best.b,best.s,best.c};errors[a]=best.error;
  if(best.error>=ldexp(1.,-24))failed.store(true);
 }
 ofstream out(argv[7]);out.precision(17);out<<count<<'\n';for(int i=0;i<count;i++){auto&r=tab[i];out<<r[0]<<' '<<r[1]<<' '<<r[2];if(allerrors)out<<' '<<errors[i];out<<'\n';}
 cout.precision(17);cout<<*max_element(errors.begin(),errors.end())<<'\n';
}
