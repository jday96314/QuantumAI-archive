// Fit a cubic table against the exact radix-four row rounding model.
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
struct Fit {double error;array<long long,4> c;};
int main(int argc,char**argv){
 if(argc!=9)return 2;int f=stoi(argv[1]),lf=stoi(argv[2]);array<int,4>w={stoi(argv[3]),stoi(argv[4]),stoi(argv[5]),stoi(argv[6])};
 ifstream in(argv[7]);int count;in>>count;vector<array<long long,4>> initial(count);for(auto &c:initial)for(auto &v:c)in>>v;
 int L=21-f,span=1<<L,leaves=1<<(lf-f),unit=span/leaves;vector<array<long long,4>> output;double worst=0;int failed=0;
 for(int addr=0;addr<count;addr++)for(int leaf=0;leaf<leaves;leaf++){
  int lo=leaf*unit,hi=lo+unit,base=((1<<f)+(addr% (1<<f)))*span;double scale=addr>=(1<<f)?sqrt(2.):1.;
  auto eval=[&](array<long long,4> coeff,bool full){
   int stride=full?1:max(1,unit/256);vector<long long> values;vector<double> exact;double minbias=1e100,maxbias=-1e100;
   for(int offset=lo;offset<hi;offset+=stride){
    int nn=base+offset;bool valid=addr>=(1<<f)?nn%2==0:(nn<=3133452||nn%4==0);if(!valid)continue;
    array<long long,4> c=coeff;int z=offset^(1<<(L-1));
    for(int j=2;j>=0;j--)for(int g=0;g<(L+1)/2;g++){
     auto bit=[&](int k){return k<0?0:(z>>min(k,L-1))&1;};int digit=bit(2*g-1)+bit(2*g)-2*bit(2*g+1);if(!digit)continue;
     int e=w[j]-w[j+1]-21+2*g+(abs(digit)==2);long long v=e>=0?c[j+1]<<e:(c[j+1]>>-e)+((c[j+1]>>(-e-1))&1);c[j]+=digit<0?v:-v;
    }
    double target=ldexp(scale,w[0]+10)/sqrt(nn);values.push_back(c[0]);exact.push_back(target);minbias=min(minbias,target-c[0]);maxbias=max(maxbias,target-c[0]);
   }
   long long center=llround((minbias+maxbias)/2),bestbias=0;double best=1e100;
   for(long long bias=center-2;bias<=center+2;bias++){
    double err=0;for(size_t i=0;i<values.size();i++)err=max(err,abs((values[i]+bias)/exact[i]-1));if(err<best){best=err;bestbias=bias;}
   }
   coeff[0]+=bestbias;return Fit{best,coeff};
  };
  vector<Fit> options;auto c=initial[addr];for(int d1=-12;d1<=12;d1++)for(int d2=-3;d2<=3;d2++)for(int d3=-2;d3<=2;d3++){
   auto trial=c;trial[1]+=d1;trial[2]+=d2;trial[3]+=d3;if(*min_element(trial.begin(),trial.end())<0)continue;options.push_back(eval(trial,false));
  }
  sort(options.begin(),options.end(),[](auto a,auto b){return a.error<b.error;});Fit best{1e100,{}};
  for(size_t i=0;i<min(size_t(12),options.size());i++){auto v=eval(options[i].c,true);if(v.error<best.error)best=v;}
  worst=max(worst,best.error);failed+=best.error>=ldexp(1.,-24);output.push_back(best.c);
 }
 ofstream out(argv[8]);out<<output.size()<<'\n';for(auto c:output){for(auto v:c)out<<v<<' ';out<<'\n';}
 cout.precision(17);cout<<"{\"worst_relative_error\":"<<worst<<",\"failed_leaves\":"<<failed<<",\"leaves\":"<<output.size()<<"}"<<endl;
}
