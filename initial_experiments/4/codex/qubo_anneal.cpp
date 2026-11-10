#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>
// Minimize b.x + sum_{i<j} H_ij x_i x_j. Diagonal H is ignored.
extern "C" double solve_qubo(int n,const double* H,const double* b,
                            int restarts,int cap,uint64_t seed,int* answer){
 auto rnd=[&](){seed^=seed<<13;seed^=seed>>7;seed^=seed<<17;return seed;};
 auto uniform=[&](){return (rnd()>>11)*0x1.0p-53;};
 std::vector<double> scale(n);for(int i=0;i<n;i++)scale[i]=std::abs(b[i]);
 std::nth_element(scale.begin(),scale.begin()+n/2,scale.end());
 double temp0=std::max(1e-18,scale[n/2]*2),best=0;
 std::fill(answer,answer+n,0);
 for(int run=0;run<restarts;run++){
  std::vector<int>x(n,0);std::vector<double>field(b,b+n);double value=0;int count=0;
  auto flip=[&](int i){int sign=x[i]?-1:1;value+=sign*field[i];
   count+=sign;x[i]=1-x[i];for(int j=0;j<n;j++)if(i!=j)field[j]+=sign*H[j*n+i];
   if(value<best){best=value;std::copy(x.begin(),x.end(),answer);}};
  for(int stage=0;stage<48;stage++){
   double temp=temp0*std::pow(1e-4,stage/47.0);
   for(int t=0;t<8*n;t++){
    int i=rnd()%n;if(!x[i]&&count>=cap)continue;
    double delta=(x[i]?-1:1)*field[i];
    if(delta<=0 || uniform()<std::exp(-delta/temp))flip(i);
   }
  }
  for(int t=0;t<4*n;t++){
   int choice=-1;double gain=-1e-18;
   for(int i=0;i<n;i++)if(x[i]||count<cap){double v=(x[i]?-1:1)*field[i];if(v<gain){gain=v;choice=i;}}
   if(choice<0)break;flip(choice);
  }
 }
 return best;
}
