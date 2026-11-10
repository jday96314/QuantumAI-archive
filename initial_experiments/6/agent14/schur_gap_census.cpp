// Independent fixed-law census of the quadratic-node Schur representation.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
using namespace std;
uint16_t ex[12285],lg[4096];unsigned K;
int mul(int a,int b){return a&&b?ex[lg[a]+lg[b]]:0;}
int inv(int a){return a?ex[4095-lg[a]]:0;}
unsigned scale(unsigned v,int c){return mul(v&4095,c)|(mul(v>>12,c)<<12);}
unsigned iext(unsigned v){int a=v&4095,b=v>>12,n=mul(a,a)^mul(a,b)^mul(K,mul(b,b)),iv=inv(n);return mul(a^b,iv)|(mul(b,iv)<<12);}
using Poly=vector<int>;
void trim(Poly &a){while(a.size()&&!a.back())a.pop_back();}
Poly add(Poly a,const Poly&b){a.resize(max(a.size(),b.size()));for(int i=0;i<(int)b.size();i++)a[i]^=b[i];trim(a);return a;}
Poly linear(const Poly&a,int u,int v){if(a.empty())return {};Poly z(a.size()+1);for(int i=0;i<(int)a.size();i++){z[i]^=mul(u,a[i]);z[i+1]^=mul(v,a[i]);}trim(z);return z;}
Poly times(Poly a,int c){for(auto &v:a)v=mul(v,c);trim(a);return a;}
int eval(const Poly&a,int x){int v=0;for(int i=(int)a.size()-1;i>=0;i--)v=mul(v,x)^a[i];return v;}
int deriv(const Poly&a,int x){int v=0;for(int i=(int)a.size()-2;i>=0;i--)v=mul(v,x)^((i%2==0)?a[i+1]:0);return v;}
template<class T>void read(ifstream&f,T*p,size_t n){f.read(reinterpret_cast<char*>(p),n*sizeof(T));}
struct Step{int a,b,k,kind,olda;};
Poly quadratic(const Poly&p,int k){return add(linear(linear(p,0,1),0,1),linear(p,k,1));}
int main(int argc,char**argv){if(argc!=4)return 2;ifstream f("parameters/schur_census.bin",ios::binary);uint32_t N;read(f,&N,1);read(f,&K,1);if(N!=70)return 3;double prob[36];read(f,prob,36);read(f,ex,4095);read(f,lg,4096);for(int i=4095;i<12285;i++)ex[i]=ex[i%4095];uint16_t kk[35];uint32_t maps[35][70];read(f,kk,35);read(f,&maps[0][0],2450);uint16_t pp[71];read(f,pp,71);Poly P(pp,pp+71);long long count=stoll(argv[1]);unsigned long long seed=stoull(argv[2]);bool swap=stoi(argv[3]);mt19937_64 rng(seed);discrete_distribution<int> law(prob,prob+36);long long badpivot=0,poles=0,false_complete=0,wrong=0,magwrong=0,swaps=0,weight_sum=0,passed=0;array<long long,36> weights{};
for(long long shot=0;shot<count;shot++){
 int d=law(rng);weights[d]++;weight_sum+=d;vector<int> pos,val;array<bool,4096>used{};for(int i=0;i<d;i++){int a;do{a=ex[rng()%4095];}while(used[a]);used[a]=true;pos.push_back(a);val.push_back(1+rng()%4095);}
 int S[70]={};for(int i=0;i<d;i++){int v=val[i];for(int j=0;j<70;j++){v=mul(v,pos[i]);S[j]^=v;}}
 unsigned ys[35];int nodes[35],kap[35];for(int i=0;i<35;i++){nodes[i]=2*i;kap[i]=kk[i];unsigned v=0;for(int j=0;j<70;j++)v^=scale(maps[i][j],S[j]);ys[i]=v;}
 int g=1,A=0,B=0;vector<Step>stages;bool failed=false;
 for(int i=0;i<35;i++){
  if(i==0){A=B=0;}else if(!stages.back().b){A=B=0;}else if(i>1&&!stages[i-2].b){A=stages[i-2].a;B=0;}else{A=stages.back().a^stages.back().b;B=stages.back().b;}
  auto pv=[&](int j){int b=ys[j]>>12,a=(ys[j]&4095)^mul(nodes[j],b);return array<int,4>{a,b,mul(g,a)^A,mul(g,b)^B};};
  auto v=pv(i);
  if(!v[3]&&!v[2]){bool z=true;for(int j=i;j<35;j++)if(scale(ys[j],g)^A^(mul(nodes[j],B)|(B<<12)))z=false;if(!z){false_complete++;failed=true;}break;}
  for(int j=i+1;j<35;j++){unsigned den=ys[j]^v[0]^(mul(nodes[j],v[1])|(v[1]<<12));if(!den){poles++;failed=true;break;}ys[j]=scale(iext(den),kap[i]^kap[j]);}if(failed)break;
  bool prior_zero=i&&!stages.back().b;int kind=prior_zero?2:!v[3]?1:0;int olda=i?stages.back().a:0;
  if(kind==2&&!v[3]){badpivot++;failed=true;break;}
  if(kind==1)swaps++;stages.push_back({v[2],v[3],kap[i],kind,olda});g=mul(prior_zero||!v[3]?1:mul(v[3],v[3]),inv(g));
 }
 if(!failed&&!stages.empty()&&stages.back().kind==1){badpivot++;failed=true;}
 if(failed)continue;
 Poly num,den{1};for(int j=(int)stages.size()-1;j>=0;j--){auto t=stages[j];Poly nn,bb;
  if(!t.kind){int h=inv(t.b),r=mul(t.a,h),c=mul(mul(t.a,t.a)^mul(t.a,t.b)^mul(t.k,mul(t.b,t.b)),h);nn=add(linear(num,r,1),times(den,c));bb=add(times(num,h),linear(den,r^1,1));}
  else if(t.kind==1){nn=add(times(num,t.a),quadratic(den,t.k));bb=num;}
  else{nn=add(linear(num,t.a,t.b),add(quadratic(den,t.k),linear(den,mul(t.olda,t.a),mul(t.olda,t.b))));bb=add(num,times(den,t.olda));}
  num=nn;den=bb;
 }
 int sc=inv(den.back());num=times(num,sc);den=times(den,sc);Poly truth{1};for(int a:pos)truth=linear(truth,a,1);if(den!=truth){wrong++;continue;}
 bool mbad=false;for(int i=0;i<d;i++){int res=mul(eval(num,pos[i]),inv(deriv(den,pos[i])));int value=mul(res,inv(mul(pos[i],eval(P,pos[i]))));if(value!=val[i]){mbad=true;break;}}
 if(mbad)magwrong++;else passed++;
 if((shot+1)%10000==0)cout<<"{\"done\":"<<shot+1<<",\"passed\":"<<passed<<",\"pivot_exceptions\":"<<badpivot<<",\"extension_poles\":"<<poles<<",\"scalar_charts\":"<<swaps<<"}"<<endl;
}
cout<<"{\"status\":\"CENSUS_COMPLETE\",\"cases\":"<<count<<",\"seed\":"<<seed<<",\"unused_cli_swap\":"<<(swap?"true":"false")<<",\"passed\":"<<passed<<",\"unsupported_scalar_charts\":"<<badpivot<<",\"extension_poles\":"<<poles<<",\"false_completion\":"<<false_complete<<",\"locator_errors_on_supported\":"<<wrong<<",\"magnitude_errors_on_supported\":"<<magwrong<<",\"scalar_charts\":"<<swaps<<",\"mean_weight\":"<<double(weight_sum)/count<<"}"<<endl;
}
