// Independent GF(4096) reference census; no verifier samples or circuit hashes.
#include <algorithm>
#include <cmath>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>
using P=std::vector<int>;
int ex[8190],lg[4096];
int raw(int a,int b){int r=0;while(b){if(b&1)r^=a;b>>=1;a<<=1;if(a&4096)a^=0x1053;}return r;}
int mul(int a,int b){return a&&b?ex[lg[a]+lg[b]]:0;}
int inv(int a){return a?ex[4095-lg[a]]:0;}
int power(int a,int k){return a?ex[(int64_t(lg[a])*k)%4095]:0;}
void trim(P& p){while(!p.empty()&&!p.back())p.pop_back();}
P rem(P a,const P& b){int bi=inv(b.back());while(a.size()>=b.size()){int d=a.size()-b.size(),q=mul(a.back(),bi);for(size_t j=0;j<b.size();j++)a[d+j]^=mul(q,b[j]);trim(a);}return a;}
int eval(const P& p,int x){int v=0;for(auto i=p.rbegin();i!=p.rend();++i)v=mul(v,x)^*i;return v;}
P product(const P& a,const P& b){P c(a.size()+b.size()-1);for(size_t i=0;i<a.size();i++)for(size_t j=0;j<b.size();j++)c[i+j]^=mul(a[i],b[j]);trim(c);return c;}
P armed_sylvester(P f,P g,P c,P d){
 int m=f.size()-1,n=g.size()-1;if(!g[n])g[n]=1;
 if(!n){for(int& x:d)x=mul(x,g[0]);return d;}
 P w=f;int gi=inv(g[n]);for(int i=m;i>=n;i--){w[i]=mul(w[i],gi);for(int j=0;j<n;j++)w[i-n+j]^=mul(w[i],g[j]);}
 P q(w.begin()+n,w.end()),r(w.begin(),w.begin()+n);int k=n-1;
 for(size_t i=0;i<q.size();i++)for(size_t j=0;j<c.size();j++)d[i+j]^=mul(q[i],c[j]);
 P low=armed_sylvester(g,r,P(d.begin(),d.begin()+k),c);low.resize(m+n);
 for(int i=k;i<m;i++)for(int j=0;j<=n;j++)low[i+j]^=mul(d[i],g[j]);
 return low;
}
P hprepare(const P& A,const P& B){
 P C,D,E,O;for(int i=0;i<35;i++)(i%2?D:C).push_back(power(i<int(A.size())?A[i]:0,2048));
 for(int i=0;i<=35;i++)(i%2?O:E).push_back(power(B[i],2048));
 P h=armed_sylvester(E,O,D,P(C.begin(),C.end()-1));h.push_back(C.back());
 for(int j=0;j<17;j++)h[17+j]^=mul(C.back(),O[j]);return h;
}
bool hchart(const P& B){
 int d=B.size()-1;if(d<=1)return true;
 P f,g;for(int i=0;i<=d;i++){if(i%2)g.push_back(power(B[i],2048));else f.push_back(power(B[i],2048));}
 if(d%2==0 && !g.back())return false;
 // Keep the first declared f degree even if its top is zero: division by
 // monic g is still valid and the generic first remainder must drop one.
 while(g.size()>1){P r=rem(f,g);if(r.size()!=g.size()-1)return false;f=g;g=r;}
 return !g.empty()&&g[0];
}
std::pair<P,P> decode(P W,bool& supported){
 const int N=W.size();supported=true;
 for(int k=0;k<N;k+=2){
  int a=W[k],u=W[k+1],L=N-k-2,ui=inv(u?u:1);
  bool nonzero=false,tail_nonzero=false;for(int j=k;j<N;j++)nonzero|=W[j]!=0;for(int j=k+2;j<N;j++)tail_nonzero|=W[j]!=0;
  if((!a&&nonzero)||(!u&&tail_nonzero))supported=false;
  P x(L+1),y(L+1);x[0]=y[0]=1;
  for(int i=1;i<=L;i++)x[i]=mul(W[k+1+i],ui);
  for(int i=1;i<=L;i++){int v=x[i];for(int j=1;j<i;j++)v^=mul(x[j],y[i-j]);y[i]=v;}
  int b=mul(u,inv(a?a:1));W[k+1]=b;
  for(int i=0;i<L;i++)W[k+2+i]=y[i+1]^(i==0?b:0);
 }
 P p{W[N-2]},q{W[N-1]};
 for(int k=N-4;k>=0;k-=2){int a=W[k],b=W[k+1],L=p.size();P r(L),np(L+1),nq(L+1);for(int i=0;i<L;i++)r[i]=p[i]^q[i];np[0]=a;nq[0]=b^r[0];for(int i=0;i<L;i++){np[i+1]=mul(r[i],a?a:1);nq[i+1]=mul(q[i],b?b:1)^(i+1<L?r[i+1]:0);}p=np;q=nq;}
 P B{1};B.insert(B.end(),q.begin(),q.end());trim(B);trim(p);return {B,p};
}
int main(int argc,char**argv){
 int shots=argc>1?std::stoi(argv[1]):10000;uint64_t seed=argc>2?std::stoull(argv[2]):20261024071ULL;
 int a=1;for(int i=0;i<4095;i++){ex[i]=a;lg[a]=i;a=raw(a,2);}if(a!=1)return 2;for(int i=4095;i<8190;i++)ex[i]=ex[i-4095];
 std::vector<double> weights(36);std::ifstream wf("parameters/host_weights.txt");double sum=0;for(double& w:weights){wf>>w;sum+=w;}if(!wf||std::abs(sum-1)>1e-12)return 3;
 std::mt19937_64 rng(seed),mrng(seed^0x9234567abcULL);std::discrete_distribution<int> wd(weights.begin(),weights.end());
 P quadratics(17);std::ifstream qf("parameters/host_quadratics.txt");for(int& q:quadratics)qf>>q;if(!qf)return 4;
 int mismatches=0,unsupported=0,phase_errors=0,hgaps=0,padded_hgaps=0,h_wrong=0,h_phase=0,chart_bad=0;std::array<int,36> count{},hgap_by_degree{};
 for(int shot=0;shot<shots;shot++){
  int d=wd(rng);count[d]++;std::set<int> set;while(int(set.size())<d)set.insert(rng()%4095);P pos(set.begin(),set.end()),val(d),S(70),B{1};
  for(int i=0;i<d;i++){int aa=ex[pos[i]],e=1+rng()%4095;val[i]=e;int ap=aa;for(int k=0;k<70;k++){S[k]^=mul(e,ap);ap=mul(ap,aa);}P nB(B.size()+1);for(size_t j=0;j<B.size();j++){nB[j]^=B[j];nB[j+1]^=mul(aa,B[j]);}B=nB;}
  if(!hchart(B)){hgaps++;hgap_by_degree[d]++;}
  P padding{1};if((35-d)%2)padding.insert(padding.begin(),0);for(int k=1;k<35;k+=2)if(k>=d)padding=product(padding,{quadratics[(k-1)/2],1,1});
  P paddedB=product(B,padding);
  if(paddedB.size()!=36)return 5;
  bool hgap=!hchart(paddedB);padded_hgaps+=hgap;
  if(hgap||shot<64){
   int bi=inv(B.back());P bm=B;for(int& x:bm)x=mul(x,bi);P A(d),H(d);
   for(int i=0;i<d;i++){
    int beta=ex[(4095-pos[i])%4095];P quotient(d);quotient[d-1]=1;for(int j=d-2;j>=0;j--)quotient[j]=bm[j+1]^mul(beta,quotient[j+1]);
    int den=eval(quotient,beta),v=power(beta,35),ac=mul(mul(v,v),inv(mul(val[i],mul(den,den)))),hc=mul(v,inv(mul(power(val[i],2048),den)));
    for(int j=0;j<d;j++){A[j]^=mul(ac,quotient[j]);H[j]^=mul(hc,quotient[j]);}
   }
   P ap=product(A,padding),bp=product(bm,padding),expected=product(H,padding),hh=hprepare(ap,bp);expected.resize(35);
   if(!hgap&&hh!=expected){std::cerr<<"H supported mismatch\n";return 6;}
   // Literal inverse padding network, then cyclic monomial unpadding.
   for(int k=33;k>=1;k-=2)if(k>=d)for(int i=34;i>=std::max(0,35-k-2);i--){if(i+1<35)hh[i]^=hh[i+1];if(i+2<35)hh[i]^=mul(quadratics[(k-1)/2],hh[i+2]);}
   std::rotate(hh.begin(),hh.begin()+(35-d),hh.end());bool bad=false;int ph=0;
   for(int i=0;i<d;i++){int beta=ex[(4095-pos[i])%4095],hv=eval(hh,beta),got=mul(power(beta,70),inv(mul(hv,hv)));bad|=got!=val[i];ph^=__builtin_parity(unsigned((got^val[i])&(mrng()%4096)));}
   h_wrong+=bad;h_phase+=ph;
  }
  P omega(35);for(size_t i=0;i<B.size();i++)for(int j=0;j<35-int(i);j++)omega[i+j]^=mul(B[i],S[j]);trim(omega);
  bool supported;auto [b,p]=decode(S,supported);unsupported+=!supported;
  bool bad=b!=B||p!=omega;if(supported&&bad)chart_bad++;if(!bad)continue;mismatches++;
  P deriv(b.size()-1);for(size_t i=1;i<b.size();i+=2)deriv[i-1]=b[i];trim(deriv);P truth(4095);for(int i=0;i<d;i++)truth[pos[i]]=val[i];int phase=0;
  for(int j=0;j<4095;j++){int beta=ex[(4095-j)%4095],e=0;if(!eval(b,beta))e=mul(eval(p,beta),inv(eval(deriv,beta)));int mask=mrng()%4096;phase^=__builtin_parity(unsigned((e^truth[j])&mask));}phase_errors+=phase;
 }
 std::cout<<"{\"status\":\""<<(chart_bad?"FAIL":"PASS_CLASSICAL_MODEL")<<"\",\"shots\":"<<shots<<",\"seed\":"<<seed<<",\"cf_unsupported\":"<<unsupported<<",\"cf_wrong_pairs\":"<<mismatches<<",\"cf_wrong_supported_pairs\":"<<chart_bad<<",\"cf_random_mask_phase_errors\":"<<phase_errors<<",\"h_generic_chart_failures\":"<<hgaps<<",\"h_squarefree_padded_chart_failures\":"<<padded_hgaps<<",\"h_wrong_vectors\":"<<h_wrong<<",\"h_random_mask_phase_errors\":"<<h_phase<<",\"degree_counts\":[";
 for(int i=0;i<36;i++)std::cout<<(i?",":"")<<count[i];std::cout<<"],\"h_chart_failures_by_degree\":[";for(int i=0;i<36;i++)std::cout<<(i?",":"")<<hgap_by_degree[i];std::cout<<"],\"scope\":\"independent classical host-law census, not native verification or accepted confidence score\"}\n";
 return chart_bad?1:0;
}
