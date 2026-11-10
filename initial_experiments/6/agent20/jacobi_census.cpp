// Independent exact GF(4096) model of the emitted blocked Jacobi map.
// Constructed RS errors are used only to check the resulting rational pair.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>
using V=std::vector<int>;
std::array<int,8190> ex;
std::array<int,4096> lg;
int mul(int a,int b){return a&&b?ex[lg[a]+lg[b]]:0;}
int inv(int a){return a?ex[4095-lg[a]]:0;}
int power(int a,int n){int r=1;while(n){if(n&1)r=mul(r,a);a=mul(a,a);n>>=1;}return r;}
V product(const V&a,const V&b,int n){V c(n);for(int i=0;i<(int)a.size()&&i<n;i++)for(int j=0;j<(int)b.size()&&i+j<n;j++)c[i+j]^=mul(a[i],b[j]);return c;}
V divide(const V&a,const V&b,int n){V c(n);int iv=inv(b[0]);for(int i=0;i<n;i++){int x=i<(int)a.size()?a[i]:0;for(int j=1;j<=i&&j<(int)b.size();j++)x^=mul(b[j],c[i-j]);c[i]=mul(x,iv);}return c;}
void pair(V&w,int offset){int a=w[offset],ai=inv(a?a:1),n=w.size()-offset;V f(n);f[0]=1;for(int i=1;i<n;i++)f[i]=mul(w[offset+i],ai);V q=divide({1},f,n);for(int i=1;i<n;i++)w[offset+i]=q[i];}
void block(V&w,int offset){int n=w.size()-offset;if(n==2){pair(w,offset);return;}V h(w.begin()+offset,w.begin()+offset+4);pair(h,0);pair(h,2);int a=h[0],b=h[1],c=h[2],d=h[3];int e=b^d,eta=mul(b,d)^c,p=mul(a,c);if(!p)p=1;V A={1,e,eta},D={1,b};int P0=mul(a,power(b,4)^mul(c,c)^mul(c,mul(e,e)));int P1=mul(eta,mul(a,power(b,3)^mul(c,d)));V r(w.begin()+offset+4,w.end());int len=r.size();if(len){r=product(A,r,len);r[0]^=P0;if(len>1)r[1]^=P1;V denominator={p,0};V dr=product(D,r,std::max(len-2,0));denominator.insert(denominator.end(),dr.begin(),dr.end());r=divide(product(A,r,len),denominator,len);}std::copy(h.begin(),h.end(),w.begin()+offset);std::copy(r.begin(),r.end(),w.begin()+offset+4);}
std::pair<V,V> fold(const V&w,int offset){int a=w[offset],b=w[offset+1];if(offset+2==(int)w.size())return {{a},{b}};auto [p,q]=fold(w,offset+2);int t=p.size();V P={a},Q={b^q[0]};for(int i=0;i<t;i++){P.push_back(mul(a?a:1,q[i]));Q.push_back(p[i]^mul(b,q[i])^(i+1<t?q[i+1]:0));}return {P,Q};}
int eval(const V&p,int x){int v=0;for(auto i=p.rbegin();i!=p.rend();++i)v=mul(v,x)^*i;return v;}
void lookahead_block(V&w,int offset){int n=w.size()-offset;if(n==2){pair(w,offset);return;}V h(w.begin()+offset,w.begin()+offset+4);bool gap=!h[0];if(gap){int u=h[1],e=mul(h[2],inv(u?u:1)),eta=mul(h[3],inv(u?u:1))^mul(e,e);h={0,u,eta,e};}else{pair(h,0);pair(h,2);h[2]^=mul(h[1],h[3]);h[3]^=h[1];}int a=h[0],b=h[1],eta=h[2],e=h[3],bg=gap?b:0,c=eta^mul(b,e)^mul(b,b),d=e^b,p=mul(a,c)^bg;if(!p)p=1;int p0=mul(a,power(b,4)^mul(c,c)^mul(c,mul(e,e)))^mul(bg,power(e,3));int p1=mul(eta,mul(a,power(b,3)^mul(c,d))^mul(bg,eta^mul(e,e)));V r(w.begin()+offset+4,w.end());int len=r.size(),shift=gap?3:2;V A={1,e,eta},D={1,b^bg};if(len){r=product(A,r,len);r[0]^=p0;if(len>1)r[1]^=p1;V denominator(shift);denominator[0]=p;V dr=product(D,r,std::max(len-shift,0));denominator.insert(denominator.end(),dr.begin(),dr.end());r=divide(product(A,r,len),denominator,len);}std::copy(h.begin(),h.end(),w.begin()+offset);std::copy(r.begin(),r.end(),w.begin()+offset+4);}
std::pair<V,V> fold_one(int a,int b,const V&p,const V&q){V P={a},Q={b^(q.empty()?0:q[0])};for(int i=0;i<(int)p.size();i++){P.push_back(mul(a?a:1,q[i]));Q.push_back(p[i]^mul(b,q[i])^(i+1<(int)q.size()?q[i+1]:0));}return {P,Q};}
std::pair<V,V> lookahead_fold(const V&w,int offset){if(offset+2==(int)w.size())return {{w[offset]},{w[offset+1]}};V p,q;if(offset+4<(int)w.size())std::tie(p,q)=lookahead_fold(w,offset+4);int a=w[offset],b=w[offset+1],h=w[offset+2],e=w[offset+3];if(a){int d=e^b,c=h^mul(b,d);auto mid=fold_one(c,d,p,q);return fold_one(a,b,mid.first,mid.second);}auto get=[&](int i){return i<(int)q.size()?q[i]:0;};V P={0,b},Q={e^get(0),h^mul(e,get(0))^get(1)};for(int i=0;i<(int)p.size();i++){P.push_back(mul(b?b:1,q[i]));Q.push_back(p[i]^mul(h,get(i))^mul(e,get(i+1))^get(i+2));}return {P,Q};}
int main(int argc,char**argv){int shots=argc>1?std::stoi(argv[1]):20000;uint64_t seed=argc>2?std::stoull(argv[2]):2026102701;lg.fill(-1);int x=1;for(int i=0;i<4095;i++){ex[i]=x;lg[x]=i;x<<=1;if(x&4096)x^=0x1053;}for(int i=4095;i<8190;i++)ex[i]=ex[i-4095];
 std::ifstream weights("parameters/dqi_weight_probabilities.txt");std::vector<double> probabilities;double v;while(weights>>v)probabilities.push_back(v);if(probabilities.size()!=36)return 2;
 bool repair=argc>3&&std::string(argv[3])=="lookahead";std::mt19937_64 rng(seed);std::discrete_distribution<int> weight(probabilities.begin(),probabilities.end());std::uniform_int_distribution<int> field(1,4095),position(0,4094);int wrong=0,serial_wrong=0;std::array<int,36> totals{},bad{};
 for(int shot=0;shot<shots;shot++){int d=weight(rng);totals[d]++;V pos,values,s(70);while((int)pos.size()<d){int p=position(rng);if(std::find(pos.begin(),pos.end(),p)==pos.end())pos.push_back(p);}for(int p:pos){int a=ex[p],e=field(rng);values.push_back(e);int pp=a;for(int i=0;i<70;i++){s[i]^=mul(e,pp);pp=mul(pp,a);}}
  auto check=[&](bool blocked){V w=s;if(blocked){for(int i=0;i<70;i+=4)(repair?lookahead_block:block)(w,i);}else for(int i=0;i<70;i+=2)pair(w,i);auto [P,Q]=repair&&blocked?lookahead_fold(w,0):fold(w,0);Q.insert(Q.begin(),1);int degree=Q.size()-1;while(degree&&Q[degree]==0)degree--;if(degree!=d)return false;V derivative(degree);for(int i=1;i<=degree;i+=2)derivative[i-1]=Q[i];for(int i=0;i<d;i++){int beta=ex[(4095-pos[i])%4095];if(eval(Q,beta))return false;int denominator=eval(derivative,beta);if(!denominator||mul(eval(P,beta),inv(denominator))!=values[i])return false;}return true;};
  if(!check(true)){wrong++;bad[d]++;}if(!check(false))serial_wrong++;
 }
 std::cout<<"{\"shots\":"<<shots<<",\"seed\":"<<seed<<",\"wrong_magnitude_vectors\":"<<wrong<<",\"serial_wrong_magnitude_vectors\":"<<serial_wrong<<",\"by_degree\":[";for(int d=0;d<36;d++){if(d)std::cout<<',';std::cout<<"{\"d\":"<<d<<",\"total\":"<<totals[d]<<",\"wrong\":"<<bad[d]<<'}';}std::cout<<"]}"<<std::endl;
}
