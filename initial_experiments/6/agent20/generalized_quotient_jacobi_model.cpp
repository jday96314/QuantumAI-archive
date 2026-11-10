// Exact variable-degree Jacobi/Euclidean model. No quantum resource claim.
// For valuation v, f=x^v*a/(1+b1*x+...+bL*x^L+x^(L+1)*g), L=v+1.
// Thus f(1/t)/t=a/(t^L+b1*t^(L-1)+...+bL+g(1/t)/t).
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>
using V=std::vector<int>;
std::array<int,8190> ex;std::array<int,4096> lg;
int mul(int a,int b){return a&&b?ex[lg[a]+lg[b]]:0;}
int inv(int a){if(!a)throw std::runtime_error("zero divisor");return ex[4095-lg[a]];}
V prod(const V&a,const V&b){if(a.empty()||b.empty())return {};V c(a.size()+b.size()-1);for(size_t i=0;i<a.size();i++)for(size_t j=0;j<b.size();j++)c[i+j]^=mul(a[i],b[j]);return c;}
void add(V&a,const V&b){a.resize(std::max(a.size(),b.size()));for(size_t i=0;i<b.size();i++)a[i]^=b[i];}
void trim(V&a){while(a.size()>1&&!a.back())a.pop_back();}
int eval(const V&a,int x){int z=0;for(auto i=a.rbegin();i!=a.rend();++i)z=mul(z,x)^*i;return z;}
V series(const V&a,const V&b,int n){V c(n);int iv=inv(b[0]);for(int i=0;i<n;i++){int z=i<(int)a.size()?a[i]:0;for(int j=1;j<=i&&j<(int)b.size();j++)z^=mul(b[j],c[i-j]);c[i]=mul(z,iv);}return c;}
struct Packet{int degree,a;V monic,q;};
std::vector<Packet> extract(V f){std::vector<Packet> packets;int prev=1;while(!f.empty()){int v=0;while(v<(int)f.size()&&!f[v])v++;if(v==(int)f.size())break;int L=v+1;if(2*L>(int)f.size())throw std::runtime_error("not enough rational moments");int a=f[v];V unit(f.begin()+v,f.end());V p=series({a},unit,unit.size());V mono(L+1);mono[L]=1;for(int j=1;j<=L;j++)mono[L-j]=p[j];int k=inv(mul(a,prev));V q=mono;for(int&z:q)z=mul(z,k);packets.push_back({L,a,mono,q});prev=k;f=V(p.begin()+L+1,p.end());}return packets;}
std::pair<V,V> put(const std::vector<Packet>&p,bool reverse=false){V F={1},G;for(int step=0;step<(int)p.size();step++){int j=reverse?step:(int)p.size()-1-step;V next=prod(p[j].q,F);add(next,G);G=F;F=next;}trim(F);trim(G);return {F,G};}
std::pair<V,V> independent_jacobi(const std::vector<Packet>&p){V F={1},G;for(auto i=p.rbegin();i!=p.rend();++i){V nf=prod(i->monic,F);add(nf,G);G=F;for(int&z:G)z=mul(z,i->a);F=nf;}trim(F);trim(G);return {F,G};}
V moments(const V&pos,const V&values){V s(70);for(size_t i=0;i<pos.size();i++){int a=ex[pos[i]],pw=a;for(int j=0;j<70;j++){s[j]^=mul(values[i],pw);pw=mul(pw,a);}}return s;}
// Select the last v magnitudes to force the first v moments to zero.
bool force(V&vals,const V&pos,int v){int d=vals.size();std::vector<V>A(v,V(v+1));for(int r=0;r<v;r++){for(int j=0;j<d;j++){int pw=1;for(int k=0;k<=r;k++)pw=mul(pw,ex[pos[j]]);if(j<d-v)A[r][v]^=mul(vals[j],pw);else A[r][j-(d-v)]=pw;}}for(int c=0;c<v;c++){int r=c;while(r<v&&!A[r][c])r++;if(r==v)return false;std::swap(A[r],A[c]);int iv=inv(A[c][c]);for(int j=c;j<=v;j++)A[c][j]=mul(A[c][j],iv);for(int i=0;i<v;i++)if(i!=c){int f=A[i][c];for(int j=c;j<=v;j++)A[i][j]^=mul(f,A[c][j]);}}for(int i=0;i<v;i++){vals[d-v+i]=A[i][v];if(!vals[d-v+i])return false;}return true;}
int main(int argc,char**argv){int shots=argc>1?std::stoi(argv[1]):100000;uint64_t seed=2026102526;lg.fill(-1);int x=1;for(int i=0;i<4095;i++){ex[i]=x;lg[x]=i;x<<=1;if(x&4096)x^=0x1053;}for(int i=4095;i<8190;i++)ex[i]=ex[i-4095];std::ifstream file("parameters/dqi_weight_probabilities.txt");std::vector<double>prob;double z;while(file>>z)prob.push_back(z);if(prob.size()!=36)return 2;std::mt19937_64 rng(seed);std::discrete_distribution<int>weight(prob.begin(),prob.end());std::uniform_int_distribution<int>field(1,4095),position(0,4094);std::array<int,36>hist{};int wrong=0,forced=0,gapped=0,maxdegree=0;
 auto check=[&](const V&pos,const V&vals,int expected_first){auto s=moments(pos,vals);auto p=extract(s);auto [F,G]=put(p);auto [RF,RG]=put(p,true);auto [JF,JG]=independent_jacobi(p);V c1=prod(G,JF),c2=prod(JG,F);trim(c1);trim(c2);if(c1!=c2)return false;if(F!=RF)return false;if((int)F.size()!=int(pos.size())+1)return false;int total=0;bool gap=false;for(auto&i:p){hist[i.degree]++;total+=i.degree;gap|=i.degree>1;maxdegree=std::max(maxdegree,i.degree);}if(total!=(int)pos.size())return false;if(expected_first&&p[0].degree!=expected_first)return false;gapped+=gap;V D(F.size()-1);for(int i=1;i<(int)F.size();i+=2)D[i-1]=F[i];for(size_t i=0;i<pos.size();i++){int a=ex[pos[i]],der=eval(D,a);if(eval(F,a)||!der)return false;int direct=mul(eval(G,a),inv(mul(a,der)));int endpoint=inv(mul(a,mul(eval(RG,a),der)));if(direct!=vals[i]||endpoint!=vals[i])return false;}return true;};
 for(int shot=0;shot<shots;shot++){int d=weight(rng);V pos,vals;while((int)pos.size()<d){int pp=position(rng);if(std::find(pos.begin(),pos.end(),pp)==pos.end())pos.push_back(pp);}for(int i=0;i<d;i++)vals.push_back(field(rng));if(!check(pos,vals,0))wrong++;}
 for(int d=0;d<=35;d++)for(int trial=0;trial<8;trial++){V pos,vals;while((int)pos.size()<d){int pp=position(rng);if(std::find(pos.begin(),pos.end(),pp)==pos.end())pos.push_back(pp);}for(int i=0;i<d;i++)vals.push_back(field(rng));if(!check(pos,vals,0))wrong++;}
 for(int v=1;v<=4;v++)for(int trial=0;trial<64;trial++){int d=5+trial%31;V pos,vals;while((int)pos.size()<d){int pp=position(rng);if(std::find(pos.begin(),pos.end(),pp)==pos.end())pos.push_back(pp);}for(int i=0;i<d;i++)vals.push_back(field(rng));if(!force(vals,pos,v)){trial--;continue;}if(!check(pos,vals,v+1))wrong++;forced++;}
 std::cout<<"{\"status\":\""<<(wrong?"FAIL":"PASS_CLASSICAL_MODEL")<<"\",\"law_shots\":"<<shots<<",\"seed\":"<<seed<<",\"completion_degree_cases\":288,\"forced_valuation_cases\":"<<forced<<",\"wrong_vectors\":"<<wrong<<",\"vectors_with_gap\":"<<gapped<<",\"maximum_quotient_degree\":"<<maxdegree<<",\"quotient_degree_histogram\":[";for(int i=0;i<36;i++){if(i)std::cout<<',';std::cout<<hist[i];}std::cout<<"],\"scope\":\"variable-degree extraction, independent Jacobi recurrence, reversed quotient cofactor inverse and RS magnitudes; no emitted quantum circuit or resource claim\"}"<<std::endl;return wrong?1:0;
}
