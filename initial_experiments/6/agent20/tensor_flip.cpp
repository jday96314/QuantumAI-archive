// Exact GF(2) rank-one tensor flips and plus moves. Local research only.
#include <array>
#include <vector>
#include <random>
#include <fstream>
#include <iostream>
#include <chrono>
#include <algorithm>
using T=std::array<unsigned short,3>; using V=std::vector<T>;
std::mt19937_64 rng(2026102501); size_t pick(size_t n){return rng()%n;}
void remove(V &v,size_t k){v[k]=v.back();v.pop_back();}
void reduce(V &v){bool again=true;while(again){again=false;for(size_t i=0;i<v.size();i++)if(!v[i][0]||!v[i][1]||!v[i][2]){remove(v,i);again=true;break;}if(again)continue;for(size_t i=0;i<v.size()&&!again;i++)for(size_t j=i+1;j<v.size();j++){int same=0,axis=-1;for(int k=0;k<3;k++)if(v[i][k]==v[j][k])same++;else axis=k;if(same>=2){if(axis<0)axis=0;v[i][axis]^=v[j][axis];remove(v,j);again=true;break;}}}}
struct F{int i,j,k;};
std::vector<F> flips(const V &v){std::vector<F> f;for(size_t i=0;i<v.size();i++)for(size_t j=i+1;j<v.size();j++)for(int k=0;k<3;k++)if(v[i][k]==v[j][k])f.push_back({(int)i,(int)j,k});return f;}
void flip(V &v,F f){int a=f.i,b=f.j;if(pick(2))std::swap(a,b);int j=(f.k+1)%3,k=(f.k+2)%3;if(pick(2))std::swap(j,k);v[a][j]^=v[b][j];v[b][k]^=v[a][k];reduce(v);}
void plus(V &v){int i=pick(v.size()),j=pick(v.size()-1);if(j>=i)j++;int x=pick(3),y=(x+1)%3,z=(x+2)%3;if(pick(2))std::swap(y,z);T a=v[i],b=v[j],p=a,q=b,r=b;p[x]^=b[x];q[y]^=a[y];r[y]=a[y];r[z]^=a[z];v[i]=p;v[j]=q;v.push_back(r);reduce(v);}
std::array<unsigned short,144> tensor(const V &v){std::array<unsigned short,144> t{};for(auto r:v)for(int i=0;i<12;i++)if((r[0]>>i)&1)for(int j=0;j<12;j++)if((r[1]>>j)&1)t[12*i+j]^=r[2];return t;}
void save(const V &v,const char *path){std::ofstream f(path);for(auto t:v)f<<t[0]<<' '<<t[1]<<' '<<t[2]<<'\n';}
int main(int argc,char **argv){int seconds=argc>1?std::stoi(argv[1]):120;std::ifstream in("parameters/mul42_tensor.txt");V initial;unsigned a,b,c;while(in>>a>>b>>c)initial.push_back(T{(unsigned short)a,(unsigned short)b,(unsigned short)c});const auto want=tensor(initial);V best=initial,v=initial;auto start=std::chrono::steady_clock::now();uint64_t nflip=0,nplus=0,restarts=0,steps=0;int age=0;while(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<seconds){steps++;auto f=flips(v);if((f.empty()||pick(100)<3)&&v.size()<best.size()+6){plus(v);nplus++;}else if(!f.empty()){flip(v,f[pick(f.size())]);nflip++;}else age=10000;if(v.size()<best.size()){if(tensor(v)!=want){std::cerr<<"INVALID TENSOR\n";return 2;}best=v;save(best,"parameters/mul_tensor_best.txt");std::cout<<"best "<<best.size()<<" steps "<<steps<<std::endl;age=0;}if(++age>=10000){v=best;age=0;restarts++;}}if(tensor(v)!=want||tensor(best)!=want)return 2;save(best,"parameters/mul_tensor_best.txt");std::cout<<"{\"status\":\"EXACT_TENSOR_SEARCH_COMPLETE\",\"initial_rank\":"<<initial.size()<<",\"best_rank\":"<<best.size()<<",\"steps\":"<<steps<<",\"flips\":"<<nflip<<",\"plus\":"<<nplus<<",\"restarts\":"<<restarts<<",\"seconds\":"<<seconds<<",\"seed\":2026102501}"<<std::endl;}
