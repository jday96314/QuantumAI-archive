// Quadratic (not restricted bilinear) XOR-AND decomposition search.
// Each AND may mix both 12-bit operands. Products are represented modulo
// affine terms by (a wedge b) tensor c. GL(2,2) changes inside each plane
// are therefore free; final affine differences are recovered exactly.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
using T=std::array<uint32_t,3>;using V=std::vector<T>;
std::mt19937_64 rng(2026102616);size_t pick(size_t n){return rng()%n;}
void canonical(T&t){uint32_t a=t[0],b=t[1];t[0]=std::min({a,b,a^b});t[1]=std::min(a==t[0]?b:a,b==t[0]?a:b);if(t[1]==t[0])t[1]=a^b;}
bool zero(T t){return !t[0]||!t[1]||t[0]==t[1]||!t[2];}
void remove(V&v,size_t i){v[i]=v.back();v.pop_back();}
uint32_t common(T a,T b){for(uint32_t x:{a[0],a[1],a[0]^a[1]})for(uint32_t y:{b[0],b[1],b[0]^b[1]})if(x==y)return x;return 0;}
void orient(T&t,uint32_t a){uint32_t b=t[0]==a?t[1]:t[0];t[0]=a;t[1]=b;}
void random_basis(T&t){uint32_t a[3]={t[0],t[1],t[0]^t[1]};int i=pick(3),j=pick(2);if(j>=i)j++;t[0]=a[i];t[1]=a[j];}
void reduce(V&v){bool again=true;while(again){again=false;for(size_t i=0;i<v.size();i++){if(zero(v[i])){remove(v,i);again=true;break;}canonical(v[i]);}if(again)continue;
 for(size_t i=0;i<v.size()&&!again;i++)for(size_t j=i+1;j<v.size();j++){
  if(v[i][0]==v[j][0]&&v[i][1]==v[j][1]){v[i][2]^=v[j][2];remove(v,j);again=true;break;}
  if(v[i][2]==v[j][2]){uint32_t a=common(v[i],v[j]);if(a){orient(v[i],a);orient(v[j],a);v[i][1]^=v[j][1];remove(v,j);again=true;break;}}
 }
}}
struct Move{int i,j;uint32_t a;};
std::vector<Move> flips(const V&v){std::vector<Move> out;for(int i=0;i<(int)v.size();i++)for(int j=i+1;j<(int)v.size();j++){uint32_t a=common(v[i],v[j]);if(a)out.push_back({i,j,a});if(v[i][2]==v[j][2])out.push_back({i,j,0});}return out;}
void flip(V&v,Move m){int i=m.i,j=m.j;if(pick(2))std::swap(i,j);if(m.a){orient(v[i],m.a);orient(v[j],m.a);if(pick(2)){v[i][1]^=v[j][1];v[j][2]^=v[i][2];}else{v[i][2]^=v[j][2];v[j][1]^=v[i][1];}}
 else {random_basis(v[i]);random_basis(v[j]);v[i][0]^=v[j][0];v[j][1]^=v[i][1];}reduce(v);}
void plus(V&v){int i=pick(v.size()),j=pick(v.size()-1);if(j>=i)j++;random_basis(v[i]);random_basis(v[j]);int x=pick(3),y=(x+1)%3,z=(x+2)%3;if(pick(2))std::swap(y,z);T a=v[i],b=v[j],p=a,q=b,r=b;p[x]^=b[x];q[y]^=a[y];r[y]=a[y];r[z]^=a[z];v[i]=p;v[j]=q;v.push_back(r);reduce(v);}
std::array<uint16_t,276> tensor(const V&v){std::array<uint16_t,276> out{};for(auto t:v){int k=0;for(int i=0;i<24;i++)for(int j=i+1;j<24;j++,k++)if(((t[0]>>i&1)&(t[1]>>j&1))^((t[0]>>j&1)&(t[1]>>i&1)))out[k]^=t[2];}return out;}
void save(const V&v,const char*path){std::ofstream f(path);for(auto t:v)f<<t[0]<<' '<<t[1]<<' '<<t[2]<<'\n';}
int main(int argc,char**argv){int seconds=argc>1?std::stoi(argv[1]):120;bool direct=argc>2;const char* output=argc>3?argv[3]:"parameters/quadratic_multiplier_best.txt";std::ifstream in(direct?argv[2]:"parameters/mul42_tensor.txt");V initial;uint32_t a,b,c;while(in>>a>>b>>c)initial.push_back({a,direct?b:b<<12,c});auto want=tensor(initial);V best=initial,v=initial;reduce(v);uint64_t steps=0,nplus=0,nflip=0,restarts=0;int age=0;auto start=std::chrono::steady_clock::now();
 while(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<seconds){steps++;auto f=flips(v);if((f.empty()||pick(100)<5)&&v.size()<best.size()+8){plus(v);nplus++;}else if(!f.empty()){flip(v,f[pick(f.size())]);nflip++;}else age=20000;
  if(v.size()<best.size()){if(tensor(v)!=want)return 2;best=v;save(best,output);std::cout<<"best "<<best.size()<<" steps "<<steps<<std::endl;age=0;}
  if(++age>=20000){v=best;age=0;restarts++;}
 }
 if(tensor(v)!=want||tensor(best)!=want)return 3;save(best,output);std::cout<<"{\"status\":\"EXACT_QUADRATIC_SEARCH_COMPLETE\",\"initial_rank\":"<<initial.size()<<",\"best_rank\":"<<best.size()<<",\"steps\":"<<steps<<",\"flips\":"<<nflip<<",\"plus\":"<<nplus<<",\"restarts\":"<<restarts<<",\"seconds\":"<<seconds<<",\"seed\":2026102616,\"scope\":\"full 276-coefficient quadratic tensor equality; affine corrections still required before native use\"}"<<std::endl;
}
