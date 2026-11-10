// Polar quadratic covers: wedge-product flips plus rank-increasing moves.
// Adapted from round-four quadratic_flip.cpp and tensor_flip.cpp.
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
using namespace std; using T=array<unsigned,3>;using U=unsigned __int128;
U polar(unsigned a,unsigned b){U p=0;int k=0;for(int i=0;i<12;i++)for(int j=i+1;j<12;j++,k++)if((((a>>i)&1)&((b>>j)&1))^(((a>>j)&1)&((b>>i)&1)))p|=U(1)<<k;return p;}
void tidy(vector<T>&v){for(size_t i=0;i<v.size();i++)for(size_t j=i+1;j<v.size();)if(polar(v[i][0],v[i][1])==polar(v[j][0],v[j][1])){v[i][2]^=v[j][2];v.erase(v.begin()+j);}else j++;v.erase(remove_if(v.begin(),v.end(),[](T t){return !t[2]||!polar(t[0],t[1]);}),v.end());}
void save(const vector<T>&v,const char*p){ofstream f(p);f<<v.size()<<'\n';for(auto t:v)f<<t[0]<<' '<<t[1]<<' '<<t[2]<<'\n';}
int main(int argc,char**argv){if(argc!=6)return 2;ifstream f(argv[1]);int n;f>>n;vector<T>best(n);for(auto&t:best)f>>t[0]>>t[1]>>t[2];tidy(best);mt19937 rng(stoul(argv[4]));long long count=stoll(argv[3]);int cap=stoi(argv[5]);save(best,argv[2]);cout<<"initial "<<best.size()<<endl;auto v=best;
for(long long it=0;it<count;it++){if(it%100000==0)v=best;if(v.size()<2)break;int i=rng()%v.size(),j=rng()%v.size();if(i==j)continue;
 // Every nonzero vector in span(a,b) can be either wedge operand.
 for(int at:{i,j}){auto &t=v[at];unsigned a=t[0],b=t[1];switch(rng()%6){case 0:break;case 1:t[0]=b;t[1]=a;break;case 2:t[1]=a^b;break;case 3:t[0]=a^b;break;case 4:t[0]=b;t[1]=a^b;break;case 5:t[0]=a^b;t[1]=a;break;}}
 bool changed=false;
 if(v[i][0]==v[j][0]){v[i][1]^=v[j][1];v[j][2]^=v[i][2];changed=true;}
 else if(v[i][2]==v[j][2]){v[i][0]^=v[j][0];v[j][1]^=v[i][1];changed=true;}
 else if(v.size()<best.size()+cap&&rng()%100==0){auto s=v[i],t=v[j];v[i]={s[0]^t[0],s[1],s[2]};v[j]={t[0],s[1]^t[1],t[2]};v.push_back({t[0],s[1],s[2]^t[2]});changed=true;}
 if(!changed)continue;tidy(v);if(v.size()<best.size()){best=v;save(best,argv[2]);cout<<"{\"iteration\":"<<it<<",\"ANDs\":"<<best.size()<<"}"<<endl;}}
save(best,argv[2]);}
