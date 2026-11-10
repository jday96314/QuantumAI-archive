// Quadratic AND-cover flips. Linear terms are repaired separately after search.
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
using namespace std;
struct T{int a,b,o;};
int polar(int a,int b){int p=0,k=0;for(int i=0;i<8;i++)for(int j=i+1;j<8;j++,k++)p|=((((a>>i)&1)&((b>>j)&1))^(((a>>j)&1)&((b>>i)&1)))<<k;return p;}
void tidy(vector<T>&v){for(size_t i=0;i<v.size();i++)for(size_t j=i+1;j<v.size();)if(polar(v[i].a,v[i].b)==polar(v[j].a,v[j].b)){v[i].o^=v[j].o;v.erase(v.begin()+j);}else j++;v.erase(remove_if(v.begin(),v.end(),[](T t){return !t.o||!polar(t.a,t.b);}),v.end());}
int main(int argc,char**argv){if(argc!=5)return 2;ifstream f(argv[1]);int n;f>>n;vector<T> initial(n);for(auto&t:initial)f>>t.a>>t.b>>t.o;long long iterations=stoll(argv[3]);int chains=stoi(argv[4]);vector<T> best=initial;mt19937 rng(791361);
for(int chain=0;chain<chains;chain++){vector<T>v=best;for(long long it=0;it<iterations;it++){if(v.size()<2)break;if((rng()%100)==0 && v.size()<best.size()+2){int q=rng()%v.size(),u=1+(rng()%255);auto t=v[q];v[q].b^=u;v.push_back({t.a,u,t.o});tidy(v);}int i=rng()%v.size(),j=rng()%v.size();if(i==j)continue;array<int,3>ai={v[i].a,v[i].b,v[i].a^v[i].b},aj={v[j].a,v[j].b,v[j].a^v[j].b};int g=0;for(int x:ai)for(int y:aj)if(x==y)g=x;if(!g)continue;int u=0,w=0;for(int x:ai)if(x!=g){u=x;break;}for(int x:aj)if(x!=g){w=x;break;}v[i].a=g;v[i].b=u^w;v[j].a=g;v[j].b=w;v[j].o^=v[i].o;tidy(v);if(v.size()<best.size()){best=v;ofstream o(argv[2]);o<<best.size()<<'\n';for(auto t:best)o<<t.a<<' '<<t.b<<' '<<t.o<<'\n';cout<<"{\"chain\":"<<chain<<",\"iteration\":"<<it<<",\"ANDs\":"<<best.size()<<"}"<<endl;}}}
ofstream o(argv[2]);o<<best.size()<<'\n';for(auto t:best)o<<t.a<<' '<<t.b<<' '<<t.o<<'\n';}
