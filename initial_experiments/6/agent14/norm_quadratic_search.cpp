// Quadratic-cover search modulo linear maps; canonical AND planes.
// Local generalization of round-four quadratic_tensor_flip.cpp.
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
using namespace std;using T=array<unsigned,3>;
void canon(T&t){unsigned a=t[0],b=t[1],c=a^b;t[0]=min({a,b,c});t[1]=a^b^c; if(t[0]==a)t[1]=min(b,c);else if(t[0]==b)t[1]=min(a,c);else t[1]=min(a,b);}
void tidy(vector<T>&v){for(auto&t:v)canon(t);for(size_t i=0;i<v.size();i++)for(size_t j=i+1;j<v.size();)if(v[i][0]==v[j][0]&&v[i][1]==v[j][1]){v[i][2]^=v[j][2];v.erase(v.begin()+j);}else j++;v.erase(remove_if(v.begin(),v.end(),[](T t){return !t[2]||!t[0]||!t[1]||t[0]==t[1];}),v.end());}
int main(int argc,char**argv){if(argc!=7)return 2;ifstream in(argv[1]);int n;in>>n;vector<T>best(n);for(auto&t:best)in>>t[0]>>t[1]>>t[2];long long iterations=stoll(argv[3]);int chains=stoi(argv[4]),cap=stoi(argv[6]);mt19937 rng(stoul(argv[5]));for(int chain=0;chain<chains;chain++){auto v=best;long long stale=0;for(long long it=0;it<iterations;it++){if(++stale>500000){v=best;stale=0;}if(v.empty())return 3;int i=rng()%v.size(),j=rng()%v.size();if(i==j)continue; // Random equivalent bases in each quadratic plane.
for(int k:{i,j}){if(rng()&1)swap(v[k][0],v[k][1]);if(rng()&1)v[k][1]^=v[k][0];}
int axis=rng()%3,b=(axis+1)%3,c=(axis+2)%3;if(v[i][axis]==v[j][axis]){v[i][b]^=v[j][b];v[j][c]^=v[i][c];}else if(v.size()<best.size()+cap&&rng()%80==0){auto s=v[i],t=v[j];T A=s,B=t,C=t;A[axis]^=t[axis];B[b]^=s[b];B[c]=s[c];C[c]^=s[c];v[i]=A;v[j]=B;v.push_back(C);}else continue;tidy(v);if(v.size()<best.size()){best=v;stale=0;ofstream out(argv[2]);out<<best.size()<<'\n';for(auto t:best)out<<t[0]<<' '<<t[1]<<' '<<t[2]<<'\n';cout<<"{\"chain\":"<<chain<<",\"iteration\":"<<it<<",\"ANDs\":"<<best.size()<<"}"<<endl;}}}ofstream out(argv[2]);out<<best.size()<<'\n';for(auto t:best)out<<t[0]<<' '<<t[1]<<' '<<t[2]<<'\n';}
