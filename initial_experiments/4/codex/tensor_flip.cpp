// Exact three-factor bilinear tensor flips and two-term plus transitions.
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
using namespace std;using T=array<unsigned,3>;
void tidy(vector<T>&v){bool changed=true;while(changed){changed=false;for(size_t i=0;i<v.size();i++)for(size_t j=i+1;j<v.size();j++){int same=(v[i][0]==v[j][0])+(v[i][1]==v[j][1])+(v[i][2]==v[j][2]);if(same<2)continue;for(int k=0;k<3;k++)if(v[i][k]!=v[j][k]||same==3){v[i][k]^=v[j][k];break;}v.erase(v.begin()+j);changed=true;break;}v.erase(remove_if(v.begin(),v.end(),[](T t){return !t[0]||!t[1]||!t[2];}),v.end());}}
int main(int argc,char**argv){if(argc!=7)return 2;ifstream in(argv[1]);int n;in>>n;vector<T>best(n);for(auto&t:best)in>>t[0]>>t[1]>>t[2];long long iterations=stoll(argv[3]);int chains=stoi(argv[4]),cap=stoi(argv[6]);mt19937 rng(stoul(argv[5]));for(int chain=0;chain<chains;chain++){auto v=best;for(long long it=0;it<iterations;it++){int i=rng()%v.size(),j=rng()%v.size();if(i==j)continue;int axis=rng()%3;int b=(axis+1)%3,c=(axis+2)%3;if(v[i][axis]==v[j][axis]){v[i][b]^=v[j][b];v[j][c]^=v[i][c];}else if(v.size()<best.size()+cap && (rng()%200)==0){auto s=v[i],t=v[j];T A=s,B=t,C=t;A[axis]^=t[axis];B[b]^=s[b];B[c]=s[c];C[c]^=s[c];v[i]=A;v[j]=B;v.push_back(C);}else continue;tidy(v);if(v.size()<best.size()){best=v;ofstream out(argv[2]);out<<best.size()<<'\n';for(auto t:best)out<<t[0]<<' '<<t[1]<<' '<<t[2]<<'\n';cout<<"{\"chain\":"<<chain<<",\"iteration\":"<<it<<",\"rank\":"<<best.size()<<"}"<<endl;}}}ofstream out(argv[2]);out<<best.size()<<'\n';for(auto t:best)out<<t[0]<<' '<<t[1]<<' '<<t[2]<<'\n';}
