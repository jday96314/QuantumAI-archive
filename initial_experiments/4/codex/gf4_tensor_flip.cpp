// Field-linear tensor flip/plus search over GF4, with projective alignment.
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
using namespace std;using T=array<unsigned,3>;
unsigned scale(unsigned v,unsigned s){if(!s)return 0;if(s==1)return v;unsigned twice=((v&0xaaa)>>1)^((v&0x555)<<1)^(v&0xaaa);return s==2?twice:twice^v;}
unsigned first(unsigned v){while(!(v&3))v>>=2;return v&3;}
unsigned inv(unsigned v){return v==1?1:v^1;}
void canonical(T&t){if(!t[0]||!t[1]||!t[2]){t={0,0,0};return;}unsigned a=first(t[0]),b=first(t[1]);t[0]=scale(t[0],inv(a));t[1]=scale(t[1],inv(b));t[2]=scale(scale(t[2],a),b);}
void tidy(vector<T>&v){for(auto&t:v)canonical(t);for(size_t i=0;i<v.size();i++)for(size_t j=i+1;j<v.size();)if(v[i][0]==v[j][0]&&v[i][1]==v[j][1]){v[i][2]^=v[j][2];v.erase(v.begin()+j);}else j++;v.erase(remove_if(v.begin(),v.end(),[](T t){return !t[0]||!t[1]||!t[2];}),v.end());}
int main(int argc,char**argv){if(argc!=7)return 2;ifstream in(argv[1]);int n;in>>n;vector<T>best(n);for(auto&t:best)in>>t[0]>>t[1]>>t[2];long long iterations=stoll(argv[3]);int chains=stoi(argv[4]),cap=stoi(argv[6]);mt19937 rng(stoul(argv[5]));for(int chain=0;chain<chains;chain++){auto v=best;for(long long it=0;it<iterations;it++){if(it%10000==0)v=best;if(v.size()<2)return 3;int i=rng()%v.size(),j=rng()%v.size();if(i==j)continue;if(rng()%100==0)for(auto&t:v)swap(t[0],t[1]);int axis=rng()%3,b=(axis+1)%3,c=(axis+2)%3;unsigned left=v[i][axis],right=v[j][axis];if(!left||!right)continue;unsigned ratio=first(left);ratio=scale(first(right),inv(ratio));ratio=inv(ratio);bool aligned=left==scale(right,ratio);if(aligned){v[j][axis]=left;v[j][b]=scale(v[j][b],inv(ratio));unsigned factor=1+rng()%3;v[i][b]^=scale(v[j][b],factor);v[j][c]^=scale(v[i][c],factor);}else if(v.size()<best.size()+cap&&rng()%50==0){auto s=v[i],t=v[j];T A=s,B=t,C=t;A[axis]^=t[axis];B[b]^=s[b];B[c]=s[c];C[c]^=s[c];v[i]=A;v[j]=B;v.push_back(C);}else continue;tidy(v);if(v.size()<best.size()){best=v;ofstream out(argv[2]);out<<best.size()<<'\n';for(auto t:best)out<<t[0]<<' '<<t[1]<<' '<<t[2]<<'\n';cout<<"{\"chain\":"<<chain<<",\"iteration\":"<<it<<",\"rank\":"<<best.size()<<"}"<<endl;}}}ofstream out(argv[2]);out<<best.size()<<'\n';for(auto t:best)out<<t[0]<<' '<<t[1]<<' '<<t[2]<<'\n';}
