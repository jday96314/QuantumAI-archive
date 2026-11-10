// Recycle private wires at unconditional resets and their final use.
// Final-use recycling requires the input circuit's clean-workspace contract.
#include <fstream>
#include <iostream>
#include <queue>
#include <string>
#include <vector>
#include <functional>
#include <cctype>
using namespace std;
vector<pair<size_t,size_t>> refs(const string&s){
 vector<pair<size_t,size_t>> r;
 if(s.rfind("INFO_TAG",0)==0||s.rfind("REGISTER",0)==0||s.rfind("APPEND_TO_REGISTER",0)==0)return r;
 for(size_t i=0;i<s.size();i++)if(s[i]=='q'&&(i==0||isspace(s[i-1]))&&i+1<s.size()&&isdigit(s[i+1])){
  size_t j=i+1;while(j<s.size()&&isdigit(s[j]))j++;r.push_back({i,j});i=j-1;
 }return r;
}
int main(int argc,char**argv){
 if(argc!=4)return 2;int pinned=stoi(argv[3]);vector<unsigned long long> last(20000,0);
 string line;unsigned long long n=0;ifstream a(argv[1]);
 while(getline(a,line)){for(auto [i,j]:refs(line))last.at(stoi(line.substr(i+1,j-i-1)))=n;n++;}a.close();
 vector<int> map(last.size(),-1);for(int i=0;i<pinned;i++)map[i]=i;
 priority_queue<int,vector<int>,greater<int>> pool;int next=pinned,depth=0,live=pinned,peak=live;
 ifstream in(argv[1]);ofstream out(argv[2]);n=0;
 while(getline(in,line)){
  auto rr=refs(line);vector<int> ids;string nl;size_t prev=0;
  for(auto [i,j]:rr){int id=stoi(line.substr(i+1,j-i-1));ids.push_back(id);
   if(map[id]<0){if(pool.empty())map[id]=next++;else{map[id]=pool.top();pool.pop();}live++;peak=max(peak,live);}
   nl+=line.substr(prev,i-prev)+"q"+to_string(map[id]);prev=j;
  }
  nl+=line.substr(prev);out<<nl<<'\n';
  bool reset=depth==0&&(line.rfind("HMR ",0)==0||line.rfind("R ",0)==0);
  for(int id:ids)if(id>=pinned&&map[id]>=0&&(reset||last[id]==n)){pool.push(map[id]);map[id]=-1;live--;}
  if(line.rfind("PUSH_CONDITION",0)==0)depth++;
  if(line.rfind("POP_CONDITION",0)==0)depth--;
  n++;
 }
 cerr<<"qubits "<<next<<" peak_live "<<peak<<" lines "<<n<<'\n';
}
