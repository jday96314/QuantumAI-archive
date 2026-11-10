#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
int main(int argc,char**argv){
 if(argc!=9&&argc!=10)throw runtime_error("tape hook_index bitbase pre sweep post output maxout [phases]");
 ifstream f(argv[1],ios::binary);char magic[8];f.read(magic,8);
 uint32_t hlen=0;if((uint8_t)magic[6]==1){uint16_t x;f.read((char*)&x,2);hlen=x;}
 else f.read((char*)&hlen,4);
 f.seekg(hlen,ios::cur);auto start=f.tellg();f.seekg(0,ios::end);uint64_t n=(f.tellg()-start)/36;
 vector<array<uint32_t,9>> rows(n);f.seekg(start);f.read((char*)rows.data(),n*36);
 uint64_t hook=stoull(argv[2]);uint32_t base=stoul(argv[3]);
 vector<uint64_t> last(stoul(argv[8]),0);vector<int> wires(last.size(),-1),pool;
 for(uint64_t i=0;i<n;i++)for(int j:{4,5})if(!(rows[i][0]>>(7+j)&1))last.at(rows[i][j])=i;
 ofstream out(argv[7]);auto cat=[&](char*p){ifstream in(p);if(in.peek()!=ifstream::traits_type::eof())out<<in.rdbuf();};cat(argv[4]);
 auto bit=[&](uint32_t id){int &w=wires.at(id);if(w<0){if(pool.empty())w=base++;else{w=pool.back();pool.pop_back();}}return w;};
 string names[]={"X","Z","NEG","CX","CZ","SWAP","CCX","CCZ","R","HMR","CX","BIT_INVERT","BIT_STORE0","BIT_STORE1"};
 vector<pair<uint64_t,string>> marks;size_t mark=0;
 if(argc==10){ifstream mm(argv[9]);uint64_t index;string name;while(mm>>index>>name)marks.push_back({index,name});}
 for(uint64_t i=0;i<n;i++){
  while(mark<marks.size()&&marks[mark].first==i){out<<"INFO_TAG \"DBG:"<<marks[mark].second<<"\"\n";mark++;}
  if(i==hook)cat(argv[5]);auto&r=rows[i];int k=r[0]&255;
  bool ot=!(r[0]>>11&1),oc=!(r[0]>>12&1),cc=!(r[0]>>15&1);
  if(cc)throw runtime_error("unexpected classical bit condition");
  if(oc&&k!=10)out<<"PUSH_CONDITION if b"<<bit(r[5])<<'\n';
  out<<names[k];
  if(k==6||k==7)out<<" q"<<r[1];
  if(k>=3&&k<=7)out<<" q"<<r[2];
  if(k!=2&&k<10)out<<" q"<<r[3];
  if(k==9)out<<" b"<<bit(r[4]);
  if(k==10)out<<" b"<<bit(r[5])<<" b"<<bit(r[4]);
  if(k>=11)throw runtime_error("unexpected bit instruction");
  out<<'\n';if(oc&&k!=10)out<<"POP_CONDITION\n";
  for(int j:{4,5})if((j==4?ot:oc)&&last[r[j]]==i&&wires[r[j]]>=0){
   out<<"BIT_STORE0 b"<<wires[r[j]]<<'\n';pool.push_back(wires[r[j]]);wires[r[j]]=-1;
  }
 }
 cat(argv[6]);cerr<<"ops "<<n<<" bits "<<base<<'\n';
}
