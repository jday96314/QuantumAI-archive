// Reproduce the error RNG stream without rerunning the full million-case census.
#define main census_original_main
#include "carried_census.cpp"
#undef main
#include <set>
int main(){
 ifstream al("parameters/antilog.bin",ios::binary),ll("parameters/log.bin",ios::binary);al.read((char*)ex,4095*4);ll.read((char*)lg,4096*4);copy(ex,ex+4095,ex+4095);
 ifstream pf("parameters/probs.txt");vector<double>prob(36);for(auto&v:prob)pf>>v;
 ofstream out("parameters/residual_inputs.jsonl");
 for(int run=0;run<2;run++){
  uint64_t seed=2026102405+run;set<int>wanted=run?set<int>{519490,862621}:set<int>{389052,879219,958156};
  mt19937_64 rng(seed);discrete_distribution<int>wd(prob.begin(),prob.end());uniform_int_distribution<int>posd(0,4094),vald(1,4095);
  for(int shot=0;shot<=*wanted.rbegin();shot++){
   int w=wd(rng);P pos,val;while((int)pos.size()<w){int j=posd(rng);if(find(pos.begin(),pos.end(),j)==pos.end()){pos.push_back(j);val.push_back(vald(rng));}}
   if(!wanted.count(shot))continue;P S(70);for(int j=0;j<w;j++)for(int i=0;i<70;i++)S[i]^=mul(val[j],ex[((i+1)*pos[j])%4095]);
   out<<"{\"seed\":"<<seed<<",\"shot\":"<<shot<<",\"weight\":"<<w<<",\"syndromes\":[";for(int i=0;i<70;i++)out<<(i?",":"")<<S[i];out<<"],\"positions\":[";for(int i=0;i<w;i++)out<<(i?",":"")<<pos[i];out<<"],\"values\":[";for(int i=0;i<w;i++)out<<(i?",":"")<<val[i];out<<"]}\n";
  }
 }
}
