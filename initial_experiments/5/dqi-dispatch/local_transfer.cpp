// Autonomous compact-band Euclid transfers. Field tables/truth are read-only
// inputs from dqi-repair. No exact-EEA dispatch hint enters the decoder.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
using P=vector<int>;
int lg[4096],ex[8190];
int mul(int a,int b){return a&&b?ex[lg[a]+lg[b]]:0;}
int inv(int a){return a?ex[4095-lg[a]]:0;}
void trim(P&a){while(!a.empty()&&!a.back())a.pop_back();}
P add(P a,const P&b){a.resize(max(a.size(),b.size()));for(size_t i=0;i<b.size();i++)a[i]^=b[i];trim(a);return a;}
P prod(const P&a,const P&b){if(a.empty()||b.empty())return {};P c(a.size()+b.size()-1);for(size_t i=0;i<a.size();i++)for(size_t j=0;j<b.size();j++)c[i+j]^=mul(a[i],b[j]);trim(c);return c;}
P scale(P a,int v){for(int&x:a)x=mul(x,v);trim(a);return a;}
P rev(P a){reverse(a.begin(),a.end());return a;}
P sp(const P&a,const P&b,int n){P c(n);for(int i=0;i<(int)a.size()&&i<n;i++)for(int j=0;j<(int)b.size()&&i+j<n;j++)c[i+j]^=mul(a[i],b[j]);return c;}
P sd(const P&a,const P&b,int n){if(b.empty()||!b[0])throw runtime_error("zero series divisor");P c(n);int bi=inv(b[0]);for(int i=0;i<n;i++){int z=i<(int)a.size()?a[i]:0;for(int j=1;j<=i&&j<(int)b.size();j++)z^=mul(b[j],c[i-j]);c[i]=mul(z,bi);}return c;}
P monic(P a){trim(a);return a.empty()?a:scale(a,inv(a.back()));}
int eval(const P&a,int x){int z=0;for(int i=(int)a.size()-1;i>=0;i--)z=mul(z,x)^a[i];return z;}
P bm(const P&S){P c(71),b(71);c[0]=b[0]=1;int L=0,m=1,p=1;for(int n=0;n<70;n++){int d=S[n];for(int i=1;i<=L;i++)d^=mul(c[i],S[n-i]);if(!d){m++;continue;}P t=c;int v=mul(d,inv(p));for(int j=0;j+m<=70;j++)c[j+m]^=mul(v,b[j]);if(2*L<=n){L=n+1-L;b=t;p=d;m=1;}else m++;}c.resize(L+1);return monic(c);}
struct State{P bp{},b{1},c;int pi=1;};
map<int,State> exact_states(const P&S){
 // Conventional unnormalized polynomial EEA, independent of moment extraction.
 // These states are validation only and never determine a decoder branch.
 map<int,State> out;out[0]=State{{},{1},S,1};P p(71),r=S,fp{},f{1};p[70]=1;trim(r);
 while(!r.empty()){
  P rr=p;P q(max(0,(int)p.size()-(int)r.size()+1));int ri=inv(r.back());
  while(rr.size()>=r.size()){
   int g=rr.size()-r.size(),v=mul(rr.back(),ri);q[g]=v;
   for(int i=0;i<(int)r.size();i++)rr[g+i]^=mul(v,r[i]);trim(rr);
  }
  P nf=add(fp,prod(q,f));fp=f;f=nf;p=r;r=rr;int d=f.size()-1;if(d>35)break;
  P band(70-2*d);for(int j=d;j<70-d;j++)if(j<(int)r.size())band[j-d]=mul(r[j],inv(f.back()));
  out[d]=State{monic(fp),monic(f),band,mul(fp.back(),f.back())};
 }
 return out;
}
struct Transfer{P a{1},b{},c{},d{1};int beta=1,used=0;string pattern;};
struct Stats{long blocks=0,carry=0,carry_after_gap=0,completion=0,immediate_gap_stop=0,serial_gap=0,final_two=0,max_head=0,audited_states=0;map<string,long>patterns;array<long,36>gap_stages{};};
// F is reverse(R)/reverse(P), aligned so R/P=z*F. C retains R[d:N-d].
P moments(const State&s,int n){int d=s.b.size()-1,h=(int)s.bp.size()-1;P t=rev(s.c),den(d-h+1);den[0]=inv(s.pi);P z=sp(t,rev(s.bp),max(0,n-(d-h+1)));den.insert(den.end(),z.begin(),z.end());return sd(sp(t,rev(s.b),n),den,n);}
Transfer extract(P f,int budget,int stage,Stats&st,int maxgap){
 Transfer t;int left=budget;
 while(left){int zeros=0;while(zeros<(int)f.size()&&!f[zeros])zeros++;
  if(zeros==(int)f.size())break;
  int g=zeros+1;
  if(g>maxgap)throw runtime_error("unsupported quotient degree");
  if(g>left){st.carry++;break;}
  if((int)f.size()<2*g)throw runtime_error("insufficient local moments");
  int beta=f[zeros];P sh(f.begin()+zeros,f.end()),u=sd(P{beta},sh,sh.size());
  P q(u.begin(),u.begin()+g+1);q=rev(q);
  P na=t.c,nb=t.d;
  t.c=add(scale(t.a,beta),prod(q,t.c));t.d=add(scale(t.b,beta),prod(q,t.d));t.a=na;t.b=nb;
  t.beta=mul(t.beta,beta);t.used+=g;left-=g;
  t.pattern+=(t.pattern.empty()?"":",")+to_string(g);
  if(g>1){st.gap_stages[stage+t.used-g]++;if(stage>=22)st.serial_gap++;if(stage+t.used==35&&g==2)st.final_two++;}
  f=P(u.begin()+g+1,u.end());
 }
 if(!t.used)throw runtime_error("no local progress");
 st.patterns[t.pattern]++;st.blocks++;
 return t;
}
State apply_transfer(const State&s,const Transfer&t){
 State r;r.bp=add(prod(t.a,s.bp),prod(t.b,s.b));r.b=add(prod(t.c,s.bp),prod(t.d,s.b));r.pi=mul(s.pi,inv(t.beta));
 P num=sp(rev(r.b),rev(s.c),s.c.size());int n=(int)s.c.size()-2*t.used;
 if(n<0)throw runtime_error("past final degree");
 P tail(num.begin()+2*t.used,num.end());r.c=rev(sd(tail,rev(s.b),n));
 return r;
}
State decode(P S,Stats&st,int maxgap=3,const map<int,State>*audit=nullptr){
 State s;s.c=S;int previous_gap=0,carry=0;
 while((int)s.b.size()-1<35){int degree=s.b.size()-1;
  if(none_of(s.c.begin(),s.c.end(),[](int x){return x;})){st.completion++;if(previous_gap>1)st.immediate_gap_stop++;break;}
  int nominal=degree<22?3:2;int budget=min(nominal,35-degree);
  // Each quantum block must provision this maximum, independent of rarity.
  int n=min((int)s.c.size(),2*budget);st.max_head=max(st.max_head,(long)n);
  P f=moments(s,n);Transfer t=extract(f,budget,degree,st,maxgap);
  int last=t.pattern.back()-'0';
  State r=apply_transfer(s,t);
  if((int)r.b.size()-1!=degree+t.used)throw runtime_error("incorrect transfer degree");
  if(audit){
   const State&want=audit->at(degree+t.used);
   if(r.bp!=want.bp||r.b!=want.b||r.c!=want.c||r.pi!=want.pi)throw logic_error("compact/full EEA state mismatch");
   st.audited_states++;
  }
  P det=add(prod(t.a,t.d),prod(t.b,t.c));if(det!=P{t.beta})throw runtime_error("transfer determinant");
  // A final degree-two quotient leaves one fewer stored Bprev coefficient.
  // Its physical implementation requires the promised word realignment.
  if(carry&&last>1)st.carry_after_gap++;
  carry=t.used<budget;previous_gap=last;s=move(r);
 }
 return s;
}
vector<pair<int,int>> decoded(const State&s){vector<pair<int,int>>out;for(int j=0;j<4095;j++){int x=ex[(4095-j)%4095];if(eval(s.b,x))continue;int der=0;for(int i=1;i<(int)s.b.size();i+=2)der^=mul(s.b[i],ex[(lg[x]*(i-1))%4095]);int v=mul(ex[(lg[x]*70)%4095],inv(mul(mul(s.pi,eval(s.bp,x)),der)));if(v)out.push_back({j,v});}return out;}
bool correct(const State&s,const P&pos,const P&val,const P&locator){if(s.b!=locator)return false;for(int j=0;j<(int)pos.size();j++){int x=ex[(4095-pos[j])%4095],der=0;for(int i=1;i<(int)s.b.size();i+=2)der^=mul(s.b[i],ex[(lg[x]*(i-1))%4095]);int v=mul(ex[(lg[x]*70)%4095],inv(mul(mul(s.pi,eval(s.bp,x)),der)));if(v!=val[j])return false;}return true;}
void pairs(ostream&o,const vector<pair<int,int>>&v){o<<'[';for(int i=0;i<(int)v.size();i++)o<<(i?",":"")<<'['<<v[i].first<<','<<v[i].second<<']';o<<']';}
#include "fixed_schedule.hpp"
int main(int argc,char**argv){
 int N=argc>1?stoi(argv[1]):102400;uint64_t seed=argc>2?stoull(argv[2]):2026102405;string mode=argc>3?argv[3]:"host";int maxgap=argc>4?stoi(argv[4]):3;bool fixed=argc>5&&string(argv[5])=="fixed";
 ifstream al("../dqi-repair/parameters/antilog.bin",ios::binary),ll("../dqi-repair/parameters/log.bin",ios::binary);al.read((char*)ex,4095*4);ll.read((char*)lg,4096*4);if(!al||!ll)throw runtime_error("field tables");copy(ex,ex+4095,ex+4095);
 vector<double>pr(36);ifstream prob("../dqi-repair/parameters/probs.txt");for(auto&x:pr)prob>>x;
 mt19937_64 rng(seed),mrng(seed^0x473df126adULL);discrete_distribution<int>wd(pr.begin(),pr.end());uniform_int_distribution<int>pd(0,4094),vd(1,4095),md(0,4095);
 ifstream host;if(mode=="host")host.open("../dqi-repair/parameters/host_cases.bin",ios::binary);
 string stem=(fixed?"fixed_":"")+mode+"_"+to_string(N)+"_"+to_string(seed)+"_g"+to_string(maxgap);ofstream failures("results/"+stem+"_failures.jsonl");
 long bad=0,phase=0,bmchecks=0,exceptions=0,peak_words=0;Stats st;auto start=chrono::steady_clock::now();
 for(int shot=0;shot<N;shot++){P S(70),pos,val;int w;
  if(mode=="host"){uint16_t rec[141];host.read((char*)rec,sizeof rec);if(!host)throw runtime_error("host input");copy(rec,rec+70,S.begin());w=rec[70];pos=P(rec+71,rec+71+w);val=P(rec+106,rec+106+w);}
  else{w=mode=="controls"?shot%36:wd(rng);while((int)pos.size()<w){int j=pd(rng);if(find(pos.begin(),pos.end(),j)==pos.end()){pos.push_back(j);val.push_back(vd(rng));}}for(int j=0;j<w;j++)for(int i=0;i<70;i++)S[i]^=mul(val[j],ex[((i+1)*pos[j])%4095]);}
  P locator=bm(S);bmchecks++;State s;bool failed=false;
  string reason;
  try{auto audit=(fixed||shot<1024||mode=="controls")?exact_states(S):map<int,State>{};s=fixed?decode_fixed(S,st,maxgap,&audit):decode(S,st,maxgap,audit.empty()?nullptr:&audit);peak_words=max(peak_words,(long)(s.bp.size()+s.b.size()+s.c.size()-1));failed=!correct(s,pos,val,locator);}
  catch(const runtime_error&e){exceptions++;failed=true;s.b={1};s.bp={};s.pi=1;reason=e.what();}
  if(failed){bad++;auto got=decoded(s);int par=0;for(int j=0;j<4095;j++){int a=0,b=0;for(int i=0;i<w;i++)if(pos[i]==j)a=val[i];for(auto [p,v]:got)if(p==j)b=v;par^=__builtin_parity((unsigned)((a^b)&md(mrng)));}phase+=par;failures<<"{\"shot\":"<<shot<<",\"reason\":\""<<reason<<"\",\"syndromes\":[";for(int i=0;i<70;i++)failures<<(i?",":"")<<S[i];failures<<"],\"decoded\":";pairs(failures,got);failures<<"}\n";}
  if((shot+1)%100000==0)cerr<<shot+1<<" cases, "<<bad<<" wrong, "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<" seconds\n";
 }
 ostream&o=cout;if(fixed){o<<"{\"fixed_schedule\":true,\"fixed_statistics\":";fixed_report(o);o<<",";}o<<(fixed?"":"{")<<"\"mode\":\""<<mode<<"\",\"cases\":"<<N<<",\"seed\":"<<seed<<",\"max_supported_quotient\":"<<maxgap<<",\"wrong_decodes\":"<<bad<<",\"independent_mask_phase_errors\":"<<phase<<",\"unsupported\":"<<exceptions<<",\"bm_comparisons\":"<<bmchecks<<",\"full_eea_state_audits\":"<<st.audited_states<<",\"max_final_field_words\":"<<peak_words<<",\"blocks\":"<<st.blocks<<",\"boundary_carries\":"<<st.carry<<",\"carry_followed_by_gap\":"<<st.carry_after_gap<<",\"completion\":"<<st.completion<<",\"immediate_gap_completion\":"<<st.immediate_gap_stop<<",\"serial_gap\":"<<st.serial_gap<<",\"final_two_gap\":"<<st.final_two<<",\"max_head_words\":"<<st.max_head<<",\"patterns\":{";
 bool first=true;for(auto [pat,n]:st.patterns){o<<(first?"":",")<<'"'<<pat<<"\":"<<n;first=false;}o<<"},\"gap_stage_hist\":[";for(int i=0;i<36;i++)o<<(i?",":"")<<st.gap_stages[i];o<<"],\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
}
