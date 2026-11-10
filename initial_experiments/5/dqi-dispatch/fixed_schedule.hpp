// Fixed boundaries and partial quotient carry. No matrices in this path.
struct FixedStats {
 long regions=0,checks=0,encodings=0,pending1=0,pending2=0,gap_units=0;
 long deferred_beta=0,deferred_cubic_high=0,completions=0,max_words=0,max_fresh=0,max_head=0;
 array<long,8> modes{};array<long,36> boundaries{};map<string,long> patterns;
};
FixedStats fixed_stats;
struct CarryFrame {
 P A,B,tail,packed;int pi=1,nominal=0,pending=0,completion=0,mode=0;
 array<int,2> history{};int degree=0,previous_degree=-1,kind=0;
 bool A_monic=true,initial=false;
};
bool same_state(const State&a,const State&b){return a.bp==b.bp&&a.b==b.b&&a.c==b.c&&a.pi==b.pi;}
int coef(const P&a,int i){return i>=0&&i<(int)a.size()?a[i]:0;}
P shift_poly(P a,int n){if(!a.empty())a.insert(a.begin(),n,0);return a;}
void record_gap(array<int,2>&history,int stage){
 int tag=stage+1;if(history[0]==tag||history[1]==tag)return;
 if(!history[0])history[0]=tag;else if(!history[1])history[1]=tag;
 else throw runtime_error("fixed history exceeds two missing degree slots");
 fixed_stats.gap_units++;
}
int carry_mode(int pending,int gap,int kind){
 if(!pending&&gap>=1&&gap<=3)return gap-1;
 if(pending==1&&kind==2&&gap>=1&&gap<=2)return 2+gap;
 if(pending==1&&kind==3&&gap==1)return 5;
 if(pending==2&&kind==3&&gap==1)return 6;
 throw runtime_error("fixed mode requires more than two gap slots");
}
int mode_from_history(int k,int completion,array<int,2>history){
 auto has=[&](int tag){return tag>0&&(history[0]==tag||history[1]==tag);};
 int pending=has(k)?(has(k-1)?2:1):0;
 int kind=pending==1?(has(k+1)?3:2):pending==2?3:0;
 int degree=k-completion-pending;
 int gp=1+(has(degree-1)?1+(int)has(degree-2):0);
 return carry_mode(pending,gp,kind);
}
CarryFrame encode_fixed(const State&s,int nominal,int completed,int kind,array<int,2>history){
 CarryFrame f;f.degree=s.b.size()-1;f.previous_degree=(int)s.bp.size()-1;
 f.nominal=nominal;f.completion=completed;f.pending=nominal-f.degree-completed;f.kind=kind;f.history=history;
 f.mode=carry_mode(f.pending,f.degree-f.previous_degree,kind);f.initial=!f.degree&&!completed;
 if(f.mode!=mode_from_history(nominal,completed,history))throw logic_error("mode cannot be derived from sparse history");
 f.A=s.bp;f.B=s.b;f.tail=s.c;f.pi=s.pi;
 if(completed){
  if(f.pending)throw logic_error("completion overlaps pending");
  f.A=shift_poly(f.A,completed);f.B=shift_poly(f.B,completed);f.tail.resize(70-2*(f.degree+completed));
 }else if(f.pending==1&&kind==2){
  P ff=moments(s,2);int beta=ff[1];if(ff[0]||!beta)throw logic_error("bad quadratic deferral");
  f.A=scale(f.A,beta);f.A_monic=false;f.pi=mul(s.pi,inv(beta));f.tail.resize(f.tail.size()-2);fixed_stats.deferred_beta++;
 }else if(f.pending==1&&kind==3){
  P ff=moments(s,2);if(ff[0]||ff[1])throw logic_error("bad cubic first deferral");f.tail.resize(f.tail.size()-2);
 }else if(f.pending==2&&kind==3){
  P ff=moments(s,4);int beta=ff[2];if(ff[0]||ff[1]||!beta)throw logic_error("bad cubic second deferral");
  int q2=mul(ff[3],inv(beta));f.A=add(scale(s.bp,beta),scale(shift_poly(s.b,2),q2));f.A.resize(f.degree+3);
  f.A_monic=false;f.pi=mul(s.pi,inv(beta));f.tail.resize(f.tail.size()-4);fixed_stats.deferred_cubic_high++;
 }
 f.packed=f.A;if(f.A_monic&&!f.packed.empty())f.packed.pop_back();
 f.packed.insert(f.packed.end(),f.B.begin(),f.B.end()-1);
 if(nominal){if((int)f.packed.size()>2*nominal-1)throw logic_error("cofactor chart exceeds fixed prefix");f.packed.resize(2*nominal-1);}
 f.packed.insert(f.packed.end(),f.tail.begin(),f.tail.end());
 if(nominal)f.packed.push_back(f.pi);
 if(f.packed.size()>70)throw logic_error("fixed chart exceeds 70 words");
 fixed_stats.max_words=max(fixed_stats.max_words,(long)f.packed.size());f.packed.resize(70);fixed_stats.modes[f.mode]++;fixed_stats.encodings++;return f;
}
CarryFrame unpack_fixed(P words,int nominal,int completion,int mode,array<int,2>history){
 CarryFrame f;f.nominal=nominal;f.completion=completion;f.mode=mode;f.history=history;
 if(mode<0||mode>6||words.size()!=70)throw logic_error("invalid packed frame");
 f.pending=mode>=3?(mode==6?2:1):0;f.kind=mode>=5?3:mode>=3?2:0;
 int gap=mode<3?mode+1:mode==4?2:1;f.degree=nominal-completion-f.pending;f.previous_degree=f.degree-gap;
 f.initial=!f.degree&&!completion;f.A_monic=!(mode==3||mode==4||mode==6);
 int ah=f.previous_degree;
 if(completion&&f.degree)ah+=completion;
 int alen=f.A_monic?max(0,ah):ah+1;if(mode==6)alen=f.degree+3;
 int blen=f.degree+completion,pos=0;
 f.A=P(words.begin(),words.begin()+alen);pos+=alen;
 if(f.A_monic&&ah>=0)f.A.push_back(1);
 f.B=P(words.begin()+pos,words.begin()+pos+blen);f.B.push_back(1);pos+=blen;
 if(nominal){
  for(int i=pos;i<2*nominal-1;i++)if(words[i])throw logic_error("dirty cofactor ghost");
  pos=2*nominal-1;f.pi=words[69];
 }
 f.tail=P(words.begin()+pos,words.begin()+pos+70-2*nominal);f.packed=move(words);return f;
}
State decode_frame(const CarryFrame&f){
 State s;s.bp=f.A;s.b=f.B;s.c=f.tail;s.pi=f.pi;
 if(f.completion){
  if(!s.bp.empty())s.bp.erase(s.bp.begin(),s.bp.begin()+f.completion);
  s.b.erase(s.b.begin(),s.b.begin()+f.completion);s.c.resize(70-2*f.degree);
 }else if(f.pending==1&&f.kind==2){
  int r=inv(f.pi),beta=f.previous_degree<0?r:coef(f.A,f.previous_degree);
  if(!beta)throw logic_error("zero carried quadratic beta");s.bp=scale(f.A,inv(beta));s.pi=mul(f.pi,beta);s.c.push_back(r);s.c.push_back(0);
 }else if(f.pending==1&&f.kind==3){s.c.push_back(0);s.c.push_back(0);}
 else if(f.pending==2&&f.kind==3){
  int q2=coef(f.A,f.degree+2),r=inv(f.pi),sub=mul(q2^coef(f.B,f.degree-1),r);
  P original=add(f.A,scale(shift_poly(f.B,2),q2));int beta=f.previous_degree<0?r:coef(original,f.previous_degree);
  if(!beta)throw logic_error("zero carried cubic beta");s.bp=scale(original,inv(beta));s.pi=mul(f.pi,beta);
  s.c.push_back(sub);s.c.push_back(r);s.c.push_back(0);s.c.push_back(0);
 }
 trim(s.bp);trim(s.b);return s;
}
State quotient_update(const State&s,int beta,const P&q){
 int g=q.size()-1;State r;r.bp=s.b;r.b=add(scale(s.bp,beta),prod(q,s.b));r.pi=mul(s.pi,inv(beta));
 int n=(int)s.c.size()-2*g;if(n<0)throw runtime_error("quotient crosses decoder end");
 P num=sp(rev(r.b),rev(s.c),s.c.size()),tail(num.begin()+2*g,num.end());r.c=rev(sd(tail,rev(s.b),n));return r;
}
State decode_fixed(P S,Stats&st,int maxgap=3,const map<int,State>*audit=nullptr){
 State s;s.c=S;array<int,2>history{};int completed=0,kind=0;CarryFrame physical=encode_fixed(s,0,0,0,history);
 vector<pair<int,int>> schedule{{0,1}};for(int k=1;k<22;k+=3)schedule.push_back({k,3});for(int k=22;k<35;k++)schedule.push_back({k,1});
 for(auto [start,count]:schedule){
  physical=unpack_fixed(physical.packed,start,physical.completion,physical.mode,physical.history);
  if(!same_state(decode_frame(physical),s))throw logic_error("carried input mismatch");
  int end=start+count,degree=s.b.size()-1,oldpending=physical.pending;fixed_stats.regions++;fixed_stats.boundaries[start]++;
  int head=min((int)s.c.size(),2*(count+oldpending)),fresh=max(0,head-2*oldpending);
  if(fresh>2*count)throw logic_error("extra fresh head words");fixed_stats.max_fresh=max(fixed_stats.max_fresh,(long)fresh);fixed_stats.max_head=max(fixed_stats.max_head,(long)head);
  P f=moments(s,head);int available=end-degree-completed;string pattern;
  while(available>0){
   if(none_of(s.c.begin(),s.c.end(),[](int x){return x;})){completed+=available;fixed_stats.completions+=available;available=0;kind=0;break;}
   int zeros=0;while(zeros<(int)f.size()&&!f[zeros])zeros++;
   int g=zeros==(int)f.size()?max(3,zeros+1):zeros+1;if(g>maxgap)throw runtime_error("fixed quotient exceeds degree three");kind=g;
   for(int j=0;j<g-1;j++)record_gap(history,degree+j);
   // Assert each degree-slot cut, including partial quotients entirely
   // inside a block. Native slots delay all quadratic MACs to slot two;
   // cubic slot one is empty, slot two prepares beta and q2*x^2*B.
   for(int p=1;p<g;p++){
    int cut=degree+completed+p;
    if(cut<=start||cut>=end)continue;
    CarryFrame part=encode_fixed(s,cut,completed,g,history);
    part=unpack_fixed(part.packed,cut,completed,part.mode,history);
    if(!same_state(decode_frame(part),s))throw logic_error("internal partial-slot mismatch");
    if(audit){auto it=audit->find(degree);if(it==audit->end()||!same_state(it->second,decode_frame(part)))throw logic_error("internal slot/EEA mismatch");}
    fixed_stats.checks++;
   }
   if(g>available){if(available==1)fixed_stats.pending1++;else if(available==2)fixed_stats.pending2++;break;}
   if((int)f.size()<2*g)throw logic_error("insufficient fixed head");
   int beta=f[g-1];P sh(f.begin()+g-1,f.end()),u=sd(P{beta},sh,sh.size());P q(u.begin(),u.begin()+g+1);q=rev(q);
   State next=quotient_update(s,beta,q);
   Stats dummy;Transfer t=extract(moments(s,2*g),g,degree,dummy,maxgap);if(!same_state(next,apply_transfer(s,t)))throw logic_error("fixed/variable mismatch");
   degree+=g;available-=g;pattern+=(pattern.empty()?"":",")+to_string(g);
   if(audit){auto it=audit->find(degree);if(it==audit->end()||!same_state(it->second,next))throw logic_error("fixed/EEA quotient mismatch");st.audited_states++;}
   fixed_stats.checks++;s=move(next);f=P(u.begin()+g+1,u.end());
  }
  if(!pattern.empty())fixed_stats.patterns[pattern]++;int pending=end-(int)s.b.size()+1-completed;if(!pending)kind=0;
  physical=encode_fixed(s,end,completed,kind,history);
  physical=unpack_fixed(physical.packed,end,completed,physical.mode,history);State restored=decode_frame(physical);
  if(!same_state(restored,s))throw logic_error("pending frame reconstruction mismatch");
  if(audit){auto it=audit->find(s.b.size()-1);if(it==audit->end()||!same_state(it->second,restored))throw logic_error("fixed boundary/EEA mismatch");}
  fixed_stats.checks++;
 }
 if(physical.pending)throw runtime_error("pending survives final slot");return decode_frame(physical);
}
void fixed_report(ostream&o){
 o<<"{\"regions\":"<<fixed_stats.regions<<",\"intermediate_checks\":"<<fixed_stats.checks<<",\"chart_encodings\":"<<fixed_stats.encodings<<",\"pending_one\":"<<fixed_stats.pending1<<",\"pending_two\":"<<fixed_stats.pending2<<",\"deferred_beta\":"<<fixed_stats.deferred_beta<<",\"deferred_cubic_high\":"<<fixed_stats.deferred_cubic_high<<",\"gap_units\":"<<fixed_stats.gap_units<<",\"completion_slots\":"<<fixed_stats.completions<<",\"max_public_words\":"<<fixed_stats.max_words<<",\"max_fresh_head_words\":"<<fixed_stats.max_fresh<<",\"max_logical_head_words\":"<<fixed_stats.max_head<<",\"mode_hist\":[";
 for(int i=0;i<8;i++)o<<(i?",":"")<<fixed_stats.modes[i];o<<"],\"region_hist\":[";for(int i=0;i<36;i++)o<<(i?",":"")<<fixed_stats.boundaries[i];o<<"],\"patterns\":{";
 bool first=true;for(auto [p,n]:fixed_stats.patterns){o<<(first?"":",")<<'"'<<p<<"\":"<<n;first=false;}o<<"}}";
}
