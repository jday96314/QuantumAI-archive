// Independent totalized prefix/counter decoder model. No quantum acceptance.
#define main generalized_model_main
#include "generalized_quotient_jacobi_model.cpp"
#undef main

V truncprod(const V&a,const V&b,int n){V c(n);for(int i=0;i<(int)a.size()&&i<n;i++)for(int j=0;j<(int)b.size()&&i+j<n;j++)c[i+j]^=mul(a[i],b[j]);return c;}
V pair_model(V f){int a=f[0];V norm=f;norm[0]=1;int ai=inv(a?a:1);for(int j=1;j<(int)f.size();j++)norm[j]=mul(f[j],ai);V out=series({1},norm,f.size());out[0]=a;return out;}
V head_model(V h){for(int j=0;j<(int)h.size();j+=2){V x=pair_model(V(h.begin()+j,h.end()));std::copy(x.begin(),x.end(),h.begin()+j);}return h;}
struct Frame{V A,D;int p;};
Frame frame_model(const V&H,bool early){
 if(!early){V old={1},cur={H[1],1};int p=H[0];for(int j=2;j<(int)H.size();j+=2){int a=H[j],b=H[j+1];if(a){V out=prod({b,1},cur);for(int i=0;i<(int)old.size();i++)out[i]^=mul(a,old[i]);old=cur;cur=out;p=mul(p,a);}else {old.insert(old.begin(),b);cur.insert(cur.begin(),0);}}V A={1},D={1};A.insert(A.end(),cur.rbegin()+1,cur.rend());D.insert(D.end(),old.rbegin()+1,old.rend());return {A,D,p};}
 V P,B={H[1]};int p=H[0];
 for(int j=2;j<(int)H.size();j+=2){int beta=H[j],alpha=H[j+1];V out(P.size()+2);
  if(early||beta){int be=beta?beta:1;p=mul(p,be);for(int i=0;i<(int)P.size();i++)out[i]=mul(be,P[i])^mul(alpha,B[i])^(i?B[i-1]:0);out[P.size()]=beta^mul(alpha,B.back())^(B.size()>1?B[B.size()-2]:0);out[P.size()+1]=alpha^B.back();P=B;B=out;}
  else {out[0]=alpha;std::copy(P.begin(),P.end(),out.begin()+1);P=B;B=out;}
 }
 V A={1},D={1};A.insert(A.end(),B.rbegin(),B.rend());D.insert(D.end(),P.rbegin(),P.rend());return {A,D,p};
}
struct Extraction{V W;bool early_zero=false,diagonal_zero=false;};
Extraction prefix_extract(V W,int prefix,int through,int suffix){Extraction result;int k=0;
 while(k<(int)W.size()/2){int b=k<through?std::min({prefix,through-k,int(W.size()/2)-k}):std::min(suffix,int(W.size()/2)-k);int start=2*k,n=W.size()-start-2*b;
  if(b==1){V v=pair_model(V(W.begin()+start,W.end()));std::copy(v.begin(),v.end(),W.begin()+start);k++;continue;}
  V raw(W.begin()+start,W.begin()+start+2*b),H=head_model(raw),R(W.begin()+start+2*b,W.end());
  if(n){bool early=k<through;Frame f=frame_model(H,early);V P(b);if(early){for(int j=0;j<2*b;j+=2)result.early_zero|=!H[j];result.diagonal_zero|=!f.A.back();}
   for(int i=0;i<b;i++){P[i]=mul(early?(f.A.back()?f.A.back():1):f.A.back(),raw[b+i]);for(int z=i+1;z<b;z++)P[i]^=mul(f.A[b+i-z],raw[b+z]);}
   R=truncprod(f.A,R,n);for(int j=0;j<std::min(b,n);j++)R[j]^=P[j];V den={f.p?f.p:1,0};V dr=truncprod(f.D,R,std::max(0,n-2));den.insert(den.end(),dr.begin(),dr.end());R=series(truncprod(f.A,R,n),den,n);
  }
  std::copy(H.begin(),H.end(),W.begin()+start);std::copy(R.begin(),R.end(),W.begin()+start+2*b);k+=b;
 }
 result.W=W;return result;
}
std::pair<V,V> counter_put(V W){int n=W.size()/2,zeros=0,prev=1;for(int i=0;i<n;i++)zeros+=W[2*i]==0;
 for(int i=0;i<n;i++){int a=W[2*i],b=W[2*i+1];int z=mul(a,prev?prev:1);int k=z?inv(z):0;W[2*i]=k;W[2*i+1]=mul(b,k?k:1);prev=k;}
 int d=n-zeros;if(zeros&&W[2*d]==0)W[2*d]=1;V F={zeros?0:1},G;
 for(int i=n-1;i>=0;i--){int q0=W[2*i+1],q1=W[2*i];V out(F.size()+1);for(int j=0;j<(int)G.size();j++)out[j]=G[j]^mul(q0,F[j])^(j?mul(q1,F[j-1]):0);out[out.size()-2]=mul(q0,F.back()?F.back():1)^(F.size()>1?mul(q1,F[F.size()-2]):0);out.back()=mul(q1,F.back()?F.back():1);G=F;F=out;}
 return {F,G};
}
bool correct(const V&F,const V&G,const V&pos,const V&vals){V D(F.size()-1);for(int i=1;i<(int)F.size();i+=2)D[i-1]=F[i];
 for(int j=0;j<(int)pos.size();j++){int a=ex[pos[j]],den=mul(a,eval(D,a));if(eval(F,a)||!den||mul(eval(G,a),inv(den))!=vals[j])return false;}
 // If all true roots and residues pass, ensure no other nonzero field root
 // contributes a spurious residue. Usually common t padding proves this.
 int lo=0,hi=F.size()-1;while(lo<hi&&!F[lo])lo++;while(hi>lo&&!F[hi])hi--;if(hi-lo==(int)pos.size())return true;
 for(int j=0;j<4095;j++)if(std::find(pos.begin(),pos.end(),j)==pos.end()){int a=ex[j],den=mul(a,eval(D,a));if(!eval(F,a)&&den&&eval(G,a))return false;}return true;
}
#ifndef PREFIX_COUNTER_NO_MAIN
int main(int argc,char**argv){int shots=argc>1?std::stoi(argv[1]):100000,prefix=argc>2?std::stoi(argv[2]):6,through=argc>3?std::stoi(argv[3]):18,suffix=argc>4?std::stoi(argv[4]):2;uint64_t seed=2026102612;lg.fill(-1);int x=1;for(int i=0;i<4095;i++){ex[i]=x;lg[x]=i;x<<=1;if(x&4096)x^=0x1053;}for(int i=4095;i<8190;i++)ex[i]=ex[i-4095];std::ifstream file("parameters/dqi_weight_probabilities.txt");std::vector<double>prob;double z;while(file>>z)prob.push_back(z);if(prob.size()!=36)return 2;
 std::mt19937_64 rng(seed);std::discrete_distribution<int>weight(prob.begin(),prob.end());std::uniform_int_distribution<int>field(1,4095),position(0,4094);int wrong=0,zero=0,diag=0;
 for(int shot=0;shot<shots;shot++){int d=weight(rng);V pos,vals;while((int)pos.size()<d){int pp=position(rng);if(std::find(pos.begin(),pos.end(),pp)==pos.end())pos.push_back(pp);}for(int i=0;i<d;i++)vals.push_back(field(rng));auto r=prefix_extract(moments(pos,vals),prefix,through,suffix);auto [F,G]=counter_put(r.W);wrong+=!correct(F,G,pos,vals);zero+=r.early_zero;diag+=r.diagonal_zero;}
 int degree_wrong=0;for(int d=0;d<=35;d++)for(int trial=0;trial<8;trial++){V pos,vals;while((int)pos.size()<d){int pp=position(rng);if(std::find(pos.begin(),pos.end(),pp)==pos.end())pos.push_back(pp);}for(int i=0;i<d;i++)vals.push_back(field(rng));auto r=prefix_extract(moments(pos,vals),prefix,through,suffix);auto [F,G]=counter_put(r.W);degree_wrong+=!correct(F,G,pos,vals);}
 std::cout<<"{\"status\":\"CLASSICAL_ERROR_CENSUS\",\"law_shots\":"<<shots<<",\"seed\":"<<seed<<",\"prefix\":"<<prefix<<",\"through\":"<<through<<",\"suffix\":"<<suffix<<",\"wrong_vectors\":"<<wrong<<",\"early_zero_vectors\":"<<zero<<",\"early_diagonal_zero_vectors\":"<<diag<<",\"completion_fixture_cases\":288,\"completion_fixture_wrong_vectors\":"<<degree_wrong<<",\"scope\":\"independent totalized rational map and all nonzero-field roots/magnitudes; unsupported chart errors counted, no quantum acceptance\"}"<<std::endl;
}
#endif
