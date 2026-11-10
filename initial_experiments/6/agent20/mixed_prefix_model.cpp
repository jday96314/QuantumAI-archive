// Independent whole-law census of the totalized native mixed decoder.
#define PREFIX_COUNTER_NO_MAIN
#include "prefix_counter_model.cpp"
int square(int a){return mul(a,a);}
V look_block(V W){bool flag=!W[0];V H(W.begin(),W.begin()+4),R(W.begin()+4,W.end());int n=R.size();
 if(flag){int bi=inv(H[1]?H[1]:1),e=mul(H[2],bi);H={H[0],H[1],mul(H[3],bi)^square(e),e};}
 else {H=head_model(H);int a=H[0],b=H[1],c=H[2],d=H[3];H={a,b,c^mul(b,d),b^d};}
 int a=H[0],b=H[1],h=H[2],e=H[3],c=h^mul(b,e)^square(b),d=e^b,bg=flag?b:0,p=flag?b:mul(a,c?c:1);V A={1,e,h},D={1,b^bg};int p0=mul(a,square(square(b))^square(c)^mul(c,square(e)))^mul(bg,mul(square(e),e));int p1=mul(h,mul(a,mul(square(b),b)^mul(c,d))^mul(bg,h^square(e)));
 if(n){R=truncprod(A,R,n);R[0]^=p0;if(n>1)R[1]^=p1;V den={p?p:1};den.resize(2+flag);V dr=truncprod(D,R,std::max(0,n-2-int(flag)));den.insert(den.end(),dr.begin(),dr.end());R=series(truncprod(A,R,n),den,n);}H.insert(H.end(),R.begin(),R.end());return H;
}
V hybrid_extract(V W,int prefix,int through){int k=0,n=W.size()/2;while(k<n){int b=k<through?std::min({prefix,through-k,n-k}):std::min(2,n-k);int start=2*k;V tail(W.begin()+start,W.end());
 if(b==1)tail=pair_model(tail);
 else if(k>=through)tail=look_block(tail);
 else {V raw(tail.begin(),tail.begin()+2*b),H=head_model(raw),R(tail.begin()+2*b,tail.end());int nn=R.size();if(nn){auto f=frame_model(H,true);V P(b);for(int i=0;i<b;i++){P[i]=mul(f.A.back()?f.A.back():1,raw[b+i]);for(int z=i+1;z<b;z++)P[i]^=mul(f.A[b+i-z],raw[b+z]);}R=truncprod(f.A,R,nn);for(int i=0;i<std::min(b,nn);i++)R[i]^=P[i];V den={f.p?f.p:1,0};V dr=truncprod(f.D,R,std::max(0,nn-2));den.insert(den.end(),dr.begin(),dr.end());R=series(truncprod(f.A,R,nn),den,nn);}tail=H;tail.insert(tail.end(),R.begin(),R.end());}
 std::copy(tail.begin(),tail.end(),W.begin()+start);k+=b;
 }return W;}
std::pair<V,V> put_one(int q0,int q1,const V&f,const V&g){V out(f.size()+1);for(int j=0;j<(int)g.size();j++)out[j]=g[j]^mul(q0,f[j])^(j?mul(q1,f[j-1]):0);out[out.size()-2]=mul(q0,f.back()?f.back():1)^(f.size()>1?mul(q1,f[f.size()-2]):0);out.back()=mul(q1,f.back()?f.back():1);return {out,f};}
std::pair<V,V> hybrid_put(V W,int through){int n=W.size()/2,prev=1,zeros=0;std::vector<bool> flags;
 for(int k=0;k<n-1;k+=2){int a=W[2*k],b=W[2*k+1],h=W[2*k+2],e=W[2*k+3];bool flag=k>=through&&!a;flags.push_back(flag);int c=k>=through?(flag?e:h^mul(b,e^b)):h,d=k>=through?(flag?h:e^b):e;
  if(flag){zeros+=2*(!b);int v=mul(b,prev?prev:1),q=v?inv(v):0;W[2*k]=q;W[2*k+1]=mul(c,q?q:1);W[2*k+2]=0;W[2*k+3]=mul(d,q?q:1);prev=q;}
  else {zeros+=(!a)+(!c);int v=mul(a,prev?prev:1),q=v?inv(v):0;v=mul(c,q?q:1);int r=v?inv(v):0;W[2*k]=q;W[2*k+1]=mul(b,q?q:1);W[2*k+2]=r;W[2*k+3]=mul(d,r?r:1);prev=r;}
 }
 if(n%2){int a=W[2*n-2],b=W[2*n-1];zeros+=!a;int v=mul(a,prev?prev:1),q=v?inv(v):0;W[2*n-2]=q;W[2*n-1]=mul(b,q?q:1);}
 int degree=n-zeros;if(zeros&&!W[2*degree])W[2*degree]=1;V F={zeros?0:1},G;
 if(n%2)std::tie(F,G)=put_one(W[2*n-1],W[2*n-2],F,G);
 for(int k=n/2-1;k>=0;k--){int q2=W[4*k],q1=W[4*k+1],r1=W[4*k+2],q0=W[4*k+3];if(!flags[k]){std::tie(F,G)=put_one(q0,r1,F,G);std::tie(F,G)=put_one(q1,q2,F,G);}else {V gg=F;gg.push_back(0);V out=G;out.insert(out.end(),{q0,q1,q2});int nn=gg.size();for(int i=0;i<nn-2;i++){int q[3]={q0,q1,q2};for(int j=0;j<3;j++)if(i>=j)out[i]^=mul(q[j],gg[i-j]);}int lead=gg[nn-2]?gg[nn-2]:1;out[nn-2]=mul(q0,lead)^(nn>=3?mul(q1,gg[nn-3]):0)^(nn>=4?mul(q2,gg[nn-4]):0);out[nn-1]=mul(q1,lead)^(nn>=3?mul(q2,gg[nn-3]):0);out[nn]=mul(q2,lead);F=out;G=gg;}}
 return {F,G};}
int main(int argc,char**argv){lg.fill(-1);int x=1;for(int i=0;i<4095;i++){ex[i]=x;lg[x]=i;x<<=1;if(x&4096)x^=0x1053;}for(int i=4095;i<8190;i++)ex[i]=ex[i-4095];
 if(argc>1&&std::string(argv[1])=="--fixtures"){int n,p,t;while(std::cin>>n>>p>>t){V W(n);for(int&v:W)std::cin>>v;auto [F,G]=hybrid_put(hybrid_extract(W,p,t),t);for(int v:F)std::cout<<v<<' ';for(int v:G)std::cout<<v<<' ';std::cout<<'\n';}return 0;}
 int shots=argc>1?std::stoi(argv[1]):100000,prefix=argc>2?std::stoi(argv[2]):6,through=argc>3?std::stoi(argv[3]):18;uint64_t seed=2026102621;std::ifstream file("parameters/dqi_weight_probabilities.txt");std::vector<double>prob;double z;while(file>>z)prob.push_back(z);if(prob.size()!=36)return 2;std::mt19937_64 rng(seed);std::discrete_distribution<int>weight(prob.begin(),prob.end());std::uniform_int_distribution<int>field(1,4095),position(0,4094);int wrong=0,regular_wrong=0;
 auto check=[&](int d){V pos,vals;while((int)pos.size()<d){int pp=position(rng);if(std::find(pos.begin(),pos.end(),pp)==pos.end())pos.push_back(pp);}for(int i=0;i<d;i++)vals.push_back(field(rng));auto [F,G]=hybrid_put(hybrid_extract(moments(pos,vals),prefix,through),through);return correct(F,G,pos,vals);};
 for(int i=0;i<shots;i++)wrong+=!check(weight(rng));for(int d=0;d<=35;d++)for(int j=0;j<8;j++)regular_wrong+=!check(d);
 std::cout<<"{\"status\":\"CLASSICAL_ERROR_CENSUS\",\"law_shots\":"<<shots<<",\"seed\":"<<seed<<",\"prefix\":"<<prefix<<",\"through\":"<<through<<",\"wrong_vectors\":"<<wrong<<",\"completion_fixture_cases\":288,\"completion_fixture_wrong_vectors\":"<<regular_wrong<<",\"scope\":\"independent whole-field root and magnitude test of native totalized mixed decoder; not native task acceptance\"}"<<std::endl;
}
