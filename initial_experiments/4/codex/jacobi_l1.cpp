// Exact one-dimensional search for squared column L1 normalization.
// Orthogonal column rotations preserve A A^T, including the identity offset.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>
using namespace std;
struct Event {double theta; int row; bool first;};
struct Rotation {double theta,gain;};
static int sign(double x){return (x>0)-(x<0);}
static Rotation optimum(const vector<double>& A,int n,int c,int j,int k) {
    double a=0,b=0,cc=0,d=0;
    vector<Event> events;events.reserve(n);
    for(int i=0;i<n;i++) {
        double x=A[i*c+j],y=A[i*c+k];
        int sx=x?sign(x):sign(y),sy=y?sign(y):-sign(x);
        a+=sx*x;b+=sx*y;cc+=sy*y;d-=sy*x;
        if(x*y<0)events.push_back({atan2(abs(x),abs(y)),i,true});
        else if(x*y>0)events.push_back({atan2(abs(y),abs(x)),i,false});
    }
    sort(events.begin(),events.end(),[](const Event& x,const Event& y){return x.theta<y.theta;});
    double initial=a*a+cc*cc,best=initial,angle=0,left=0;
    auto segment=[&](double right) {
        double u=.5*(a*a+cc*cc-b*b-d*d),v=a*b+cc*d;
        double t=.5*(atan2(v,u)+acos(-1.));
        auto check=[&](double theta){
            double co=cos(theta),si=sin(theta);
            double xx=a*co+b*si,yy=cc*co+d*si,value=xx*xx+yy*yy;
            if(value<best){best=value;angle=theta;}
        };
        check(right);
        if(t>left&&t<right)check(t);
    };
    for(const auto& e:events) {
        segment(e.theta);left=e.theta;
        double x=A[e.row*c+j],y=A[e.row*c+k];
        if(e.first){int sx=sign(x);a-=2*sx*x;b-=2*sx*y;}
        else {int sy=sign(y);cc-=2*sy*y;d+=2*sy*x;}
    }
    segment(.5*acos(-1.));
    return {angle,initial-best};
}
static double objective(const vector<double>& A,int n,int c) {
    vector<double> sums(c,0);
    for(int i=0;i<n;i++)for(int j=0;j<c;j++)sums[j]+=abs(A[i*c+j]);
    double result=0;for(double x:sums)result+=x*x;
    return .25*result;
}
static double apply(vector<double>& A,int n,int c,int j,int k,double theta) {
    double co=cos(theta),si=sin(theta),oldx=0,oldy=0,newx=0,newy=0;
    for(int i=0;i<n;i++){
        double x=A[i*c+j],y=A[i*c+k];oldx+=abs(x);oldy+=abs(y);
        newx+=abs(co*x+si*y);newy+=abs(-si*x+co*y);
    }
    double gain=oldx*oldx+oldy*oldy-newx*newx-newy*newy;
    if(gain<=1e-13)return 0;
    for(int i=0;i<n;i++){
        double x=A[i*c+j],y=A[i*c+k];
        A[i*c+j]=co*x+si*y;A[i*c+k]=-si*x+co*y;
    }
    return .25*gain;
}
static void save(const string& file,const vector<double>& A,uint32_t n,uint32_t c) {
    ofstream f(file,ios::binary);f.write((char*)&n,4);f.write((char*)&c,4);
    f.write((char*)A.data(),A.size()*8);
}
int main(int argc,char**argv) {
    if(argc<5){cerr<<"input output passes pairs_per_pass [seed]\n";return 2;}
    uint32_t n,c;ifstream f(argv[1],ios::binary);f.read((char*)&n,4);f.read((char*)&c,4);
    vector<double>A(size_t(n)*c);f.read((char*)A.data(),A.size()*8);if(!f)return 3;
    int passes=stoi(argv[3]),maxpairs=stoi(argv[4]);mt19937_64 rng(argc>5?stoull(argv[5]):791291);
    vector<pair<int,int>> pairs;
    for(int j=0;j<(int)c;j++)for(int k=0;k<j;k++)pairs.push_back({j,k});
    cout<<setprecision(17)<<"{\"pass\":-1,\"lambda_h2\":"<<objective(A,n,c)<<"}"<<endl;
    for(int p=0;p<passes;p++) {
        shuffle(pairs.begin(),pairs.end(),rng);int accepted=0;double gain=0;
        int count=maxpairs>0?min(maxpairs,(int)pairs.size()):(int)pairs.size();
        for(int h=0;h<count;h++) {
            auto[j,k]=pairs[h];auto r=optimum(A,n,c,j,k);
            if(r.gain>1e-12){double g=apply(A,n,c,j,k,r.theta);gain+=g;accepted+=g>0;}
        }
        save(argv[2],A,n,c);
        cout<<"{\"pass\":"<<p<<",\"lambda_h2\":"<<objective(A,n,c)
            <<",\"accepted\":"<<accepted<<",\"gain\":"<<gain<<"}"<<endl;
        if(gain<1e-9)break;
    }
}
