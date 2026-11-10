#include <fstream>
#include <iostream>
#include <string>
#include <cctype>
void append(std::ofstream &out,const std::string &file,int shift){
 std::ifstream in(file);if(!in)throw std::runtime_error("cannot open "+file);std::string line;
 while(std::getline(in,line)){
  if(line.rfind("REGISTER",0)==0||line.rfind("APPEND_TO_REGISTER",0)==0)continue;
  if(!shift){out<<line<<'\n';continue;}
  for(size_t i=0;i<line.size();){
   if(line[i]=='b'&&i+1<line.size()&&isdigit(line[i+1])&&(i==0||!isalnum(line[i-1]))){size_t j=i+1;int v=0;while(j<line.size()&&isdigit(line[j]))v=v*10+line[j++]-'0';out<<'b'<<v+shift;i=j;}
   else out<<line[i++];
  }out<<'\n';
 }
}
int main(int argc,char **argv){try{if(argc!=14)throw std::runtime_error("usage: out headers begin ef shift sf shift scan si shift ei shift end");std::ofstream out(argv[1]);std::ifstream header(argv[2]);out<<header.rdbuf();append(out,argv[3],0);append(out,argv[4],std::stoi(argv[5]));append(out,argv[6],std::stoi(argv[7]));append(out,argv[8],0);append(out,argv[9],std::stoi(argv[10]));append(out,argv[11],std::stoi(argv[12]));append(out,argv[13],0);if(!out)throw std::runtime_error("write failed");}catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}}
