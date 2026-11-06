// Rational blocks at fixed beta, exact in the angle.

#include "exact_rf.hpp"
#include <fstream>
#include <memory>
#include <chrono>
#include <iostream>
struct RefNode{int op,a,b;std::string c;};
int main(int argc,char**argv){try{
 if(argc!=5){std::cerr<<"reference_dag DAG OUTPUT betaNumerator betaDenominator\n";return 2;}
 long bn=std::stol(argv[3]),bd=std::stol(argv[4]);if(bn<=0||bd<=bn)throw std::runtime_error("bad beta");
 std::ifstream in(argv[1]);std::string line,word;int root,md;std::getline(in,line);std::istringstream hs(line);hs>>word>>root>>word>>md;
 if(line.rfind("roots",0)==0){if(root!=1)throw std::runtime_error("requires unsplit DAG");std::getline(in,line);root=std::stoi(line);}
 std::vector<RefNode> nodes{{}};std::vector<int> refs{0};while(std::getline(in,line)){std::istringstream ss(line);int op,a=0,b=0;std::string val;ss>>op>>val>>b;if(op)a=std::stoi(val);nodes.push_back({op,a,b,val});refs.push_back(0);if(op>=3){refs.at(a)++;if(op<=6)refs.at(b)++;}}
 RFEngine e;RF xx=e.div(RF(bd-bn),RF(bd+bn)),cc;cc.n=Poly::var(1);RF zz=e.div(e.add(RF(bd*bd+bn*bn),e.mul(RF(2*bn*bd),cc)),RF(bd*bd-bn*bn));
 std::vector<std::unique_ptr<RF>> vals(nodes.size());RF result;auto start=std::chrono::steady_clock::now();size_t maxterms=0;
 for(size_t i=1;i<nodes.size();i++){auto n=nodes[i];RF r;switch(n.op){case 0:r=RF(Big(n.c));break;case 1:r=xx;break;case 2:r=zz;break;case 3:r=e.add(*vals[n.a],*vals[n.b]);break;case 4:r=e.add(*vals[n.a],*vals[n.b],-1);break;case 5:r=e.mul(*vals[n.a],*vals[n.b]);break;case 6:r=e.div(*vals[n.a],*vals[n.b]);break;case 7:r=e.pow(*vals[n.a],n.b);break;case 8:r=e.neg(*vals[n.a]);break;default:throw std::runtime_error("bad DAG opcode");}
  maxterms=std::max(maxterms,r.n.t.size());if(i==static_cast<size_t>(root))result=r;vals[i]=std::make_unique<RF>(std::move(r));if(n.op>=3){if(--refs[n.a]==0)vals[n.a].reset();if(n.op<=6&&--refs[n.b]==0)vals[n.b].reset();}if(refs[i]==0)vals[i].reset();if(i%100000==0)std::cerr<<"reference nodes="<<i<<"/"<<nodes.size()-1<<" maxterms="<<maxterms<<" bases="<<e.bases.size()<<" sec="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n";
 }
 result=e.reduce(result);std::ofstream out(argv[2]);if(!out)throw std::runtime_error("cannot create reference output");out<<e.str(result,"UNUSED","costh")<<"\n";out.close();if(!out)throw std::runtime_error("reference output write failed");std::cerr<<"DONE seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<" terms="<<result.n.t.size()<<"\n";return 0;
}catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<"\n";return 1;}}
