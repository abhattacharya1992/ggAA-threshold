// Truncated series arithmetic.

#include "exact_rf.hpp"
#include <fstream>
#include <memory>
#include <chrono>
#include <iostream>
#include <climits>
struct JNode{int op,a,b;std::string c;};
using Jet=std::map<int,RF>;
struct Jets{
 RFEngine e;
 Jet add(Jet a,const Jet&b,int sign=1){for(auto&v:b){auto it=a.find(v.first);if(it==a.end()){auto q=v.second;if(sign<0)q=e.neg(q);if(!q.zero())a.emplace(v.first,std::move(q));}else{it->second=e.add(it->second,v.second,sign);if(it->second.zero())a.erase(it);}}return a;}
 Jet mul(const Jet&a,const Jet&b,int hi){Jet out;for(auto&u:a)for(auto&v:b)if(u.first+v.first<=hi){auto q=e.mul(u.second,v.second);auto it=out.find(u.first+v.first);if(it==out.end()){if(!q.zero())out.emplace(u.first+v.first,std::move(q));}else{it->second=e.add(it->second,q);if(it->second.zero())out.erase(it);}}return out;}
 Jet power(const Jet&a,int n,int hi,int lo){if(n<0)throw std::runtime_error("negative jet power");Jet out{{0,RF(1)}};for(int i=0;i<n;i++)out=mul(out,a,hi-(n-i-1)*lo);return out;}
 Jet poly(const Poly&p,int shift,int hi){Jet out;for(auto&v:p.t){int k=v.first/1024-shift;if(k>hi)continue;out[k].n.add(v.first%1024,v.second);}return out;}
 Jet invregular(const Poly&p,int shift,int hi){Jet d=poly(p,shift,hi);if(!d.count(0)||d.at(0).zero())throw std::runtime_error("zero regular denominator");Jet out;out[0]=e.inv(d.at(0));for(int k=1;k<=hi;k++){RF sum;for(auto&v:d)if(v.first>0&&v.first<=k){auto it=out.find(k-v.first);if(it!=out.end())sum=e.add(sum,e.mul(v.second,it->second));}auto a=e.neg(e.div(sum,d.at(0)));if(!a.zero())out.emplace(k,std::move(a));}return out;}
 static int pv(const Poly&p){if(p.zero())return INT_MAX/4;return p.t.begin()->first/1024;}
 static int valuation(const RF&r,const RFEngine&re){if(r.zero())return INT_MAX/4;int v=pv(r.n);for(auto&f:r.f)v-=f.second*pv(re.bases[f.first]);return v;}
 Jet rational(const RF&r,const RFEngine&re,int hi){if(r.zero())return{};int nv=pv(r.n),v=valuation(r,re),depth=hi-v;if(depth<0)return{};if(depth>100)throw std::runtime_error("required jet depth exceeds safety cap100");Jet q=poly(r.n,nv,depth);for(auto&f:r.f){Jet den=invregular(re.bases[f.first],pv(re.bases[f.first]),depth);q=mul(q,power(den,f.second,depth,0),depth);}RF scalar=e.inv(RF(r.d));Jet out;for(auto&u:q){auto rr=e.mul(u.second,scalar);if(!rr.zero())out.emplace(u.first+v,std::move(rr));}return out;}
};
int main(int argc,char**argv){try{
 if(argc!=4){std::cerr<<"usage jet_dag DAG OUTPUT HI\n";return 2;}
 int hi=std::stoi(argv[3]);if(hi<0||hi>24)throw std::runtime_error("unsupported order");
 std::ifstream in(argv[1]);if(!in)throw std::runtime_error("no DAG");std::string line,word;int root,md;
 std::getline(in,line);std::istringstream head(line);head>>word>>root>>word>>md;
 if(line.rfind("roots",0)==0){if(root!=1)throw std::runtime_error("jet backend requires unsplit DAG");std::getline(in,line);root=std::stoi(line);}
 std::vector<JNode> n{{}};while(std::getline(in,line)){std::istringstream is(line);int op,a=0,b=0;std::string v;is>>op>>v>>b;if(op)a=std::stoi(v);n.push_back({op,a,b,v});}
 size_t N=n.size();std::vector<bool> small(N,false);for(size_t i=N-1;i>0;i--){if(n[i].op==7&&n[i].b<0)small[i]=true;if(small[i]&&n[i].op>=3){small[n[i].a]=true;if(n[i].op<=6)small[n[i].b]=true;}}
 RFEngine re;RF b,c,one(1);b.n=Poly::var(1024);c.n=Poly::var(1);RF b2=re.mul(b,b),bc=re.mul(b,c);
 for(auto p:std::vector<Poly>{re.add(one,bc).n,re.add(bc,one,-1).n,re.add(re.add(one,b2),re.mul(RF(2),bc)).n,re.add(re.add(one,b2),re.mul(RF(2),bc),-1).n}){Big co=content(p);if(p.t.rbegin()->second<0)co=-co;re.reg(divscalar(p,co));}
 RF xx=re.div(re.add(one,b,-1),re.add(one,b)),zz=re.div(re.add(re.add(one,b2),re.mul(RF(2),bc)),re.add(one,b2,-1));
 std::vector<std::unique_ptr<RF>> exact(N);std::vector<int> lo(N,0),need(N,-10000),refs(N,0);auto start=std::chrono::steady_clock::now();size_t sc=0;
 for(size_t i=1;i<N;i++){auto q=n[i];if(small[i]){RF r;
  if(q.op==0)r=RF(Big(q.c));else if(q.op==1)r=xx;else if(q.op==2)r=zz;else if(q.op==3)r=re.add(*exact[q.a],*exact[q.b]);else if(q.op==4)r=re.add(*exact[q.a],*exact[q.b],-1);else if(q.op==5)r=re.mul(*exact[q.a],*exact[q.b]);else if(q.op==6)r=re.div(*exact[q.a],*exact[q.b]);else if(q.op==7)r=re.pow(*exact[q.a],q.b);else if(q.op==8)r=re.neg(*exact[q.a]);exact[i]=std::make_unique<RF>(std::move(r));sc++;}
  if(q.op==0)lo[i]=q.c=="0"?1000:0;else if(q.op<=2)lo[i]=0;else if(q.op==3||q.op==4)lo[i]=std::min(lo[q.a],lo[q.b]);else if(q.op==5)lo[i]=std::min(1000,lo[q.a]+lo[q.b]);else if(q.op==7){lo[i]=q.b<0?Jets::valuation(*exact[i],re):std::min(1000,lo[q.a]*q.b);}else if(q.op==8)lo[i]=lo[q.a];else throw std::runtime_error("un-normalized division");
 }
 std::cerr<<"small rational ancestors="<<sc<<" nodes="<<N-1<<" lower_bound="<<lo[root]<<" sec="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n";
 need[root]=hi;for(size_t i=N-1;i>0;i--){if(need[i]<lo[i])continue;auto q=n[i];if(q.op<3||(q.op==7&&q.b<0))continue;int a=need[i],bb=need[i];if(q.op==5){a-=lo[q.b];bb-=lo[q.a];}if(q.op==7)a-=(q.b-1)*lo[q.a];need[q.a]=std::max(need[q.a],a);refs[q.a]++;if(q.op<=6){need[q.b]=std::max(need[q.b],bb);refs[q.b]++;}}
 Jets j;std::vector<std::unique_ptr<Jet>> vals(N);Jet result;
 for(size_t i=1;i<N;i++){if(need[i]<lo[i]){vals[i]=std::make_unique<Jet>();continue;}auto q=n[i];Jet v;
 if(q.op==0){RF z(Big(q.c));if(!z.zero())v.emplace(0,std::move(z));}
 else if(q.op==1){v[0]=RF(1);for(int k=1;k<=need[i];k++)v[k]=RF(k%2?-2:2);}
 else if(q.op==2){v[0]=RF(1);for(int k=1;k<=need[i];k++){RF z(2);if(k%2)z.n=scale(Poly::var(1),2);v[k]=z;}}
 else if(q.op==7&&q.b<0)v=j.rational(*exact[i],re,need[i]);
 else if(q.op==3||q.op==4)v=j.add(*vals[q.a],*vals[q.b],q.op==3?1:-1);
 else if(q.op==5)v=j.mul(*vals[q.a],*vals[q.b],need[i]);
 else if(q.op==7)v=j.power(*vals[q.a],q.b,need[i],lo[q.a]);
 else if(q.op==8){v=*vals[q.a];for(auto&u:v)u.second=j.e.neg(u.second);}
 else throw std::runtime_error("bad opcode");
 for(auto it=v.begin();it!=v.end();)if(it->first>need[i])it=v.erase(it);else++it;
 vals[i]=std::make_unique<Jet>(std::move(v));
 if(q.op>=3&&!(q.op==7&&q.b<0)){if(--refs[q.a]==0)vals[q.a].reset();if(q.op<=6&&--refs[q.b]==0)vals[q.b].reset();}
 if(i==static_cast<size_t>(root))result=*vals[i];if(refs[i]==0)vals[i].reset();
 if(i%50000==0)std::cerr<<"jet nodes "<<i<<"/"<<N-1<<" sec="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<" bases="<<j.e.bases.size()<<"\n";
 }
 std::ofstream out(argv[2]);if(!out)throw std::runtime_error("cannot create output");out<<"<|\"Order\"->"<<hi<<",\"LowerBound\"->"<<lo[root]<<",\"Jet\"-><|";bool first=true;for(auto&v:result){if(!first)out<<",\n";first=false;out<<v.first<<"->"<<j.e.str(j.e.reduce(v.second),"UNUSED","costh");}out<<"|>|>\n";out.close();if(!out)throw std::runtime_error("output write failed");
 std::cerr<<"DONE sec="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<" coefficients="<<result.size()<<"\n";return 0;
}catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<"\n";return 1;}}
