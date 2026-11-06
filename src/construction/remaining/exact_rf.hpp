// Exact rational functions with big integers.

#pragma once
#ifndef GGAA_USE_GMP
#include <boost/multiprecision/cpp_int.hpp>
#endif
#include <map>
#include <vector>
#include <string>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <cstdint>
#include <cmath>
#ifdef GGAA_USE_GMP
#include <gmpxx.h>
using Big=mpz_class;
Big powbig(const Big&a,unsigned n){Big v;mpz_pow_ui(v.get_mpz_t(),a.get_mpz_t(),n);return v;}
#else
using Big=boost::multiprecision::cpp_int;
Big powbig(const Big&a,unsigned n){return boost::multiprecision::pow(a,n);}
#endif
Big absb(Big a){return a<0?-a:a;} Big gcdb(Big a,Big b){
#ifdef GGAA_USE_GMP
 Big v;mpz_gcd(v.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());return v;
#else
 a=absb(a);b=absb(b);while(b!=0){Big r=a%b;a=b;b=r;}return a;
#endif
}
struct Poly {
 std::map<int,Big> t;
 Poly(){} Poly(long a){if(a)t[0]=a;} Poly(const Big&a){if(a!=0)t[0]=a;}
 bool zero()const{return t.empty();} bool constant()const{return t.empty()||(t.size()==1&&t.begin()->first==0);}
 void add(int k,const Big &a){if(a==0)return;auto it=t.find(k);if(it==t.end())t.emplace(k,a);else{it->second+=a;if(it->second==0)t.erase(it);}}
 bool operator==(const Poly&o)const{return t==o.t;}
 bool operator<(const Poly&o)const{return t<o.t;}
 static Poly var(int k){Poly p;p.t[k]=1;return p;}
};
Poly scale(Poly p,const Big&a){if(a==0)return Poly();if(a==1)return p;for(auto&v:p.t)v.second*=a;return p;}
Poly addp(Poly p,const Poly&q,int sign=1){for(auto&v:q.t)p.add(v.first,sign*v.second);return p;}
Poly mulp(const Poly&a,const Poly&b){
 Poly p;if(a.zero()||b.zero())return p;
 if(a.constant())return scale(b,a.t.begin()->second);if(b.constant())return scale(a,b.t.begin()->second);
 int ax=a.t.rbegin()->first/1024,az=0,bx=b.t.rbegin()->first/1024,bz=0;
 for(auto&v:a.t)az=std::max(az,v.first%1024);for(auto&v:b.t)bz=std::max(bz,v.first%1024);
 int dx=ax+bx,dz=az+bz;if(dx>500||dz>500)throw std::runtime_error("polynomial degree cap exceeded");

#ifdef GGAA_USE_GMP
 if(a.t.size()*b.t.size()>100000){
  size_t abits=0,bbits=0;
  for(auto&v:a.t)abits=std::max(abits,mpz_sizeinbase(v.second.get_mpz_t(),2));
  for(auto&v:b.t)bbits=std::max(bbits,mpz_sizeinbase(v.second.get_mpz_t(),2));
  size_t count=std::min(a.t.size(),b.t.size()),carrybits=0;while(count){carrybits++;count>>=1;}
  size_t limbcount=(abits+bbits+carrybits+3+63)/64,bits=64*limbcount;
  size_t stride=dz+1, maxidx=dx*stride+dz;
  if((maxidx+2)*limbcount>64000000)throw std::runtime_error("packed multiplication exceeds 512MB cap");
  auto pack=[&](const Poly&q){Big z=0;size_t prev=0;bool first=true;
   for(auto it=q.t.rbegin();it!=q.t.rend();++it){size_t idx=(it->first/1024)*stride+it->first%1024;
    if(!first)mpz_mul_2exp(z.get_mpz_t(),z.get_mpz_t(),bits*(prev-idx));
    z+=it->second;prev=idx;first=false;
   }if(prev)mpz_mul_2exp(z.get_mpz_t(),z.get_mpz_t(),bits*prev);return z;
  };
  Big aa=pack(a),bb=pack(b),prod;mpz_mul(prod.get_mpz_t(),aa.get_mpz_t(),bb.get_mpz_t());
  int sign=mpz_sgn(prod.get_mpz_t());if(sign<0)prod=-prod;
  size_t len=(mpz_sizeinbase(prod.get_mpz_t(),2)+63)/64;
  std::vector<uint64_t> limbs(len+limbcount+1,0);size_t actual=0;
  mpz_export(limbs.data(),&actual,-1,sizeof(uint64_t),0,0,prod.get_mpz_t());
  Big base;mpz_setbit(base.get_mpz_t(),bits);Big half=base/2;int carry=0;
  for(size_t idx=0;idx<=maxidx;idx++){
   Big c;size_t offset=idx*limbcount;
   if(offset<actual)mpz_import(c.get_mpz_t(),std::min(limbcount,actual-offset),-1,sizeof(uint64_t),0,0,limbs.data()+offset);
   c+=carry;carry=0;if(c>=half){c-=base;carry=1;}
   if(c!=0){if(sign<0)c=-c;p.t.emplace((idx/stride)*1024+idx%stride,std::move(c));}
  }
  if(carry)throw std::runtime_error("packed multiplication leftover carry");
  return p;
 }
#endif
 int stride=dz+1;std::vector<Big> data((dx+1)*stride);
 for(auto&u:a.t)for(auto&v:b.t){int key=(u.first/1024+v.first/1024)*stride+u.first%1024+v.first%1024;
#ifdef GGAA_USE_GMP
 mpz_addmul(data[key].get_mpz_t(),u.second.get_mpz_t(),v.second.get_mpz_t());
#else
 data[key]+=u.second*v.second;
#endif
 }
 for(int x=0;x<=dx;x++)for(int z=0;z<=dz;z++)if(data[x*stride+z]!=0)p.t.emplace(x*1024+z,std::move(data[x*stride+z]));
 return p;
}
Poly powp(Poly a,int n){if(n<0||n>1000)throw std::runtime_error("bad polynomial power");Poly p(1);while(n){if(n&1)p=mulp(p,a);n>>=1;if(n)a=mulp(a,a);}return p;}
Big content(const Poly&p){Big g=0;for(auto&v:p.t){g=gcdb(g,v.second);if(g==1)break;}return g;}
Poly divscalar(Poly p,const Big&a){if(a==0)throw std::runtime_error("zero scalar divisor");for(auto&v:p.t){if(v.second%a!=0)throw std::runtime_error("nonexact scalar division");v.second/=a;}return p;}
bool dividep(const Poly&a,const Poly&b,Poly&result){
 if(b.zero())throw std::runtime_error("zero polynomial divisor");if(a.zero()){result=Poly();return true;}
 if(b.constant()){Big q=b.t.begin()->second;for(auto&v:a.t)if(v.second%q!=0)return false;result=divscalar(a,q);return true;}
 int ax=a.t.rbegin()->first/1024,az=0,bx=b.t.rbegin()->first/1024,bz=0;
 for(auto&v:a.t)az=std::max(az,v.first%1024);for(auto&v:b.t)bz=std::max(bz,v.first%1024);
 if(ax<bx||az<bz)return false;
 int qxmax=ax-bx,qzmax=az-bz,stride=az+1;
 auto lead=*b.t.rbegin();int lx=lead.first/1024,lz=lead.first%1024;
 std::vector<Big> rem((ax+1)*stride);for(auto&v:a.t)rem[(v.first/1024)*stride+v.first%1024]=v.second;
 int top=rem.size()-1;Poly quotient;
 while(top>=0){if(rem[top]==0){--top;continue;}int ex=top/stride-lx,ez=top%stride-lz;
  if(ex<0||ez<0||ex>qxmax||ez>qzmax)return false;
  Big c;
#ifdef GGAA_USE_GMP
  if(!mpz_divisible_p(rem[top].get_mpz_t(),lead.second.get_mpz_t()))return false;
  mpz_divexact(c.get_mpz_t(),rem[top].get_mpz_t(),lead.second.get_mpz_t());
#else
  if(rem[top]%lead.second!=0)return false;c=rem[top]/lead.second;
#endif
  quotient.add(ex*1024+ez,c);
  for(auto&v:b.t){int k=(ex+v.first/1024)*stride+ez+v.first%1024;
#ifdef GGAA_USE_GMP
   mpz_submul(rem[k].get_mpz_t(),c.get_mpz_t(),v.second.get_mpz_t());
#else
   rem[k]-=c*v.second;
#endif
  }
 }
 result=std::move(quotient);return true;
}
std::string polystr(const Poly&p,const std::string&x="x",const std::string&z="z"){
 if(p.zero())return "0";std::ostringstream s;bool first=true;
 for(auto it=p.t.rbegin();it!=p.t.rend();++it){Big c=it->second;if(!first&&c>=0)s<<'+';first=false;int xx=it->first/1024,zz=it->first%1024;
  if(c==-1&&(xx||zz))s<<'-';else if(c!=1||!(xx||zz))s<<c;
  bool coeff=absb(c)!=1||!(xx||zz);if(xx){if(coeff)s<<'*';s<<x;if(xx!=1)s<<'^'<<xx;}if(zz){if(xx||coeff)s<<'*';s<<z;if(zz!=1)s<<'^'<<zz;}
 }return s.str();
}
struct Factors{Big scalar=1;std::map<int,int> f;};
struct RF {Poly n;Big d=1;std::map<int,int> f;RF(){}RF(long a):n(a){}RF(const Big&a):n(a){}bool zero()const{return n.zero();}};
struct RFEngine {
 std::vector<Poly> bases;std::map<Poly,int> ids;std::map<Poly,Factors> factorcache;size_t nreduce=0;
 int reg(Poly p){auto it=ids.find(p);if(it!=ids.end())return it->second;int id=bases.size();bases.push_back(p);ids.emplace(std::move(p),id);if(bases.size()>5000)throw std::runtime_error("too many distinct denominator bases");return id;}
 RFEngine(){Poly x=Poly::var(1024),z=Poly::var(1),one(1),x2=mulp(x,x),xz=mulp(x,z);
  for(auto p:std::vector<Poly>{x,z,addp(x,one),addp(x,one,-1),addp(z,one),addp(z,one,-1),addp(x,z,-1),addp(xz,one,-1),addp(x2,one),addp(addp(x2,xz,-1),one),addp(addp(addp(x2,xz,-1),x),one),addp(addp(addp(x2,xz,-1),x,-1),one)})reg(p);
 }
 Factors factor(Poly p){if(p.zero())throw std::runtime_error("division by identically zero polynomial");auto it=factorcache.find(p);if(it!=factorcache.end())return it->second;Poly original=p;Factors out;out.scalar=content(p);if(p.t.rbegin()->second<0)out.scalar=-out.scalar;p=divscalar(p,out.scalar);
  if(!p.constant()){for(size_t i=0;i<bases.size();i++){Poly q;while(!p.constant()&&dividep(p,bases[i],q)){out.f[i]++;p=std::move(q);}if(p.constant())break;}}
  if(!p.constant())out.f[reg(p)]++;else if(!p.zero())out.scalar*=p.t.begin()->second;
  factorcache.emplace(std::move(original),out);return out;
 }
 RF reduce(RF a){nreduce++;if(a.zero())return RF(0);if(a.d<0){a.d=-a.d;a.n=scale(a.n,-1);}if(a.d==0)throw std::runtime_error("zero denominator");
  if(a.d!=1){Big g=gcdb(content(a.n),a.d);if(g!=1){a.n=divscalar(a.n,g);a.d/=g;}}
  for(auto it=a.f.begin();it!=a.f.end();){Poly q;while(it->second&&dividep(a.n,bases[it->first],q)){a.n=std::move(q);it->second--;}if(it->second==0)it=a.f.erase(it);else++it;}
  return a;
 }
 RF add(const RF&a,const RF&b,int sign=1){if(a.zero()){auto out=b;out.n=scale(out.n,sign);return out;}if(b.zero())return a;
  RF r;r.f=a.f;for(auto&v:b.f)r.f[v.first]=std::max(r.f[v.first],v.second);
  Big g=gcdb(a.d,b.d);r.d=(a.d/g)*b.d;Poly an=scale(a.n,r.d/a.d),bn=scale(b.n,sign*(r.d/b.d));
  for(auto&v:r.f){auto ai=a.f.find(v.first),bi=b.f.find(v.first);int da=v.second-(ai==a.f.end()?0:ai->second),db=v.second-(bi==b.f.end()?0:bi->second);if(da)an=mulp(an,powp(bases[v.first],da));if(db)bn=mulp(bn,powp(bases[v.first],db));}
  r.n=addp(std::move(an),bn);return reduce(std::move(r));
 }
 void cancelAcross(Poly &num, Big &scalar, std::map<int,int> &factors){
 if(scalar!=1){Big g=gcdb(content(num),scalar);if(g!=1){num=divscalar(num,g);scalar/=g;}}
 for(auto it=factors.begin();it!=factors.end();){Poly q;while(it->second&&dividep(num,bases[it->first],q)){num=std::move(q);it->second--;}if(it->second==0)it=factors.erase(it);else ++it;}
 }
 RF mul(const RF&a,const RF&b){if(a.zero()||b.zero())return RF(0);RF aa=a,bb=b;cancelAcross(aa.n,bb.d,bb.f);cancelAcross(bb.n,aa.d,aa.f);RF r;r.n=mulp(aa.n,bb.n);r.d=aa.d*bb.d;r.f=aa.f;for(auto&v:bb.f)r.f[v.first]+=v.second;return r;}

 RF inv(const RF&a){Factors f=factor(a.n);RF r;r.n=Poly(a.d);for(auto&v:a.f)r.n=mulp(r.n,powp(bases[v.first],v.second));r.d=f.scalar;r.f=f.f;return reduce(std::move(r));}
 RF div(const RF&a,const RF&b){if(b.zero())throw std::runtime_error("zero rational divisor");if(a.zero())return RF(0);
 Factors bf=factor(b.n);RF r=a;auto numeratorFactors=b.f;
 for(auto &v:numeratorFactors){auto it=r.f.find(v.first);if(it!=r.f.end()){int n=std::min(it->second,v.second);it->second-=n;v.second-=n;if(it->second==0)r.f.erase(it);}}
 cancelAcross(r.n,bf.scalar,bf.f);
 Big g=gcdb(r.d,b.d);r.n=scale(r.n,b.d/g);r.d/=g;r.d*=bf.scalar;
 for(auto &v:numeratorFactors)if(v.second)r.n=mulp(r.n,powp(bases[v.first],v.second));
 for(auto &v:bf.f)r.f[v.first]+=v.second;return reduce(std::move(r));
 }

 RF pow(const RF&a,int n){if(n==0)return RF(1);if(n<0)return pow(inv(a),-n);RF r;r.n=powp(a.n,n);r.d=powbig(a.d,n);for(auto&v:a.f)r.f[v.first]=v.second*n;return r;}
 RF neg(const RF&a){RF r=a;r.n=scale(r.n,-1);return r;}
 std::string str(const RF&a,const std::string&x="x",const std::string&z="z")const{std::ostringstream s;s<<'('<<polystr(a.n,x,z)<<")/("<<a.d;for(auto&v:a.f){s<<"*("<<polystr(bases[v.first],x,z)<<')';if(v.second!=1)s<<'^'<<v.second;}s<<')';return s.str();}
};
