// Exact truncated Laurent series of rational blocks in beta.

#ifdef GGAA_USE_GMP
#include <gmpxx.h>
#else
#include <boost/multiprecision/cpp_int.hpp>
#endif
#include <fstream>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>
#include <string>
#include <cctype>
#include <climits>
#include <stdexcept>
#include <sstream>
#include <chrono>
#include <functional>
#ifdef GGAA_USE_GMP
using LargeInt=mpz_class;
std::string decimal(const LargeInt&v){return v.get_str();}
#else
using LargeInt=boost::multiprecision::cpp_int;
std::string decimal(const LargeInt&v){return v.convert_to<std::string>();}
#endif
struct Node {int op; int a=0,b=0; std::string c;};
struct Value {int id; int mass;};
struct DAG {
 std::vector<Node> nodes; std::map<std::string,int> constants;
 std::map<std::tuple<int,int,int>,int> memo; int zero,one,x,z;
 DAG(){zero=con("0");one=con("1");x=put({1});z=put({2});}
 int put(Node n){if(nodes.size()>8000000)throw std::runtime_error("DAG node cap exceeded");nodes.push_back(n);return nodes.size()-1;}
 int con(const std::string &s){std::string clean=s;size_t sign=(clean.size()>0&&clean[0]=='-')?1:0;while(clean.size()>sign+1&&clean[sign]=='0')clean.erase(sign,1);std::string v=decimal(LargeInt(clean));auto it=constants.find(v);if(it!=constants.end())return it->second;int k=put({0,0,0,v});constants[v]=k;return k;}
 bool iscon(int a){return nodes[a].op==0;}
 int inverse(int a){
   if(a==zero)throw std::runtime_error("division by zero");
   Node n=nodes[a];
   if(n.op==5)return make(5,inverse(n.a),inverse(n.b));
   if(n.op==7)return make(7,n.a,-n.b);
   if(n.op==8)return make(8,inverse(n.a));
   return make(7,a,-1);
 }
 int make(int op,int a,int b=0){
   if(op==6){if(b==zero)throw std::runtime_error("zero denominator");return make(5,a,inverse(b));}
   if(op==7&&b<0&&nodes[a].op==5){Node n=nodes[a];return make(5,make(7,n.a,b),make(7,n.b,b));}
   if(op==3||op==5){if(a>b)std::swap(a,b);}
   if(op==3){if(a==zero)return b;if(b==zero)return a;}
   if(op==4){if(b==zero)return a;if(a==b)return zero;}
   if(op==5){if(a==zero||b==zero)return zero;if(a==one)return b;if(b==one)return a;}
   if(op==6){if(b==zero)throw std::runtime_error("literal division by zero");if(b==one)return a;if(a==zero)return zero;}
   if(op==7){if(b==0){if(a==zero)throw std::runtime_error("undefined zero to zero power");return one;}if(b==1)return a; if(a==zero){if(b<0)throw std::runtime_error("zero to negative power");return zero;}if(a==one)return one;}
   if(op==8){if(a==zero)return zero;if(nodes[a].op==8)return nodes[a].a;}
   if(iscon(a) && (op==8 || (op!=7&&iscon(b)))){
    LargeInt aa(nodes[a].c), ans=0,bb=0;if(op!=8)bb=LargeInt(nodes[b].c);
    if(op==3)ans=aa+bb;else if(op==4)ans=aa-bb;else if(op==5)ans=aa*bb;else if(op==8)ans=-aa;else if(op==6){if(aa%bb!=0)goto skipcon;ans=aa/bb;}else goto skipcon;
    return con(decimal(ans));
   }
 skipcon:
   auto k=std::make_tuple(op,a,b);auto it=memo.find(k);if(it!=memo.end())return it->second;int i=put({op,a,b,{}});memo[k]=i;return i;
 }
};
class Parser {
 const std::string&s;size_t p=0;DAG&d;size_t depth=0;
 void ws(){while(p<s.size()&&isspace((unsigned char)s[p]))p++;}
 bool eat(char c){ws();if(p<s.size()&&s[p]==c){p++;return true;}return false;}
 [[noreturn]] void err(const std::string&m){throw std::runtime_error(m+" at byte "+std::to_string(p));}
 Value val(int id,int mass){return {id,id==d.zero?INT_MIN:mass};}
 Value add(Value a,Value b,int op){if(a.mass!=INT_MIN&&b.mass!=INT_MIN&&a.mass!=b.mass)err("mixed mass degrees");int m=a.mass==INT_MIN?b.mass:a.mass;return val(d.make(op,a.id,b.id),m);}
 Value mul(Value a,Value b,int op){if(b.id==d.zero&&op==6)err("zero divisor");int m=a.mass==INT_MIN||b.mass==INT_MIN?INT_MIN:a.mass+(op==5?b.mass:-b.mass);return val(d.make(op,a.id,b.id),m);}
 Value expression(){auto a=term();for(;;){if(eat('+'))a=add(a,term(),3);else if(eat('-'))a=add(a,term(),4);else break;}return a;}
 Value term(){auto a=unary();for(;;){if(eat('*'))a=mul(a,unary(),5);else if(eat('/'))a=mul(a,unary(),6);else break;}return a;}
 Value unary(){if(++depth>10000)err("parser nesting limit");Value a;if(eat('+'))a=unary();else if(eat('-')){a=unary();a.id=d.make(8,a.id);}else a=power();depth--;return a;}
 int exponent(){ws();if(eat('(')){int n=exponent();if(!eat(')'))err("exponent )");return n;}int sign=1;if(eat('-'))sign=-1;else eat('+');ws();if(p>=s.size()||!isdigit((unsigned char)s[p]))err("integer exponent required");int n=0;while(p<s.size()&&isdigit((unsigned char)s[p])){n=10*n+s[p++]-'0';if(n>1000)err("power limit");}return sign*n;}
 Value power(){auto a=atom();if(eat('^')){int n=exponent();a=val(d.make(7,a.id,n),a.mass==INT_MIN?INT_MIN:a.mass*n);}return a;}
 Value atom(){ws();if(eat('(')){auto a=expression();if(!eat(')'))err("missing )");return a;}if(p<s.size()&&isdigit((unsigned char)s[p])){size_t start=p;while(p<s.size()&&isdigit((unsigned char)s[p]))p++;return val(d.con(s.substr(start,p-start)),0);}if(p<s.size()&&isalpha((unsigned char)s[p])){size_t start=p;while(p<s.size()&&isalnum((unsigned char)s[p]))p++;auto name=s.substr(start,p-start);if(name=="x")return val(d.x,0);if(name=="z")return val(d.z,0);if(name=="M2")return val(d.one,1);err("unrecognized symbol "+name);}err("unexpected character");}
public:Parser(const std::string&t,DAG&g):s(t),d(g){}Value run(){auto a=expression();ws();if(p!=s.size())err("trailing input");return a;}
};
int main(int argc,char**argv){try{
 if(argc!=3){std::cerr<<"Usage: rational_dag SOURCE OUTPUT_PREFIX\n";return 2;}
 std::ifstream f(argv[1],std::ios::binary);if(!f)throw std::runtime_error("cannot read source");std::string s((std::istreambuf_iterator<char>(f)),{});DAG d;auto t=std::chrono::steady_clock::now();Value a=Parser(s,d).run();if(a.mass!=INT_MIN&&a.mass!=1)throw std::runtime_error("source mass degree is not one");
 std::vector<int> roots{a.id};
 std::vector<char> live(d.nodes.size(),0);std::vector<int> stack=roots;
 while(!stack.empty()){int i=stack.back();stack.pop_back();if(live[i])continue;live[i]=1;auto n=d.nodes[i];if(n.op>=3){stack.push_back(n.a);if(n.op<=6)stack.push_back(n.b);}}
 std::vector<int> ids(d.nodes.size(),-1);int count=0;for(size_t i=0;i<live.size();i++)if(live[i])ids[i]=++count;
 std::ofstream o(std::string(argv[2])+".dag.tsv");if(!o)throw std::runtime_error("cannot write DAG");
 o<<"roots\t"<<roots.size()<<"\tmass_degree\t"<<(a.mass==INT_MIN?0:a.mass)<<"\n";
 for(size_t i=0;i<roots.size();i++){if(i)o<<'\t';o<<ids[roots[i]];}o<<'\n';
 for(size_t i=0;i<live.size();i++)if(live[i]){auto n=d.nodes[i];o<<n.op<<'\t';if(n.op==0)o<<n.c<<"\t0";else if(n.op<=2)o<<"0\t0";else o<<ids[n.a]<<'\t'<<(n.op==7?n.b:(n.op==8?0:ids[n.b]));o<<'\n';}o.close();if(!o)throw std::runtime_error("DAG output write failed");
 std::cerr<<"roots="<<roots.size()<<" ";
 std::cerr<<"source_bytes="<<s.size()<<" interned_nodes="<<d.nodes.size()<<" live_nodes="<<count<<" mass_degree="<<a.mass<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count()<<"\n";return 0;
}catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}}
