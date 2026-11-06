// Checks the fixed-beta references at finite-field points.

#include <fstream>
#include <string>
#include <iostream>
#include <stdexcept>
#include <cctype>
#include <vector>
#include <cstdint>
using U=uint64_t;
U modpow(U a,U n,U p){U r=1;while(n){if(n&1)r=(__uint128_t)r*a%p;a=(__uint128_t)a*a%p;n>>=1;}return r;}
struct Point {U p,b,c,m,x,z;};
struct Parse {
 const std::string &s;size_t pos=0;Point v;
 void ws(){while(pos<s.size()&&std::isspace((unsigned char)s[pos]))pos++;}
 bool eat(char c){ws();if(pos<s.size()&&s[pos]==c){pos++;return true;}return false;}
 U add(U a,U b){return(a+b)%v.p;} U neg(U a){return a?v.p-a:0;} U mul(U a,U b){return(__uint128_t)a*b%v.p;}
 U power(U a,long n){if(n<0){if(!a)throw std::runtime_error("pole at check point");a=modpow(a,v.p-2,v.p);n=-n;}return modpow(a,n,v.p);}
 long exponent(){ws();if(eat('(')){long n=exponent();if(!eat(')'))throw std::runtime_error("bad exponent bracket");return n;}int sign=eat('-')?-1:1;eat('+');ws();long n=0;if(pos>=s.size()||!isdigit((unsigned char)s[pos]))throw std::runtime_error("noninteger exponent");while(pos<s.size()&&isdigit((unsigned char)s[pos])){n=n*10+s[pos++]-'0';if(n>1000)throw std::runtime_error("exponent too large");}return sign*n;}
 U expr(){U a=term();for(;;){if(eat('+'))a=add(a,term());else if(eat('-'))a=add(a,neg(term()));else return a;}}
 U term(){U a=unary();for(;;){if(eat('*'))a=mul(a,unary());else if(eat('/'))a=mul(a,power(unary(),-1));else return a;}}
 U unary(){if(eat('+'))return unary();if(eat('-'))return neg(unary());U a=atom();if(eat('^'))a=power(a,exponent());return a;}
 U atom(){ws();if(eat('(')){U a=expr();if(!eat(')'))throw std::runtime_error("bad bracket");return a;}
 if(pos<s.size()&&isdigit((unsigned char)s[pos])){U a=0;while(pos<s.size()&&isdigit((unsigned char)s[pos]))a=(10*a+s[pos++]-'0')%v.p;return a;}
 if(pos<s.size()&&isalpha((unsigned char)s[pos])){size_t first=pos;while(pos<s.size()&&isalnum((unsigned char)s[pos]))pos++;std::string t=s.substr(first,pos-first);if(t=="x")return v.x;if(t=="z")return v.z;if(t=="M2")return v.m;if(t=="beta")return v.b;if(t=="costh")return v.c;throw std::runtime_error("bad variable "+t);}
 throw std::runtime_error("bad atom at "+std::to_string(pos));}
 U all(){U r;if(eat('{')){r=0;if(!eat('}')){do{r=add(r,expr());}while(eat(','));if(!eat('}'))throw std::runtime_error("list not closed");}}else r=expr();ws();if(pos!=s.size())throw std::runtime_error("trailing text");return r;}
};
std::string read(const char*name){std::ifstream f(name);if(!f)throw std::runtime_error("cannot open input");return std::string((std::istreambuf_iterator<char>(f)),{});}
int main(int argc,char**argv){try{if(argc!=5){std::cerr<<"Usage: check_fixed_reference RAW_M2_BLOCK REFERENCE betaNumerator betaDenominator\n";return 2;}auto raw=read(argv[1]),reduced=read(argv[2]);std::vector<U> primes{998244353,1000000007};int count=0;
 for(U p:primes)for(U seed:{17,43,101}){Point v;v.p=p;v.b=(std::stoull(argv[3])*modpow(std::stoull(argv[4]),p-2,p))%p;v.c=((seed+29)*modpow(seed+113,p-2,p))%p;v.m=seed+7;v.x=(1+p-v.b)*modpow(1+v.b,p-2,p)%p;v.z=(1+v.b*v.b%p+2*v.b%p*v.c)%p*modpow((1+p-v.b*v.b%p)%p,p-2,p)%p;
 U a=Parse{raw,0,v}.all(),b=Parse{reduced,0,v}.all();bool ok=a==b*v.m%p;std::cout<<"prime="<<p<<" seed="<<seed<<" source="<<a<<" reduced_times_M2="<<b*v.m%p<<" pass="<<(ok?"true":"false")<<"\n";if(!ok)return 1;count++;}
 std::cout<<"CHECKS_PASSED="<<count<<" (algebraic sampling only; reduction itself uses exact integer arithmetic)\n";return 0;
}catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}}
