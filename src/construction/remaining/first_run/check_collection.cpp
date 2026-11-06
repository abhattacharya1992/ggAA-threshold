// Checks the collected blocks against the source.
#include <array>
#include <sstream>
#include <filesystem>
#include <vector>
#include <map>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
using U=uint64_t;
constexpr U MOD=998244353;
U pw(U a, U n) {
    U r=1;
    for(; n; n>>=1, a=a*a%MOD)if(n&1)r=r*a%MOD;
    return r;
}
struct V {
    std::array<U, 3>a;
};
V same(U v) {
    return V {
        {
            v, v, v
        }
    };
}
V add(V x, V y) {
    for(int i=0; i<3; i++)x.a[i]=(x.a[i]+y.a[i])%MOD;
    return x;
}
V neg(V x) {
    for(U&v:x.a)v=v?MOD-v:0;
    return x;
}
V mul(V x, V y) {
    for(int i=0; i<3; i++)x.a[i]=x.a[i]*y.a[i]%MOD;
    return x;
}
V powv(V x, long n) {
    for(U&v:x.a) {
        if(n<0) {
            if(v==0)throw std::runtime_error("zero denominator");
            v=pw(v, MOD-2);
        }
        v=pw(v, n<0?-n:n);
    }
    return x;
}
V divv(V x, V y) {
    return mul(x, powv(y, -1));
}
std::string read(const char*p) {
    std::ifstream f(p);
    if(!f)throw std::runtime_error("unreadable input");
    std::string s((std::istreambuf_iterator<char>(f)), {
    }
    ), o;
    o.reserve(s.size());
    for(size_t i=0; i<s.size(); i++) {
        if(s[i]=='\\'&&i+1<s.size()&&s[i+1]=='\n') {
            i++;
            continue;
        }
        if(s[i]=='\\'&&i+2<s.size()&&s[i+1]=='\r'&&s[i+2]=='\n') {
            i+=2;
            continue;
        }
        o+=s[i];
    }
    size_t j=o.find('=');
    if(j!=std::string::npos)o.erase(0, j+1);
    return o;
}
class Eval {
    const std::string&s;
    size_t p=0;
    U seed;
    void ws() {
        while(p<s.size()&&isspace((unsigned char)s[p]))p++;
    }
    bool eat(char x) {
        ws();
        if(p<s.size()&&s[p]==x) {
            p++;
            return true;
        }
        return false;
    }
    U h(size_t a, size_t b) {
        U out=seed+41;
        for(size_t j=a; j<b; j++)if(!isspace((unsigned char)s[j]))out=(out*257+(unsigned char)s[j]+1)%MOD;
        return out?out:17;
    }
    V expression() {
        V x=term();
        for(; ; ) {
            ws();
            if(eat('+'))x=add(x, term());
            else if(eat('-'))x=add(x, neg(term()));
            else break;
        }
        return x;
    }
    V term() {
        V x=unary();
        for(; ; ) {
            if(eat('*'))x=mul(x, unary());
            else if(eat('/'))x=divv(x, unary());
            else break;
        }
        return x;
    }
    V unary() {
        if(eat('+'))return unary();
        if(eat('-'))return neg(unary());
        return power();
    }
    long exponent() {
        ws();
        if(eat('(')) {
            long x=exponent();
            if(!eat(')'))throw std::runtime_error("bad exponent )");
            return x;
        }
        long sign=1;
        if(eat('-'))sign=-1;
        else eat('+');
        ws();
        if(!isdigit((unsigned char)s[p]))throw std::runtime_error("noninteger power");
        long v=0;
        while(p<s.size()&&isdigit((unsigned char)s[p])) {
            v=v*10+s[p++]-'0';
            if(v>10000)throw std::runtime_error("large power");
        }
        return sign*v;
    }
    V power() {
        V x=atom();
        if(eat('^'))x=powv(x, exponent());
        return x;
    }
    V atom() {
        ws();
        size_t start=p;
        if(eat('(')) {
            V x=expression();
            if(!eat(')'))throw std::runtime_error("missing )");
            return x;
        }
        if(p<s.size()&&isdigit((unsigned char)s[p])) {
            U v=0;
            while(p<s.size()&&isdigit((unsigned char)s[p]))v=(10*v+s[p++]-'0')%MOD;
            return same(v);
        }
        if(p<s.size()&&isalpha((unsigned char)s[p])) {
            while(p<s.size()&&isalnum((unsigned char)s[p]))p++;
            std::string name=s.substr(start, p-start);
            ws();
            if(p<s.size()&&s[p]=='[') {
                int level=0;
                do {
                    char c=s[p++];
                    if(c=='[')level++;
                    if(c==']')level--;
                }
                while(p<s.size()&&level);
                if(level)throw std::runtime_error("unclosed func");
                return same(h(start, p));
            }
            if(name=="Nf")return V {
                {
                    MOD-1, 0, 1
                }
            };
            if(name=="I")return same(pw(3, (MOD-1)/4));
            return same(h(start, p));
        }
        throw std::runtime_error("bad atom near "+s.substr(p, 80));
    }
public:
    Eval(const std::string&ss, U se):s(ss), seed(se) {
    }
    V run() {
        V x=expression();
        ws();
        if(p!=s.size())throw std::runtime_error("trailing input");
        return x;
    }
};
int main(int argc, char**argv) {
    try {
        if(argc!=3) {
            std::cerr<<"Usage: check_collection INPUT BLOCK_DIR\n";
            return 2;
        }
        const std::string source=read(argv[1]);
        std::filesystem::path root(argv[2]);
        std::ifstream ff(root/"functions.txt");
        if(!ff)throw std::runtime_error("Missing function inventory");
        std::vector<std::string> funcs;
        std::string line;
        while(std::getline(ff, line))funcs.push_back(line);
        std::ifstream bf(root/"blocks.tsv");
        if(!bf)throw std::runtime_error("Missing block index");
        std::getline(bf, line);
        struct Block {
            int id, nf;
            std::vector<std::pair<int, int>>mon;
        };
        std::vector<Block>blocks;
        while(std::getline(bf, line)) {
            std::istringstream iss(line);
            std::string sid, snf, sm;
            std::getline(iss, sid, '\t');
            std::getline(iss, snf, '\t');
            std::getline(iss, sm);
            Block b {
                std::stoi(sid), std::stoi(snf), {
                }
            };
            std::istringstream mi(sm);
            std::string kv;
            while(std::getline(mi, kv, ',')) {
                auto pos=kv.find(':');
                if(pos==std::string::npos)throw std::runtime_error("Bad monomial");
                b.mon.push_back( {
                    std::stoi(kv.substr(0, pos)), std::stoi(kv.substr(pos+1))
                }
                );
            }
            blocks.push_back(b);
        }
        for(U seed: {
            19, 71, 113
        }
        ) {
            V reference=Eval(source, seed).run();
            std::array<U, 3> coeff {
                {
                    0, 0, 0
                }
            };
            std::vector<U>fv;
            for(auto&f:funcs)fv.push_back(Eval(f, seed).run().a[1]);
            for(auto&block:blocks) {
                char name[80];
                snprintf(name, sizeof(name), "rational_%05d.wl", block.id);
                std::string expr=read((root/name).c_str());
                V val=Eval(expr, seed).run();
                if(val.a[0]!=val.a[1]||val.a[1]!=val.a[2])throw std::runtime_error("Flavour in rational remainder");
                U v=val.a[1];
                for(auto&m:block.mon) {
                    if(m.first<0||m.first>=(int)fv.size())throw std::runtime_error("Bad atom id");
                    U p=powv(same(fv[m.first]), m.second).a[1];
                    v=v*p%MOD;
                }
                coeff[block.nf]=(coeff[block.nf]+v)%MOD;
            }
            U c0=reference.a[1], c1=(reference.a[2]+MOD-reference.a[0])%MOD*pw(2, MOD-2)%MOD;
            U c2=(reference.a[0]+reference.a[2]+2*MOD-2*reference.a[1])%MOD*pw(2, MOD-2)%MOD;
            bool ok=coeff[0]==c0&&coeff[1]==c1&&coeff[2]==c2;
            std::cout<<"seed="<<seed<<" complete polynomial reconstruction="<<ok<<" nf0="<<coeff[0]<<" nf1="<<coeff[1]<<" nf2="<<coeff[2]<<'\n';
            if(!ok)return 1;
        }
        return 0;
    } catch(const std::exception&e) {
        std::cerr<<"ERROR "<<e.what()<<'\n';
        return 1;
    }
}
