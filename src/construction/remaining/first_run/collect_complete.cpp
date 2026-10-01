// Collects the source expressions into blocks.
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <filesystem>
namespace fs=std::filesystem;
struct Expr;
using Ptr=std::shared_ptr<Expr>;
struct Expr {
    char kind;
    size_t a=0, b=0;
    std::vector<Ptr> args;
    explicit Expr(char k):kind(k) {
    }
};
Ptr unit() {
    static Ptr p=std::make_shared<Expr>('1');
    return p;
}
Ptr raw(size_t a, size_t b) {
    auto p=std::make_shared<Expr>('r');
    p->a=a;
    p->b=b;
    return p;
}
Ptr minus(Ptr a) {
    if(!a)return {
    };
    if(a->kind=='-')return a->args[0];
    auto p=std::make_shared<Expr>('-');
    p->args= {
        a
    };
    return p;
}
Ptr plus(Ptr a, Ptr b) {
    if(!a)return b;
    if(!b)return a;
    auto p=std::make_shared<Expr>('+');
    p->args= {
        a, b
    };
    return p;
}
Ptr times(Ptr a, Ptr b) {
    if(!a||!b)return {
    };
    if(a->kind=='1')return b;
    if(b->kind=='1')return a;
    auto p=std::make_shared<Expr>('*');
    p->args= {
        a, b
    };
    return p;
}
Ptr divide(Ptr a, Ptr b) {
    if(!a)return {
    };
    if(!b)throw std::runtime_error("Zero denominator");
    if(b->kind=='1')return a;
    auto p=std::make_shared<Expr>('/');
    p->args= {
        a, b
    };
    return p;
}
Ptr ipow(Ptr a, int n) {
    if(n==0)return unit();
    if(n<0)return divide(unit(), ipow(a, -n));
    Ptr out=unit();
    while(n) {
        if(n&1)out=times(out, a);
        n>>=1;
        if(n)a=times(a, a);
    }
    return out;
}
using Key=std::map<int, int>;
using Poly=std::map<Key, Ptr>;
struct Node {
    size_t a, b;
    bool special=false;
    Poly poly;
};
void normalize(Key &k, Ptr &c) {
    int ii=k.count(3)?k[3]:0;
    int r=((ii%4)+4)%4;
    if(r>=2) {
        c=minus(c);
        r-=2;
    }
    k[3]=r;
    for(auto it=k.begin(); it!=k.end(); )if(it->second==0)it=k.erase(it);
    else ++it;
}
class Parser {
public:
    const std::string&s;
    size_t i=0;
    long atoms=0;
    std::map<std::string, int> ids;
    std::vector<std::string> functions= {
        "N", "Nf", "Pi", "I"
    };
    explicit Parser(const std::string&t):s(t) {
        for(int j=0; j<4; ++j)ids[functions[j]]=j;
    }
    [[noreturn]]void fail(std::string msg) {
        throw std::runtime_error(msg+" at "+std::to_string(i)+" near "+s.substr(i, 80));
    }
    void ws() {
        while(i<s.size()&&isspace((unsigned char)s[i]))++i;
    }
    bool eat(char c) {
        ws();
        if(i<s.size()&&s[i]==c) {
            ++i;
            return true;
        }
        return false;
    }
    Poly aspoly(const Node&n) {
        return n.special?n.poly:Poly {
            {
                Key {
                }
                , raw(n.a, n.b)
            }
        };
    }
    Poly multiply(const Poly&a, const Poly&b) {
        Poly out;
        for(auto&u:a)for(auto&v:b) {
            Key k=u.first;
            for(auto&e:v.first)k[e.first]+=e.second;
            Ptr c=times(u.second, v.second);
            normalize(k, c);
            out[k]=plus(out[k], c);
        }
        if(out.size()>50000)fail("Too many intermediate monomials");
        return out;
    }
    Node add(Node a, Node b, char op) {
        Node n {
            a.a, b.b, a.special||b.special, {
            }
        };
        if(n.special) {
            n.poly=a.special?std::move(a.poly):aspoly(a);
            for(auto&v:aspoly(b))n.poly[v.first]=plus(n.poly[v.first], op=='-'?minus(v.second):v.second);
        }
        return n;
    }
    Node mul(Node a, Node b, char op) {
        Node n {
            a.a, b.b, a.special||b.special, {
            }
        };
        if(n.special) {
            if(op=='*')n.poly=multiply(aspoly(a), aspoly(b));
            else {
                Poly bp=aspoly(b);
                if(bp.size()!=1)fail("Non-monomial denominator in transcendental/color/flavour variables");
                Key inv;
                for(auto&v:bp.begin()->first)inv[v.first]=-v.second;
                Ptr c=divide(unit(), bp.begin()->second);
                normalize(inv, c);
                n.poly=multiply(aspoly(a), Poly {
                    {
                        inv, c
                    }
                }
                );
            }
        }
        return n;
    }
    int integer(size_t a, size_t b) {
        std::string t;
        for(size_t j=a; j<b; ++j)if(!isspace((unsigned char)s[j])&&s[j]!='('&&s[j]!=')')t+=s[j];
        if(t.empty())fail("Empty exponent");
        size_t k=(t[0]=='-'||t[0]=='+')?1:0;
        for(size_t j=k; j<t.size(); ++j)if(!isdigit((unsigned char)t[j]))fail("Noninteger special exponent");
        int n=std::stoi(t);
        if(n < -8 || n>8)fail("Special exponent outside safety cap");
        return n;
    }
    Node expression() {
        auto a=term();
        while(true) {
            ws();
            if(i==s.size()||(s[i]!='+'&&s[i]!='-'))break;
            char op=s[i++];
            a=add(std::move(a), term(), op);
        }
        return a;
    }
    Node term() {
        auto a=unary();
        while(true) {
            ws();
            if(i==s.size()||(s[i]!='*'&&s[i]!='/'))break;
            char op=s[i++];
            a=mul(std::move(a), unary(), op);
        }
        return a;
    }
    Node unary() {
        ws();
        size_t start=i;
        if(eat('+')) {
            auto a=unary();
            a.a=start;
            return a;
        }
        if(eat('-')) {
            auto a=unary();
            if(a.special)for(auto&v:a.poly)v.second=minus(v.second);
            a.a=start;
            return a;
        }
        return power();
    }
    Node power() {
        auto a=atom();
        if(eat('^')) {
            auto b=unary();
            if(b.special)fail("Special-function exponent");
            Node n {
                a.a, b.b, a.special, {
                }
            };
            if(a.special) {
                int p=integer(b.a, b.b);
                if(p<0) {
                    if(a.poly.size()!=1)fail("Negative power of a sum of special monomials");
                    Key k;
                    for(auto&v:a.poly.begin()->first)k[v.first]=p*v.second;
                    Ptr c=ipow(a.poly.begin()->second, p);
                    normalize(k, c);
                    n.poly= {
                        {
                            k, c
                        }
                    };
                } else {
                    n.poly= {
                        {
                            Key {
                            }
                            , unit()
                        }
                    };
                    for(int j=0; j<p; ++j)n.poly=multiply(n.poly, a.poly);
                }
            }
            a=std::move(n);
        }
        return a;
    }
    Node atom() {
        ws();
        size_t start=i;
        ++atoms;
        if(atoms%1000000==0)std::cerr<<"Progress atoms "<<atoms<<" position "<<i<<"/"<<s.size()<<" primitives "<<functions.size()<<"\n";
        if(eat('(')) {
            auto a=expression();
            if(!eat(')'))fail("Expected )");
            a.a=start;
            a.b=i;
            return a;
        }
        if(i<s.size()&&isdigit((unsigned char)s[i])) {
            while(i<s.size()&&isdigit((unsigned char)s[i]))++i;
            return Node {
                start, i, false, {
                }
            };
        }
        if(i<s.size()&&isalpha((unsigned char)s[i])) {
            while(i<s.size()&&isalnum((unsigned char)s[i]))++i;
            std::string name=s.substr(start, i-start);
            ws();
            if(i<s.size()&&s[i]=='[') {
                if(name!="Log"&&name!="PolyLog"&&name!="Zeta"&&name!="Li22")fail("Unexpected function "+name);
                int depth=0;
                do {
                    char c=s[i++];
                    if(c=='[')++depth;
                    if(c==']')--depth;
                }
                while(i<s.size()&&depth);
                if(depth)fail("Unclosed function");
                name=s.substr(start, i-start);
                std::string t;
                for(char c:name)if(!isspace((unsigned char)c))t+=c;
                name=t;
                if(name.find("Nf")!=std::string::npos)fail("Nf in a special-function argument");
            } else if(name=="M2"||name=="muR2"||name=="x"||name=="z")return Node {
                start, i, false, {
                }
            };
            else if(name!="N"&&name!="Nf"&&name!="Pi"&&name!="I")fail("Unexpected symbol "+name);
            int id;
            if(ids.count(name))id=ids[name];
            else {
                id=functions.size();
                ids[name]=id;
                functions.push_back(name);
            }
            return Node {
                start, i, true, Poly {
                    {
                        Key {
                            {
                                id, 1
                            }
                        }
                        , unit()
                    }
                }
            };
        }
        fail("Expected atom");
    }
    Node parse() {
        Node n=expression();
        ws();
        if(i!=s.size())fail("Trailing input");
        return n;
    }
};
void render(std::ostream&o, const Ptr&p, const std::string&s) {
    if(!p) {
        o<<'0';
        return;
    }
    switch(p->kind) {
        case '1':o<<'1';
        break;
        case 'r':o<<'('<<s.substr(p->a, p->b-p->a)<<')';
        break;
        case '-':o<<"(-";
        render(o, p->args[0], s);
        o<<')';
        break;
        default:o<<'(';
        render(o, p->args[0], s);
        o<<p->kind;
        render(o, p->args[1], s);
        o<<')';
    }
}
int main(int argc, char**argv) {
    try {
        if(argc!=3) {
            std::cerr<<"Usage: collect_complete INPUT OUTPUT_DIR\n";
            return 2;
        }
        fs::create_directories(argv[2]);
        std::ifstream f(argv[1]);
        if(!f)throw std::runtime_error("Cannot read input");
        std::string s((std::istreambuf_iterator<char>(f)), {
        }
        ), clean;
        clean.reserve(s.size());
        for(size_t i=0; i<s.size(); ++i) {
            if(s[i]=='\\'&&i+1<s.size()&&s[i+1]=='\n') {
                ++i;
                continue;
            }
            if(s[i]=='\\'&&i+2<s.size()&&s[i+1]=='\r'&&s[i+2]=='\n') {
                i+=2;
                continue;
            }
            clean+=s[i];
        }
        s.clear();
        s.shrink_to_fit();
        auto eq=clean.find('=');
        if(eq!=std::string::npos) {
            if(clean.find('=', eq+1)!=std::string::npos)throw std::runtime_error("Multiple labels/assignments");
            clean.erase(0, eq+1);
        }
        Parser p(clean);
        auto node=p.parse();
        auto poly=p.aspoly(node);
        std::cerr<<"Parsed "<<p.atoms<<" atoms, "<<poly.size()<<" monomials, "<<p.functions.size()<<" primitives\n";
        std::ofstream fun(fs::path(argv[2])/"functions.txt");
        for(auto&a:p.functions)fun<<a<<'\n';
        std::ofstream ind(fs::path(argv[2])/"blocks.tsv");
        ind<<"block\tnf_power\tmonomial\n";
        int id=0;
        for(auto&v:poly) {
            int nf=v.first.count(1)?v.first.at(1):0;
            if(nf<0||nf>2)throw std::runtime_error("Nonpolynomial or unexpected flavour order");
            char name[64];
            snprintf(name, sizeof(name), "rational_%05d.wl", id);
            std::ofstream rat(fs::path(argv[2])/name);
            render(rat, v.second, clean);
            rat<<'\n';
            ind<<id<<'\t'<<nf<<'\t';
            bool first=true;
            for(auto&e:v.first) {
                if(e.first==1)continue;
                if(!first)ind<<',';
                ind<<e.first<<':'<<e.second;
                first=false;
            }
            ind<<'\n';
            ++id;
        }
        return 0;
    } catch(const std::exception&e) {
        std::cerr<<"ERROR: "<<e.what()<<'\n';
        return 1;
    }
}
