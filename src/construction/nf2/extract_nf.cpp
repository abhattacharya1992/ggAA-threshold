// Extracts the Nf^2 terms of the two-loop finite remainders.
#include <algorithm>
#include <array>
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
struct E;
using P=std::shared_ptr<E>;
struct E {
    char kind;
    size_t a=0, b=0;
    std::vector<P> v;
    E(char k):kind(k) {
    }
};
P one() {
    static P x=std::make_shared<E>('1');
    return x;
}
P raw(size_t a, size_t b) {
    P x=std::make_shared<E>('r');
    x->a=a;
    x->b=b;
    return x;
}
P neg(P a) {
    if(!a)return {
    };
    if(a->kind=='n')return a->v[0];
    P x=std::make_shared<E>('n');
    x->v= {
        a
    };
    return x;
}
P add(P a, P b) {
    if(!a)return b;
    if(!b)return a;
    P x=std::make_shared<E>('+');
    x->v= {
        a, b
    };
    return x;
}
P mul(P a, P b) {
    if(!a||!b)return {
    };
    if(a->kind=='1')return b;
    if(b->kind=='1')return a;
    P x=std::make_shared<E>('*');
    x->v= {
        a, b
    };
    return x;
}
P divi(P a, P b) {
    if(!a)return {
    };
    if(!b)throw std::runtime_error("zero denominator");
    if(b->kind=='1')return a;
    P x=std::make_shared<E>('/');
    x->v= {
        a, b
    };
    return x;
}
using Pol=std::map<int, P>;
struct Node {
    size_t a, b;
    bool dep=false;
    Pol p;
};
class Parser {
public:
    const std::string &s;
    size_t i=0;
    long count=0;
    std::set<std::string> allowed= {
        "N", "Nf", "M2", "muR2", "x", "z", "I", "Pi", "Y", "Z"
    };
    Parser(const std::string&in):s(in) {
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
    [[noreturn]]void fail(std::string x) {
        throw std::runtime_error(x+" at byte "+std::to_string(i)+" near "+s.substr(i, 80));
    }
    Pol poly(const Node&x) {
        return x.dep?x.p:Pol {
            {
                0, raw(x.a, x.b)
            }
        };
    }
    Pol product(const Pol&a, const Pol&b) {
        Pol out;
        for(auto &u:a)for(auto&v:b) {
            if(u.first+v.first>8)fail("polynomial degree exceeds safety limit");
            out[u.first+v.first]=add(out[u.first+v.first], mul(u.second, v.second));
        }
        return out;
    }
    Node plus(Node a, Node b, char op) {
        Node out {
            a.a, b.b, a.dep||b.dep, {
            }
        };
        if(out.dep) {
            out.p=poly(a);
            for(auto&v:poly(b))out.p[v.first]=add(out.p[v.first], op=='-'?neg(v.second):v.second);
        }
        return out;
    }
    Node times(Node a, Node b, char op) {
        Node out {
            a.a, b.b, a.dep||b.dep, {
            }
        };
        if(out.dep) {
            if(op=='*')out.p=product(poly(a), poly(b));
            else {
                if(b.dep)fail("denominator depends on Nf");
                for(auto&v:poly(a))out.p[v.first]=divi(v.second, raw(b.a, b.b));
            }
        }
        return out;
    }
    Node expression() {
        auto a=term();
        while(true) {
            ws();
            if(i==s.size()||(s[i]!='+'&&s[i]!='-'))break;
            char op=s[i++];
            a=plus(a, term(), op);
        }
        return a;
    }
    Node term() {
        auto a=unary();
        while(true) {
            ws();
            if(i==s.size()||(s[i]!='*'&&s[i]!='/'))break;
            char op=s[i++];
            a=times(a, unary(), op);
        }
        return a;
    }
    Node unary() {
        ws();
        size_t st=i;
        if(eat('+')) {
            auto a=unary();
            a.a=st;
            return a;
        }
        if(eat('-')) {
            auto a=unary();
            if(a.dep)for(auto&v:a.p)v.second=neg(v.second);
            a.a=st;
            return a;
        }
        return power();
    }
    Node power() {
        auto a=atom();
        ws();
        if(eat('^')) {
            auto b=unary();
            if(b.dep)fail("Nf in exponent");
            Node out {
                a.a, b.b, a.dep, {
                }
            };
            if(a.dep) {
                std::string t;
                for(size_t j=b.a; j<b.b; ++j)if(!isspace((unsigned char)s[j])&&s[j]!='('&&s[j]!=')')t+=s[j];
                if(t.empty()||!std::all_of(t.begin(), t.end(), [](char c) {
                    return isdigit((unsigned char)c);
                }
                ))fail("noninteger/negative Nf-dependent power");
                int n=std::stoi(t);
                if(n>8)fail("large polynomial power");
                out.p= {
                    {
                        0, one()
                    }
                };
                for(int k=0; k<n; ++k)out.p=product(out.p, a.p);
            }
            a=out;
        }
        return a;
    }
    Node atom() {
        ws();
        size_t st=i;
        count++;
        if(eat('(')) {
            auto a=expression();
            if(!eat(')'))fail("expected )");
            a.a=st;
            a.b=i;
            return a;
        }
        if(i<s.size()&&isdigit((unsigned char)s[i])) {
            while(i<s.size()&&isdigit((unsigned char)s[i]))++i;
            return Node {
                st, i, false, {
                }
            };
        }
        if(i<s.size()&&isalpha((unsigned char)s[i])) {
            while(i<s.size()&&isalnum((unsigned char)s[i]))++i;
            std::string name=s.substr(st, i-st);
            ws();
            if(i<s.size()&&s[i]=='[') {
                if(name!="Log"&&name!="PolyLog"&&name!="Zeta"&&name!="Li22")fail("unknown function "+name);
                size_t fstart=i;
                int level=0;
                do {
                    char c=s[i++];
                    if(c=='[')++level;
                    if(c==']')--level;
                }
                while(i<s.size()&&level>0);
                if(level)fail("unclosed function");
                std::string args=s.substr(fstart, i-fstart);
                if(args.find("Nf")!=std::string::npos)fail("Nf-dependent function");
                return Node {
                    st, i, false, {
                    }
                };
            }
            if(!allowed.count(name))fail("unexpected symbol "+name);
            if(name=="Nf")return Node {
                st, i, true, {
                    {
                        1, one()
                    }
                }
            };
            return Node {
                st, i, false, {
                }
            };
        }
        fail("expected atom");
    }
    Node parse() {
        auto out=expression();
        ws();
        if(i!=s.size())fail("trailing input");
        return out;
    }
};
void printE(std::ostream&o, const P&p, const std::string&s) {
    if(!p) {
        o<<'0';
        return;
    }
    switch(p->kind) {
        case '1':o<<'1';
        break;
        case 'r':o<<'('<<s.substr(p->a, p->b-p->a)<<')';
        break;
        case 'n':o<<"(-";
        printE(o, p->v[0], s);
        o<<')';
        break;
        default:o<<'(';
        printE(o, p->v[0], s);
        o<<p->kind;
        printE(o, p->v[1], s);
        o<<')';
    }
}
int main(int argc, char**argv) {
    try {
        if(argc!=3) {
            std::cerr<<"Usage: extract_nf INPUT OUTPREFIX\n";
            return 2;
        }
        std::ifstream f(argv[1]);
        if(!f)throw std::runtime_error("input not readable");
        std::string s((std::istreambuf_iterator<char>(f)), {
        }
        );
        std::string clean;
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
            if(clean.find('=', eq+1)!=std::string::npos)throw std::runtime_error("multiple =");
            clean.erase(0, eq+1);
        }
        Parser p(clean);
        auto n=p.parse();
        auto pol=p.poly(n);
        std::cerr<<"atoms "<<p.count<<" degrees:";
        for(auto&v:pol)std::cerr<<' '<<v.first;
        std::cerr<<'\n';
        if(pol.rbegin()->first>2)throw std::runtime_error("Nf degree above two");
        for(int k=2; k<=2; ++k) {
            std::ofstream out(std::string(argv[2])+"_nf"+std::to_string(k)+".wl");
            printE(out, pol[k], clean);
            out<<'\n';
        }
        return 0;
    } catch(const std::exception&e) {
        std::cerr<<"ERROR: "<<e.what()<<'\n';
        return 1;
    }
}
