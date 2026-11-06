(**Overall operator renormalisation for the matrix element gg->AA**)

p[muR,eps_] := Sum[ ((eps*(Log[muR2])/2)^i)/(Factorial[i]) , {i,0,4} ];

ZGG = ascr^2*((484 Ca^2)/(9 eps^2)+(34 Ca^2)/(3 eps)-(176 Ca nf)/(9 eps^2)-(10 Ca nf)/(3 eps)-(2 Cf nf)/eps+(16 nf^2)/(9 eps^2)) + ascr*((22 Ca)/(3 eps)-(4 nf)/(3 eps))+1;

ZGJ = ascr*(-24*Cf/eps) + ascr^2*(((1/eps^2)*(-176*Ca*Cf+32*Cf*nf))+((1/eps)*((-284*Ca*Cf/3)+(84*Cf^2)+(8*Cf*nf/3))));

beta0 = (11*Ca - 2*nf)/3;
beta1 = (34/3*Ca^2 - 2*nf*Cf - 10/3*Ca*nf);

ascr = as;
Zas = (1 + as*2*beta0/eps + as^2*(4*beta0^2/eps^2 + beta1/eps));
asunr = as * p[muR,-eps] * Zas;

Ca = N;
Cf = (N^2-1)/2/N;

ZGG = Normal[Series[Expand[ZGG],{as,0,2}]];
ZGJ = Normal[Series[Expand[ZGJ],{as,0,2}]];

gg2AA0  = Get["ggAA01"];
gg2AA1  = Get["ggAA11"];
gg2AA2  = Get["ggAA21"];
ggjAA1  = Get["ggjAA11"];
ggjAA2  = Get["ggjAA21"];
ggjjAA1  = Get["ggjjAA11"];
ggjjAA2  = Get["ggjjAA21"];

PT = (u*t - ma2^2)/s;
u = (2*ma2 - s - t);
s = ma2*(1+x)^2/x;
t = -ma2*y;
ma2 = M2;
CJ = -as*Cf*(3/2-3*Log[muR2/mt^2])*CG;
CG=1;


gg2AA = ( gg2AA0 + asunr*gg2AA1 + asunr^2*gg2AA2 )*CG^2;
gg2AAj = ( asunr*ggjAA1 + asunr^2*ggjAA2 )*CG;
gg2AAjj = (asunr*ggjjAA1 + asunr^2*ggjjAA2 )*CJ^2;

gg2AAR = Normal[Series[(((ZGG^2)*gg2AA)+(2*ZGG*ZGJ*gg2AAj)+((ZGJ^2)*gg2AAjj)),{as,0,2}]]>>gg2AAR21;
