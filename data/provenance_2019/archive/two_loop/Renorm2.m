(**Overall operator renormalisation for the matrix element gg->AA**)

p[muR,eps_] := Sum[ ((eps*(Log[muR2])/2)^i)/(Factorial[i]) , {i,0,4} ];

Z2 = ( ascr*eps^-1*Ca * (
          + 44/3
          )

       + ascr*eps^-1*nf * (
          - 8/3
          )

       + ascr^2*eps^-2*Ca^2 * (
          + 484/9
          )

       + ascr^2*eps^-2*nf*Ca * (
          - 176/9
          )

       + ascr^2*eps^-2*nf^2 * (
          + 16/9
          )

       + ascr^2*eps^-1*Ca^2 * (
          + 136/3
          )

       + ascr^2*eps^-1*nf*Ca * (
          - 40/3
          )

       + ascr^2*eps^-1*nf*Cf * (
          - 8
          )

       + 1
         );

beta0 = (11*Ca - 2*nf)/3;
beta1 = (34/3*Ca^2 - 2*nf*Cf - 10/3*Ca*nf);

ascr = as*(1 + as*2*beta0/eps + as^2*(4*beta0^2/eps^2 + beta1/eps));
Ca = N;
Cf = (N^2-1)/2/N;

Z2 = Normal[Series[Expand[Z2],{as,0,2}]];
gg2AA0 = Get["ggAA02"];
gg2AA1 = Get["ggAA12"];
gg2AA2 = Get["ggAA22"];

y=(1+x^2)/x-z;

PT = (u*t - ma2^2)/s;
u = (2*ma2 - s -t);

gg2AA = ( gg2AA0 + as*p[muR,-eps]*gg2AA1 + as^2*p[muR,-2*eps]*(2*beta0/eps*p[muR,eps]*gg2AA1 + gg2AA2) );

gg2AAR = Normal[Series[Z2*gg2AA/gg2AA0,{as,0,2}]]>>gg2AAR22;

