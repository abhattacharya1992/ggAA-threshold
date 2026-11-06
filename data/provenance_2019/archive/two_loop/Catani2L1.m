(************************)

m[s,eps_] := m[M2,eps] * (Sum[ ((eps*(2*Log[1+x]-Log[x]))^i)/(Factorial[i]) ,{i,0,4} ] );
m[M2,eps_] := Sum[ ((eps*(Log[-M2]))^i)/(Factorial[i]) , {i,0,4} ];
p[muR,eps_] := Sum[ ((eps*(Log[muR2])/2)^i)/(Factorial[i]) , {i,0,4} ];

(* gamma[eps] = Gamma(1-eps) *)
gamma[eps_] := Exp[ Sum[ (zeta[i]/i * eps^i) , {i,2,4} ]  ];

zeta[2] = zeta2;
zeta[3] = zeta3;
zeta[4] = zeta4;
zeta4 = 2/5*zeta2^2;
nf = Nf;
ma2 = M2;
m2 = M2;

(****BornCF = 1/(N^2-1)****)
Ca = N;
Cf = (N^2-1)/2/N;

beta0 = (11*Ca - 2*nf)/3;
beta1 = (34/3*Ca^2 - 2*nf*Cf - 10/3*Ca*nf);

H2g =1/2*( Ca^2 *( -5/12 - 11/24*zeta2 - 1/2*zeta[3])
     +Ca*nf*(29/27 + zeta2/12 )
     +Cf*nf*(-1/2)
     +nf^2*(-5/27) );

I1[eps_] := (1/gamma[-eps/2])*( (- 4*Ca/eps^2+beta0/eps)*m[s, eps/2]*p[muR, -eps] );

NI1 =  Normal[Series[ I1[eps], {eps, 0, 2}]];

I2 = Normal[Series[
     - (1/2)* I1[eps] *( I1[eps] + beta0 * (-2/eps) )
     + (gamma[-eps]/gamma[-eps/2]) * ( (- beta0/eps) +  ( (67/18 - zeta2)*Ca - 5*nf/9) ) * I1[2*eps]
     + 2*H2g/eps, {eps, 0, 0}]];

gg2AAR = Get["gg2AAR21"];

lo = Coefficient[gg2AAR,as,0];
nlo = Coefficient[gg2AAR,as,1];
nnlo = Coefficient[gg2AAR,as,2];

cataniM2 = 4*I2*lo + 2*NI1*nlo;
oursM2 = nnlo;

y=(1+x^2)/x-z;
x=1/11;
z=3/11;

diff2 = Normal[Series[(oursM2-cataniM2)/.{Log[-M2]->Log[M2]-I*Pi},{eps,0,0}]];
diff10  = Coefficient[diff2,eps,0]>>Finite2L1;
diff2m1 = Coefficient[diff2,eps,-1];
diff2m2 = Coefficient[diff2,eps,-2];
diff2m3 = Coefficient[diff2,eps,-3];
diff2m4 = Coefficient[diff2,eps,-4];
