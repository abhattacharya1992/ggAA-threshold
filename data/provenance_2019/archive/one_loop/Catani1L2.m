(************************)

m[s,eps_] := m[M2,eps] * (Sum[ ((eps*(2*Log[1+x]-Log[x]))^i)/(Factorial[i]) ,{i,0,4} ] );
m[M2,eps_] := Sum[ ((eps*(Log[-M2]))^i)/(Factorial[i]) , {i,0,4} ];
p[muR,eps_] := Sum[ ((eps*(Log[muR2])/2)^i)/(Factorial[i]) , {i,0,4} ];

(* gamma[eps] = Gamma(1-eps) *)
gamma[eps_] := Exp[ Sum[ (zeta[i]/i * eps^i) , {i,2,4} ]  ];

zeta[2] = zeta2;
zeta[4] = zeta4;
zeta4 = 2/5*zeta2^2;
zeta2 = Pi^2/6;
zeta[3] = zeta3;
zeta3 = Zeta[3];
nf = Nf;
ma2 = M2;

(****BornCF = 1/(N^2-1)****)
Ca = N;
Cf = (N^2-1)/2/N;

beta0 = (11*Ca - 2*nf)/3;
beta1 = (34/3*Ca^2 - 2*nf*Cf - 10/3*Ca*nf);

I1[eps_] := (1/gamma[-eps/2])*( (- 4*Ca/eps^2+beta0/eps)*m[s, eps/2]*p[muR, -eps] );

NI1 =  Normal[Series[ I1[eps], {eps, 0, 2}]];

gg2AAR = Get["gg2AAR2"];

lo = Coefficient[gg2AAR,as,0];
nlo = Coefficient[gg2AAR,as,1];

cataniM1 = 2*NI1*lo;
oursM1 = nlo;


diff1 = Normal[Series[(oursM1-cataniM1)/.{Log[-M2]->Log[M2]-I*Pi},{eps,0,0}]];
diff10 = Coefficient[diff1,eps,0]>>Finite1L2;
diff1m1 = Coefficient[diff1,eps,-1];
diff1m2 = Coefficient[diff1,eps,-2];

