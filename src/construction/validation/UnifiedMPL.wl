(* Threshold expansions of logarithms, polylogarithms and Li22 from recurrences in beta. *)

Begin["AAUnified`"];
$HistoryLength=0;
$RecursionLimit=16384;
$MaxExtraPrecision=1000;
backendVersion="GGAA-Unified-RealMPL22-0.1";
assumptions=0<beta<1 && -1<costh<1 && m2>0 && mu2>0 && nc>1;
kinMap={x->(1-beta)/(1+beta), z->(1+beta^2+2 beta costh)/(1-beta^2)};
run=Environment["GGAA_UNIFIED_RUN"];
root=Environment["GGAA_UNIFIED_PROJECT"];
require[ok_,message_String]:=If[!TrueQ[ok],Print["STOP: ",message];Throw[Failure["GGAAUnified",<|"Reason"->message,"Stage"->currentStage|>],"UnifiedStop"]];
finiteQ[v_]:=NumberQ[v] && FreeQ[v,Indeterminate|_DirectedInfinity];
currentStage="initialize"; warnings={}; stages={};
say[s_String]:=WriteString[First[$Output],"[",DateString[{"ISODate","T","Time"}],"] ",s,"\n"];
makeDir[path_String]:=If[!DirectoryQ[path],CreateDirectory[path,CreateIntermediateDirectories->True]];
save[path_String,value_]:=Module[{tmp=path<>".tmp",ok},makeDir[DirectoryName[path]];ok=Check[Put[value,tmp];True,False];require[ok && FileExistsQ[tmp],"Cannot write "<>path];If[FileExistsQ[path],DeleteFile[path]];RenameFile[tmp,path];path];
load[path_String]:=Module[{v},require[FileExistsQ[path],"Missing saved file: "<>path];v=Check[Get[path],$Failed];require[v=!=$Failed,"Cannot read "<>path];v];
SetAttributes[stage,HoldRest];
stage[name_String,body_]:=Module[{t=AbsoluteTime[],v,msgs},currentStage=name;Export[FileNameJoin[{run,"current_stage.txt"}],name<>"\n","Text"];say["START "<>name];{v,msgs}=Block[{$MessageList={}},v=body;{v,$MessageList}];If[Length[msgs]>0,AppendTo[warnings,name<>": Wolfram messages "<>ToString[msgs,InputForm]]];require[v=!=$Failed && v=!=$Aborted && !FailureQ[v],"Stage returned failure: "<>name];AppendTo[stages,<|"Stage"->name,"Status"->"COMPLETED","Seconds"->AbsoluteTime[]-t,"Messages"->msgs|>];save[FileNameJoin[{run,"results","stages.wl"}],stages];say["COMPLETED "<>name];v];
warn[s_String]:=(AppendTo[warnings,s];say["WARNING "<>s]);

(* Restricted parser for the exact rational/MPL expressions written by the preprocessing step. *)
safeParse[text_String]:=Module[{words,residue,rules,h,e},
 words=DeleteDuplicates[StringCases[text,RegularExpression["[A-Za-z][A-Za-z0-9]*"]]];
 require[Complement[words,{"M2","muR2","N","Nf","x","z","I","Pi","Log","PolyLog","Zeta","Li22"}] === {},"Unexpected source token(s): "<>ToString[Complement[words,{"M2","muR2","N","Nf","x","z","I","Pi","Log","PolyLog","Zeta","Li22"}],InputForm]];
 residue=StringDelete[text,WhitespaceCharacter|DigitCharacter|LetterCharacter|"["|"]"|"("|")"|"+"|"-"|"*"|"/"|"^"|","];
 require[residue==="","Unexpected executable or nonexact input syntax."];
 rules=(RegularExpression["\\b"<>#[[1]]<>"\\b"]->"AAUnified`"<>#[[2]])&/@{{"M2","m2"},{"muR2","mu2"},{"N","nc"},{"Nf","nf"},{"x","x"},{"z","z"},{"Li22","MPL22"}};
 h=Check[ToExpression[StringReplace[text,rules],InputForm,HoldComplete],$Failed];require[Head[h]===HoldComplete,"Wolfram parser rejected a source block."];
 e=ReleaseHold[h];require[FreeQ[e,_Real],"Approximate source numbers are not allowed."];e
];
canon[q_]:=Cancel[Together[q]];
zeroQ[q_]:=If[FailureQ[q]||q===$Aborted||q===$Failed,False,TrueQ[canon[q]===0]];

coefficientCanonical[q_]:=Module[{atoms,syms,e,back,v},
 atoms=DeleteDuplicates[Cases[q,_Log|_PolyLog|_MB,Infinity]];
 syms=Array[ConstTag,Length[atoms]];back=Thread[syms->atoms];e=q/.Thread[atoms->syms];
 v=If[Length[syms]==0,canon[e],Collect[Expand[e],syms,canon]];
 v/.back
];

val[q0_]:=Module[{q=canon[q0],n,d},If[q===0,Return[Infinity]];{n,d}={Numerator[q],Denominator[q]};require[PolynomialQ[n,beta]&&PolynomialQ[d,beta],"Expected rational beta coefficient."];Exponent[n,beta,Min]-Exponent[d,beta,Min]];
limit0[q0_]:=Module[{q=canon[q0],v,n,d,nv,dv},If[q===0,Return[0]];v=val[q];require[v>=0,"Primitive has an unsupported infinite beta endpoint."];If[v>0,Return[0]];n=Numerator[q];d=Denominator[q];nv=Exponent[n,beta,Min];dv=Exponent[d,beta,Min];canon[Coefficient[n,beta,nv]/Coefficient[d,beta,dv]]];
rationalJet[q0_,hi_Integer]:=Module[{q=canon[q0],n,d,nv,dv,v,nn,dd,a,k,j},
 If[q===0,Return[<||>]];n=Numerator[q];d=Denominator[q];require[PolynomialQ[n,beta]&&PolynomialQ[d,beta],"Nonrational beta input."];
 nv=Exponent[n,beta,Min];dv=Exponent[d,beta,Min];v=nv-dv;require[IntegerQ[v],"Fractional beta valuation is unsupported."];If[v>hi,Return[<||>]];
 nn=Table[Coefficient[n,beta,nv+k],{k,0,hi-v}];dd=Table[Coefficient[d,beta,dv+k],{k,0,hi-v}];a=ConstantArray[0,hi-v+1];
 Do[a[[k+1]]=canon[(nn[[k+1]]-Sum[dd[[j+1]]a[[k-j+1]],{j,1,k}])/dd[[1]]],{k,0,hi-v}];
 AssociationThread[Range[v,hi],a]
];
jetAdd[a_Association,d_Association]:=Merge[{a,d},Total];
jetScale[a_Association,k_]:=Map[k #&,a];
jetMul[a_Association,d_Association,hi_Integer]:=Module[{out=<||>,k,j},Do[If[k+j<=hi,AssociateTo[out,(k+j)->(Lookup[out,k+j,0]+a[k]d[j])]],{k,Keys[a]},{j,Keys[d]}];Map[Expand,out]];
jetPower[a_Association,n_Integer,hi_Integer]:=Module[{out=<|0->1|>,j},require[n>=0,"A negative power of a special function is unsupported."];Do[out=jetMul[out,a,hi],{j,n}];out];
jetExpression[j_Association]:=Total[KeyValueMap[#2 beta^#1&,j]];
jetFromExpression[e_,hi_Integer]:=Module[{ex=Expand[e/.Log[beta]->Lb],lo},lo=Exponent[ex,beta,Min];require[IntegerQ[lo],"Input is not a Laurent/log polynomial."];AssociationThread[Range[lo,hi],Table[Coefficient[ex,beta,k],{k,lo,hi}]]];

(* Real-log conventions for the physical threshold region. *)
positive[q_]:=TrueQ[FullSimplify[q>0,Assumptions->(-1<costh<1 && m2>0 && mu2>0)]];
realLog[q0_]:=Module[{q=canon[q0],fac,total=0,f,p,sgn},
 require[positive[q],"Cannot prove positivity of a leading logarithm: "<>ToString[q,InputForm]];
 fac=FactorList[q];
 Do[{f,p}=u;If[f===1,Continue[]];Which[positive[f],total+=p Log[f],positive[-f],total+=p Log[-f],True,Return[Log[q]]],{u,fac}];
 total
];
logJet[q0_,hi_Integer]:=logJet[q0,hi]=Module[{q=canon[q0],v,lead,der,j},
 require[q=!=0,"Logarithm of an identically zero expression."];v=val[q];require[IntegerQ[v],"Noninteger logarithm valuation."];
 lead=limit0[q/beta^v];j=<|0->(v Lb+realLog[lead])|>;
 der=rationalJet[canon[D[q,beta]/q-v/beta],hi-1];
 KeyValueMap[(require[#1>=0 || zeroQ[#2],"Unaccounted negative power in d log."];If[#1>=0,AssociateTo[j,(#1+1)->(#2/(#1+1))]])&,der];j
];

integrateJet[d_Association,seed_,hi_Integer]:=Module[{out=<|0->seed|>,v,co,n,k,h,t,difference},
 difference=coefficientCanonical[Lookup[d,-1,0]-D[seed,Lb]];
 require[zeroQ[difference],"The 1/beta coefficient does not match the endpoint logarithms. Residual: "<>ToString[difference,InputForm]];
 Do[If[n< -1,require[zeroQ[coefficientCanonical[d[n]]],"Unexpected stronger-than-log primitive divergence."]];If[n<0 || n+1>hi,Continue[]];
 v=Expand[d[n]];require[PolynomialQ[v,Lb],"Primitive derivative is not polynomial in Log beta."];co=CoefficientList[v,Lb];t=0;
 Do[t+=co[[k+1]] Sum[(-1)^h Factorial[k]/Factorial[k-h] Lb^(k-h)/(n+1)^(h+1),{h,0,k}],{k,0,Length[co]-1}];AssociateTo[out,(n+1)->coefficientCanonical[t]],{n,Keys[d]}];out
];

polyConstant[n_Integer,q_]:=Which[q===0,0,q===1,Zeta[n],q===-1,-(1-2^(1-n))Zeta[n],True,PolyLog[n,q]];
polyJet[n_Integer,q0_,hi_Integer]:=polyJet[n,q0,hi]=Module[{q=canon[q0],q0lim,der,v,low,d,seed},
 If[q===0,Return[<||>]];
 If[n===0,Return[rationalJet[canon[q/(1-q)],hi]]];
 If[n===1,Return[jetScale[logJet[canon[1-q],hi],-1]]];
 require[2<=n<=4,"Only classical weights through four are implemented."];
 q0lim=limit0[q];require[TrueQ[FullSimplify[q0lim<=1,Assumptions->(-1<costh<1)]],"Classical endpoint outside supported real branch."];
 seed=polyConstant[n,q0lim];der=canon[D[q,beta]/q];If[der===0,Return[<|0->seed|>]];
 v=val[der];low=polyJet[n-1,q,hi-1-v];d=jetMul[rationalJet[der,hi-1],low,hi-1];integrateJet[d,seed,hi]
];

(* Boundary constants for the depth-two polylogarithms used in this amplitude. *)
boundary[p_Integer,q_Integer,a_,w_]:=Which[
 a===0||w===0,0,
 p===2&&q===2&&a===1&&w===1,Pi^4/120,
 p===2&&q===2&&a===-1&&w===1,-Pi^4/480,
 p===2&&q===1&&a===1&&w===1,Zeta[3],
 p===1&&q===1&&zeroQ[a-w],Log[1-a]^2/2,
 p===1&&q===1&&w===1,-PolyLog[2,-a/(1-a)],
 True,MB[p,q,a,w]
];
mpSeed[p_Integer,q_Integer,a_,w_]:=Module[{a0=limit0[a],w0=limit0[w],la,lw,va,vw,rho,c0},
 If[a0===0||w0===0,Return[0]];
 If[a0=!=1||p>1,Return[boundary[p,q,a0,w0]]];
 la=logJet[canon[1-a],0][0];
 If[q===2,c0=If[w0===1,-2 Zeta[3],-boundary[2,1,w0,w0]-PolyLog[3,w0]];Return[-PolyLog[2,w0]la+c0]];
 require[q===1,"Unsupported multiple-polylog endpoint weight."];
 If[w0=!=1,lw=Log[1-w0];Return[lw la-PolyLog[2,w0]-lw^2/2]];
 va=val[canon[1-a]];vw=val[canon[1-w]];
 Which[vw===va,rho=limit0[(1-w)/(1-a)];la^2/2+PolyLog[2,1-rho],
 vw>va,la^2/2+Pi^2/6,
 True,lw=logJet[canon[1-w],0][0];la^2/2-(lw-la)^2/2-Pi^2/6]
];

rationalTimes[rat_,provider_,hi_Integer]:=Module[{r=canon[rat],v,f,lo},If[r===0,Return[<||>]];v=val[r];f=provider[hi-v];If[Length[f]===0,Return[<||>]];lo=Min[Keys[f]];jetMul[rationalJet[r,hi-lo],f,hi]];
mpJet[p_Integer,q_Integer,a0_,w0_,hi_Integer]:=mpJet[p,q,a0,w0,hi]=Module[{a=canon[a0],w=canon[w0],A,B,sd,d1,d2},
 If[a===0||w===0,Return[<||>]];
 If[p===0,Return[rationalTimes[a/(1-a),Function[h,polyJet[q,w,h]],hi]]];
 If[q===0,
   If[zeroQ[a-w],Return[jetAdd[polyJet[p-1,a,hi],jetScale[polyJet[p,a,hi],-1]]]];
   Return[jetAdd[rationalTimes[w/(a-w),Function[h,polyJet[p,a,h]],hi],rationalTimes[-a/(a-w),Function[h,polyJet[p,w,h]],hi]]]
 ];
 require[1<=p<=2&&1<=q<=2,"Unsupported multiple-polylog indices."];
 sd=mpSeed[p,q,a,w];A=canon[D[a,beta]/a];B=canon[D[w,beta]/w-A];
 d1=rationalTimes[A,Function[h,mpJet[p-1,q,a,w,h]],hi-1];
 d2=rationalTimes[B,Function[h,mpJet[p,q-1,a,w,h]],hi-1];
 integrateJet[jetAdd[d1,d2],sd,hi]
];

mpNumber[p_Integer,q_Integer,ae_,we_,wp_Integer]:=mpNumber[p,q,ae,we,wp]=Module[{a,w,delta,pts,f,v,t,goal},
 require[FreeQ[{ae,we},_Real]||Min[Precision/@{ae,we}]>=wp,"Insufficient primitive input precision."];
 require[TrueQ[-1<=ae<=1&&-1<=we<=1],"MPL numerical point outside the implemented real domain [-1,1]."];
 If[ae===0||we===0,Return[N[0,wp]]];require[!(ae===1&&p===1),"Divergent p=1 boundary supplied to finite integral."];
 If[ae===1&&we===1&&p===2&&q===2,Return[N[Pi^4/120,wp]]];
 If[ae===1&&we===1&&p===2&&q===1,Return[N[Zeta[3],wp]]];
 a=N[ae,wp+30];w=N[we,wp+30];delta=Min[Select[{1-a,1-w},#>0&]/.{}->{1}];
 pts={0,1/2,9/10,1};If[0<delta<1/100,pts=Join[pts,{1-Sqrt[delta],1-delta}]];pts=Sort[DeleteDuplicates[pts]];
 f[t_?NumericQ]:=If[t==0||t==1,0,a(-Log[t])^(p-1)If[q==1,-Log[1-w t],PolyLog[q,w t]]/((1-a t)Factorial[p-1])];
 goal=Max[20,wp-15];
 v=TimeConstrained[Check[NIntegrate[f[t],Evaluate[Prepend[pts,t]],WorkingPrecision->wp+25,PrecisionGoal->goal,AccuracyGoal->goal,MaxRecursion->40,Method->"DoubleExponential"],$Failed],600,$Aborted];
 require[v=!=$Failed&&v=!=$Aborted&&finiteQ[v],"Numerical multiple polylog failed, indices/arguments: "<>ToString[{p,q,ae,we,wp},InputForm]];
 require[Precision[v]>=Max[15,wp-25],"Multiple-polylog quadrature returned insufficient precision."];N[v,wp]
];
MB /: N[MB[p_Integer,q_Integer,a_,w_],wp_?NumberQ]:=mpNumber[p,q,a,w,Ceiling[wp]];
finiteExpressionNumber[e_,wp_Integer]:=Module[{atoms,rules,v},
 atoms=DeleteDuplicates[Cases[e,_MB,Infinity]];
 rules=(#->mpNumber[#[[1]],#[[2]],#[[3]],#[[4]],wp+10])&/@atoms;
 v=Block[{$MaxExtraPrecision=5wp},N[e/.rules,wp]];require[finiteQ[v],"Expression retains unresolved numerical functions."];v
];
scaled[a_,d_]:=Max[Abs[Re[a-d]]/Max[1,Abs[Re[d]]],Abs[Im[a-d]]/Max[1,Abs[Im[d]]]];

functionSpec[e_]:=Which[
 e===nc,<|"Kind"->"Constant","Value"->nc|>,e===nf,<|"Kind"->"Flavour","Value"->nf|>,
 e===Pi||e===I,<|"Kind"->"Constant","Value"->e|>,
 Head[e]===Zeta,<|"Kind"->"Constant","Value"->e|>,
 Head[e]===Log,If[FreeQ[e,beta|x|z],<|"Kind"->"Constant","Value"->e|>,<|"Kind"->"Log","Argument"->canon[e[[1]]/.kinMap]|>],
 Head[e]===PolyLog,<|"Kind"->"PolyLog","Weight"->e[[1]],"Argument"->canon[e[[2]]/.kinMap]|>,
 Head[e]===MPL22,<|"Kind"->"MPL22","A"->canon[e[[1]]/.kinMap],"W"->canon[e[[1]]e[[2]]/.kinMap]|>,
 True,require[False,"Unknown primitive: "<>ToString[e,InputForm]]
];
primitiveJet[sp_Association,hi_Integer]:=Switch[sp["Kind"],
 "Constant",<|0->sp["Value"]|>,"Log",logJet[sp["Argument"],hi],
 "PolyLog",polyJet[sp["Weight"],sp["Argument"],hi],
 "MPL22",mpJet[2,2,sp["A"],sp["W"],hi],_,require[False,"Unexpected primitive kind."]];

primitiveNumber[sp_Association,bb_,cc_,mm_,uu_,nn_,wp_Integer]:=Module[{key,path,record,q,a,w,v,subs},
 subs={beta->bb,costh->cc,m2->mm,mu2->uu,nc->nn};
 If[sp["Kind"]==="Constant",v=Block[{$MaxExtraPrecision=5wp},N[sp["Value"]/.subs,wp]];require[finiteQ[v],"Nonfinite constant."];Return[v]];
 key=If[sp["Kind"]==="MPL22",{"MPL22",canon[sp["A"]/.subs],canon[sp["W"]/.subs],wp}, {sp["Kind"],Lookup[sp,"Weight",0],canon[sp["Argument"]/.subs],wp}];
 path=FileNameJoin[{run,"checkpoints","numbers",IntegerString[Hash[key,"SHA256"],16,64]<>".wl"}];
 If[FileExistsQ[path],record=load[path];require[record["Key"]===key&&record["Version"]===backendVersion,"Stale primitive cache."];Return[record["Value"]]];
 Switch[sp["Kind"],"Log",q=key[[3]];require[TrueQ[q>0],"Log numeric argument not positive."];v=N[Log[q],wp],
 "PolyLog",q=key[[3]];require[TrueQ[q<=1],"Classical argument crosses the implemented real branch."];v=N[PolyLog[sp["Weight"],q],wp],
 "MPL22",v=mpNumber[2,2,key[[2]],key[[3]],wp],_,require[False,"Unknown numeric primitive."]];
 require[finiteQ[v],"Nonfinite primitive evaluation."];save[path,<|"Key"->key,"Value"->v,"Version"->backendVersion|>];v
];

engineTests[]:=Module[{tests={},a,w,j,expr,d,ref,btest=1/1000,ctest=1/3,valnum,err,seedcases,seen,sp,path,v},
 AppendTo[tests,<|"Test"->"common denominator zero","Pass"->zeroQ[1/(costh-1)-costh/(costh-1)+1]|>];
 AppendTo[tests,<|"Test"->"nonzero not erased","Pass"->!zeroQ[1/(costh-1)-costh/(costh-1)+2]|>];
 j=mpJet[2,2,beta,beta^2,6];expr=jetExpression[j]/.Lb->Log[beta];
 ref=Sum[If[k+2 n<=6,beta^(k+2 n)/((k+n)^2 n^2),0],{k,1,6},{n,1,3}];
 AppendTo[tests,<|"Test"->"Li22 exact convergent-series coefficients","Pass"->zeroQ[expr-ref]|>];
 j=mpJet[1,1,beta,beta,6];ref=Normal[Series[Log[1-beta]^2/2,{beta,0,6}]];
 AppendTo[tests,<|"Test"->"Li11 exact repeated argument","Pass"->zeroQ[jetExpression[j]-ref]|>];

 Do[a=pair[[1]];w=pair[[2]];ref=Sum[a^k w^n/((k+n)^2 n^2),{k,1,140},{n,1,140}];v=mpNumber[2,2,a,w,55];err=Abs[v-N[ref,55]];
 AppendTo[tests,<|"Test"->"integral versus convergent nested sum","Arguments"->{a,w},"Error"->err,"Pass"->TrueQ[err<10^-35]|>],{pair,{{1/5,3/10},{-2/5,1/3},{1/2,1/4}}}];
 v=mpNumber[2,2,-1,1,55];AppendTo[tests,<|"Test"->"Li22 alternating stuffle constant","Pass"->TrueQ[Abs[v+N[Pi^4/480,55]]<10^-35]|>];
 save[FileNameJoin[{run,"results","engine_regression_tests.wl"}],tests];require[And@@Lookup[tests,"Pass"],"A multiple-polylog engine test failed."];tests
];
End[];
