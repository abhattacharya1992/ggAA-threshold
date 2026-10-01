(* Threshold expansion of the Nf coefficient of P1 through beta^4 and beta^6. *)

Begin["AA2NF1Threshold`"];
$HistoryLength = 0;
$RecursionLimit = 4096;
run = Environment["GGAA_NF1_THRESHOLD_RUN"];
payload = DirectoryName[$InputFileName];
out = FileNameJoin[{run, "results"}];
checkdir = FileNameJoin[{run, "checkpoints"}];
stages = {}; warnings = {}; current = "startup"; validationRows = {};
scope = "Coefficient of Nf in P1.";
assumptions = 0 < beta < 1 && -1 < costh < 1 && m2 > 0 && mu2 > 0 && nc > 1;
kinematicMap = {x -> (1-beta)/(1+beta), z -> (1+beta^2+2 beta costh)/(1-beta^2)};
logBack = {lbeta -> Log[beta], lminus -> Log[1-costh], lplus -> Log[1+costh],
  lmass -> Log[m2], lscale -> Log[mu2], ltwo -> Log[2], zthree -> Zeta[3], piconst -> Pi};
formalVariables = {lbeta,lminus,lplus,lmass,lscale,ltwo,zthree,piconst,nc};

say[text_String] := WriteString[First[$Output], "[", DateString[{"ISODate","T","Time"}], "] ", text, "\n"];
write[name_String, value_] := Module[{path = FileNameJoin[{out,name}],ok}, ok=Check[Put[value,path];True,False]; require[TrueQ[ok] && FileExistsQ[path],"Cannot save "<>path]; path];
require[test_, text_String] := If[!TrueQ[test], say["STOP: "<>text]; Throw[Failure["ThresholdCheck",<|"Reason"->text,"Stage"->current|>],"StopThreshold"]];
finiteQ[v_] := NumberQ[v] && FreeQ[v,Indeterminate|_DirectedInfinity];
scaled[a_,d_] := Max[Abs[Re[a-d]]/Max[1,Abs[Re[d]]], Abs[Im[a-d]]/Max[1,Abs[Im[d]]]];
SetAttributes[stage,HoldRest];
stage[name_String,body_] := Module[{t=AbsoluteTime[],v},
 current=name; Export[FileNameJoin[{run,"current_stage.txt"}],name<>"\n","Text"]; say["START "<>name];
 v=body; require[v =!= $Failed && v =!= $Aborted && !FailureQ[v],"Stage returned an unresolved failure."];
 AppendTo[stages,<|"Stage"->name,"Status"->"COMPLETED","Seconds"->(AbsoluteTime[]-t)|>];
 write["stages.wl",stages]; say["COMPLETED "<>name]; v
];

(* Parse only the small Mathematica subset used by the archived coefficient files. *)
parseText[text_String] := Module[{words,residue,rules,h,e},
 words=DeleteDuplicates[StringCases[text,RegularExpression["[A-Za-z][A-Za-z0-9]*"]]];
 require[Complement[words,{"M2","muR2","N","x","z","I","Pi","Log","PolyLog","Zeta"}] === {},"Unexpected token in coefficient input."];
 residue=StringDelete[text,WhitespaceCharacter|DigitCharacter|LetterCharacter|"["|"]"|"("|")"|"+"|"-"|"*"|"/"|"^"|","];
 require[residue === "","Unexpected arithmetic syntax in coefficient input."];
 rules=(RegularExpression["\\b"<>#[[1]]<>"\\b"] -> "AA2NF1Threshold`"<>#[[2]])& /@
   {{"M2","m2"},{"muR2","mu2"},{"N","nc"},{"x","x"},{"z","z"}};
 h=Check[ToExpression[StringReplace[text,rules],InputForm,HoldComplete],$Failed];
 require[Head[h] === HoldComplete,"Input could not be parsed."]; e=ReleaseHold[h];
 require[FreeQ[e,_Real],"Input must contain only exact arithmetic."]; e
];

valuation[q0_] := Module[{q=Cancel[Together[q0]],n,d},
 If[q === 0,Return[Infinity]]; n=Numerator[q];d=Denominator[q];
 require[PolynomialQ[n,beta] && PolynomialQ[d,beta],"Expected a rational function of beta."];
 Exponent[n,beta,Min]-Exponent[d,beta,Min]
];
rationalJet[q0_,hi_Integer] := Module[{q=Cancel[Together[q0]],n,d,nv,dv,v,aa,nn,dd,deg,k,j},
 If[q === 0,Return[<||>]];n=Numerator[q];d=Denominator[q];
 require[PolynomialQ[n,beta] && PolynomialQ[d,beta],"Non-rational primitive encountered."];
 nv=Exponent[n,beta,Min];dv=Exponent[d,beta,Min];v=nv-dv;
 require[IntegerQ[v],"Fractional beta powers require a different implementation."];
 If[v>hi,Return[<||>]];deg=hi-v;
 nn=Table[Coefficient[n,beta,nv+k],{k,0,deg}];
 dd=Table[Coefficient[d,beta,dv+k],{k,0,deg}];
 require[dd[[1]] =!= 0,"Zero leading denominator."];
 aa=ConstantArray[0,deg+1];
 Do[aa[[k+1]]=Cancel[(nn[[k+1]]-Total[Table[dd[[j+1]] aa[[k-j+1]],{j,1,k}]])/dd[[1]]],{k,0,deg}];
 AssociationThread[Range[v,hi],aa]
];
jetMul[a_Association,d_Association,hi_Integer] := Module[{r=<||>,k,j},
 Do[If[k+j<=hi,AssociateTo[r,(k+j)->(Lookup[r,k+j,0]+a[k] d[j])]],{k,Keys[a]},{j,Keys[d]}];
 Map[Expand,r]
];
jetPower[a_Association,n_Integer,hi_Integer] := Module[{r=<|0->1|>,j},
 require[n>=0,"Negative powers of special-function atoms are not supported."];
 Do[r=jetMul[r,a,hi],{j,n}];r
];

logConstant[q0_] := Module[{q=Cancel[q0],choices,hit},
 choices={{1,0},{2,ltwo},{1/2,-ltwo},{m2,lmass},{mu2,lscale},
 {1-costh,lminus},{1+costh,lplus},{2(1-costh),ltwo+lminus},
 {2(1+costh),ltwo+lplus},{1-costh^2,lminus+lplus},{2(1-costh^2),ltwo+lminus+lplus}};
 hit=Select[choices,TrueQ[Cancel[q-#[[1]]]===0]&];
 require[Length[hit]>=1,"Unresolved leading logarithm or branch: "<>ToString[q,InputForm]];
 hit[[1,2]]
];
logJet[q0_,hi_Integer] := logJet[q0,hi] = Module[{q=Cancel[q0],v,g0,der,r,k},
 require[q =!= 0,"Logarithm of an identically zero argument."];
 v=valuation[q];g0=Cancel[q/beta^v] /. beta->0;
 r=<|0->(v lbeta+logConstant[g0])|>;
 der=rationalJet[Cancel[D[q,beta]/q-v/beta],hi-1];
 Do[If[Cancel[der[k]]=!=0,
   require[k>=0,"A logarithm derivative retained an unaccounted pole."];
   AssociateTo[r,(k+1)->(der[k]/(k+1))]],{k,Keys[der]}];Map[Expand,r]
];
polyConstant[n_Integer,q0_] := Which[
 q0===0,0,
 n===2 && q0===1,piconst^2/6,
 n===2 && q0===-1,-piconst^2/12,
 n===2 && q0===1/2,piconst^2/12-ltwo^2/2,
 n===3 && q0===1,zthree,
 n===3 && q0===-1,-3 zthree/4,
 n===3 && q0===1/2,7 zthree/8-piconst^2 ltwo/12+ltwo^3/6,
 True,require[False,"Unresolved polylogarithm boundary constant."]
];
jetIntegral[a_Association] := Module[{r=<||>,v,cl,k,h,t,n},
 Do[v=Expand[a[n]];If[v===0,Continue[]];
   require[n>=0,"Finite polylogarithm derivative retained a negative-power term."];
   require[PolynomialQ[v,lbeta],"Non-polynomial logarithmic coefficient."];
   cl=CoefficientList[v,lbeta];t=0;
   Do[t=t+cl[[k+1]] Sum[(-1)^h Factorial[k]/Factorial[k-h] lbeta^(k-h)/(n+1)^(h+1),{h,0,k}],{k,0,Length[cl]-1}];
   AssociateTo[r,(n+1)->Expand[t]],{n,Keys[a]}];r
];
polyJet[n_Integer,q0_,hi_Integer] := polyJet[n,q0,hi] = Module[{q=Cancel[q0],lim,low,der,r},
 If[n===1,Return[Map[-#&,logJet[Cancel[1-q],hi]]]];
 require[MemberQ[{2,3},n],"Only classical weights 2 and 3 are implemented in this driver."];
 lim=q /. beta->0;
 require[MemberQ[{0,1,-1,1/2},lim],"Unresolved or divergent polylogarithm endpoint."];
 If[FreeQ[q,beta],Return[<|0->polyConstant[n,lim]|>]];
 low=polyJet[n-1,q,hi];der=rationalJet[Cancel[D[q,beta]/q],hi-1];
 r=jetIntegral[jetMul[der,low,hi-1]];AssociateTo[r,0->polyConstant[n,lim]];r
];
primitive[a_,hi_Integer] := Which[
 Head[a]===Log,logJet[Cancel[a[[1]] /. kinematicMap],hi],
 Head[a]===PolyLog,polyJet[a[[1]],Cancel[a[[2]] /. kinematicMap],hi],
 a===Zeta[3],<|0->zthree|>,
 True,require[False,"Unknown special-function atom: "<>ToString[a,InputForm]]
];
jetExpression[j_Association] := Total[KeyValueMap[#2 beta^#1&,j]];

reduceCoefficient[a_] := Module[{pol,cr},
 pol=Expand[a];If[pol===0,Return[0]];
 require[PolynomialQ[pol,formalVariables],"Coefficient is not polynomial in the formal constants."];
 cr=CoefficientRules[pol,formalVariables];
 Total[(Cancel[Last[#]] Times@@MapThread[Power,{formalVariables,First[#]}])& /@ cr]
];

(* Rebuild the block sum once before expanding; this catches a bad archive or kinematic map early. *)
loadAndVerify[] := Module[{path,txt,lines,rebuild,plain,hash,syms,sprules,d,j,k,bl,pref,q,v,checks={},endpoints={}},
 path=FileNameJoin[{run,"input","P1_nf1_compact.wl"}];
 hash=IntegerString[FileHash[path,"SHA256"],16,64];
 require[hash==="af362f3a8757a4335570347ccbbd781f87d060df87d163ab123bf37ac21958f4","Compact input hash mismatch."];
 txt=Import[path,"Text"];compactExpression=parseText[txt];
 rawpath=FileNameJoin[{run,"input","P1_nf1.wl"}];
 require[IntegerString[FileHash[rawpath,"SHA256"],16,64] === "0da89f4e7a22e433478b3e810b13cd720cf6c5eb29e168f771d1acd200be3815","Raw coefficient hash mismatch."];
 rawExpression=parseText[Import[rawpath,"Text"]];
 data=Get[FileNameJoin[{payload,"blocks.wl"}]];atoms=data["Atoms"];rats=data["RationalFunctionsBeta"];blocks=data["Blocks"];
 require[Length[atoms]===63 && Length[rats]===156 && Length[blocks]===526,"Unexpected decomposition dimensions."];
 lines=Select[StringTrim /@ StringSplit[txt,{"\r\n","\n"}],StringLength[#]>0&];
 require[Length[lines]===Length[blocks],"The per-line block proof cannot match this input."];
 syms=Table[Unique["special"],{Length[atoms]}];sprules=Dispatch[Thread[atoms->syms]];
 Do[bl=blocks[[j]];plain=If[StringEndsQ[lines[[j]],"+"],StringDrop[lines[[j]],-1],lines[[j]]];
  rebuild=bl[[1]] data["RationalFunctionsXZ"][[bl[[2]]]] Times@@((atoms[[#[[1]]]]^#[[2]])& /@ bl[[3]]);
  d=TimeConstrained[Cancel[(parseText[plain]/.sprules)-(rebuild/.sprules)],120,$Aborted];
  require[d===0,"Block decomposition proof failed or timed out at block "<>ToString[j]];
  If[Mod[j,100]===0,say["Proved compact block "<>ToString[j]<>"/526"]],{j,Length[blocks]}];
 Do[q=TimeConstrained[Cancel[Together[(data["RationalFunctionsXZ"][[j]] /. kinematicMap)-rats[[j]]]],120,$Aborted];
  require[q===0,"Rational kinematic map check failed at "<>ToString[j]];
  v=valuation[rats[[j]]];require[v===data["Valuations"][[j]],"Rational valuation mismatch."];
  AppendTo[checks,<|"RationalID"->j,"MapPass"->True,"Valuation"->v|>],{j,Length[rats]}];
 minVal=Min[Lookup[checks,"Valuation"]]; require[minVal===-1,"Unexpected denominator depth; do not reuse prepared depth blindly."];
 Do[If[Head[atoms[[j]]]===PolyLog,
   q=Cancel[atoms[[j,2]] /. kinematicMap];
   AppendTo[endpoints,<|"Atom"->atoms[[j]],"Limit"->(q/.beta->0),
     "OneMinusArgumentValuation"->valuation[1-q]|>]],{j,Length[atoms]}];
 write["input_and_depth_checks.wl",<|"CompactBlocksProved"->526,"RationalChecks"->checks,"MinimumValuation"->minVal,"Scope"->scope|>];
 write["endpoint_inventory.wl",endpoints];
 say["526 block identities and 156 rational maps verified. Minimum beta valuation = -1."];
 True
];

makeExpansion[p_Integer] := Module[{hi=p-minVal,aj={},rj={},pieces={},ans=<||>,bl,mon,part,pref,k,j,co,rec,sd,expr,idx,termfile},
 say["Computing primitive coefficients through beta^"<>ToString[hi]<>" for amplitude order "<>ToString[p]];
 Do[co=TimeConstrained[primitive[atoms[[j]],hi],180,$Aborted];require[AssociationQ[co],"Primitive expansion failed at atom "<>ToString[j]];
  AppendTo[aj,co];If[Mod[j,10]===0,say["Primitive "<>ToString[j]<>"/63 at order "<>ToString[p]]],{j,Length[atoms]}];
 Put[<|"Order"->p,"PrimitiveDepth"->hi,"Atoms"->atoms,"Jets"->aj|>,FileNameJoin[{checkdir,"atoms_p"<>ToString[p]<>".wl"}]];
 Do[co=TimeConstrained[rationalJet[rats[[j]],p],180,$Aborted];require[AssociationQ[co],"Rational expansion failed at "<>ToString[j]];AppendTo[rj,co],{j,Length[rats]}];
 Do[bl=blocks[[j]];pref=Cancel[bl[[1]] nc/m2] /. Pi->piconst;
  require[FreeQ[pref,m2|mu2|beta|costh|x|z],"Unexpected mass or kinematic prefactor."];
  mon=<|0->1|>;
  Do[mon=jetMul[mon,jetPower[aj[[idx[[1]]]],idx[[2]],hi],hi],{idx,bl[[3]]}];
  part=jetMul[rj[[bl[[2]]]],mon,p];
  Do[AssociateTo[ans,k->(Lookup[ans,k,0]+pref part[k])],{k,Keys[part]}];
  If[Mod[j,50]===0 || j===Length[blocks],
   Put[<|"Order"->p,"LastBlock"->j,"NormalizedPartialCoefficients"->ans|>,FileNameJoin[{checkdir,"assembly_p"<>ToString[p]<>".wl"}]];
   say["Assembled block "<>ToString[j]<>"/526 at order "<>ToString[p]]],{j,Length[blocks]}];
 Do[co=TimeConstrained[reduceCoefficient[ans[k]],300,$Aborted];require[co=!=$Aborted,"Coefficient simplification timed out at beta^"<>ToString[k]];
  AssociateTo[ans,k->co],{k,Keys[ans]}];
 ans=KeySelect[ans, !TrueQ[ans[#]===0]&];normalizedSeries[p]=ans;
 expr=(m2/nc) jetExpression[ans];
 co=Expand[expr /. logBack];
 sd=SeriesData[beta,0,Table[Expand[(m2/nc) Lookup[ans,k,0] /. logBack],{k,Min[Keys[ans]],p}],Min[Keys[ans]],p+1,1];
 rec=<|"Order"->p,"Angle"->Automatic,"Series"->sd,"Expression"->co,"FormalExpression"->expr,
   "NormalizedCoefficients"->ans,"PrimitiveDepth"->hi,"Scope"->scope,"Remainder"->"Terms beyond the requested power, with possible logarithms. No convergence radius is asserted."|>;
 write["P1_nf1_series"<>ToString[p]<>".wl",rec];
 say["Saved P1 nf1 threshold series through beta^"<>ToString[p]];rec
];

(* Exact coefficient-level checks against the independently generated Python recurrence. *)
checkAlgebra[] := Module[{refs,rows={},actual,expected,d,p,k,symm,odd,neg,logdeg,coeffrows},
 refs=Get[FileNameJoin[{payload,"python_reference_series.wl"}]];
 Do[actual=normalizedSeries[p];expected=refs["NormalizedOrder"<>ToString[p]];
  Do[d=reduceCoefficient[Lookup[actual,k,0]-Coefficient[expected,beta,k]];
   AppendTo[rows,<|"Test"->"Order "<>ToString[p]<>" beta^"<>ToString[k]<>" vs Python","Pass"->TrueQ[d===0],"Difference"->d|>];
   require[d===0,"Exact comparison to independent Python coefficients failed."],{k,0,p}];
  neg=Select[Keys[actual],#<0&];
  require[neg==={},"Inverse power of beta in the series."];
  Do[symm=Expand[actual[k] - (actual[k] /. {costh->-costh,lminus->lplus,lplus->lminus})];
   d=reduceCoefficient[symm];AppendTo[rows,<|"Test"->"Order "<>ToString[p]<>" beta^"<>ToString[k]<>" angular exchange","Pass"->TrueQ[d===0],"Difference"->d|>];
   require[d===0,"Angular exchange check failed."],{k,Keys[actual]}],{p,{4,6}}];
 Do[d=reduceCoefficient[Lookup[normalizedSeries[4],k,0]-Lookup[normalizedSeries[6],k,0]];
  AppendTo[rows,<|"Test"->"Fourth/sixth common beta^"<>ToString[k],"Pass"->TrueQ[d===0],"Difference"->d|>];
  require[d===0,"Common coefficients disagree between orders."],{k,0,4}];
 write["symbolic_checks.wl",rows];
 coeffrows=Flatten[Table[
   (<|"BetaPower"->k,"LogBetaPower"->First[First[#]],"Coefficient"->Expand[(m2/nc) Last[#] /. logBack]|>& /@
      CoefficientRules[normalizedSeries[6][k],{lbeta}]),{k,Sort[Keys[normalizedSeries[6]]]}],1];
 write["coefficients_P1_nf1.wl",<|"Order"->6,"Rows"->coeffrows,"Scope"->scope|>];
 Export[FileNameJoin[{out,"coefficients_P1_nf1.txt"}],ToString[coeffrows,InputForm],"Text"];
 write["threshold_structure.wl",<|"Powers"->Keys[normalizedSeries[6]],
  "LogPowers"->Union[Flatten[Table[Flatten[First /@ CoefficientRules[normalizedSeries[6][k],{lbeta}]],{k,Keys[normalizedSeries[6]]}]]],
  "LeadingFormalCoefficient"->(m2/nc) normalizedSeries[6][Min[Keys[normalizedSeries[6]]]],"Scope"->scope|>];
 say["Exact coefficient, angular and common-order checks passed."];True
];

pointRules[bb_,cc_,mm_,uu_,nn_,wp_Integer] := Thread[
 {x,z,m2,mu2,nc,beta,costh} -> N[{(1-bb)/(1+bb),(1+bb^2+2 bb cc)/(1-bb^2),mm,uu,nn,bb,cc},wp]];
evaluate[e_,rr_,wp_Integer] := Module[{v},v=Block[{$MaxExtraPrecision=5 wp},N[e/.rr,wp]];
 require[finiteQ[v],"A numerical expression did not evaluate to a finite value."];v];
(* Numerical grid used to judge the size of the beta^4 and beta^6 truncation errors. *)
checkNumbers[] := Module[{cases,bs,cs,par,bb,cc,r100,r140,a100,a140,compact,s4,s6,e4,e6,row,cols},
 cases={{15625,15625,3},{15625,62500,3},{22500,5625,4}};
 bs={1/1000,1/500,1/200,1/100,1/50,1/20,1/10,1/5,3/10};cs={-2/3,-1/3,0,1/3,2/3};
 Do[say["Numerical check "<>ToString[{par,bb,cc},InputForm]];
  r100=pointRules[bb,cc,par[[1]],par[[2]],par[[3]],100];r140=pointRules[bb,cc,par[[1]],par[[2]],par[[3]],140];
  a100=evaluate[rawExpression,r100,100];a140=evaluate[rawExpression,r140,140];
  require[scaled[a100,a140]<10^-30,"Raw coefficient did not stabilize with precision."];
  compact=evaluate[compactExpression,r140,140];require[scaled[compact,a140]<10^-35,"Raw and compact coefficients disagree."];
  s4=evaluate[series4["Expression"],r140,140];s6=evaluate[series6["Expression"],r140,140];
  e4=Abs[s4-a140];e6=Abs[s6-a140];
  row=<|"MassSquared"->par[[1]],"MuSquared"->par[[2]],"Nc"->par[[3]],"Beta"->bb,"CosTheta"->cc,
    "ReUnexpanded"->Re[a140],"ImUnexpanded"->Im[a140],"ReOrder4"->Re[s4],"ImOrder4"->Im[s4],"ReOrder6"->Re[s6],"ImOrder6"->Im[s6],
    "AbsoluteError4"->e4,"AbsoluteError6"->e6,"AbsoluteRealError4"->Abs[Re[s4-a140]],"AbsoluteImagError4"->Abs[Im[s4-a140]],
    "AbsoluteRealError6"->Abs[Re[s6-a140]],"AbsoluteImagError6"->Abs[Im[s6-a140]],
    "RelativeError4"->If[Abs[a140]>10^-60,e4/Abs[a140],Missing["NearZeroReference"]],
    "RelativeError6"->If[Abs[a140]>10^-60,e6/Abs[a140],Missing["NearZeroReference"]],
    "RawPrecisionChange"->scaled[a100,a140],"RawCompactDifference"->scaled[compact,a140]|>;
  AppendTo[validationRows,row];write["validation_P1_nf1.wl",<|"Rows"->validationRows,"Scope"->scope|>];
  If[bb===1/1000,require[e6/Max[1,Abs[a140]]<10^-12,"Threshold sanity check failed at beta=1/1000."]],
  {par,cases},{cc,cs},{bb,bs}];
 cols=Keys[First[validationRows]];
 Export[FileNameJoin[{out,"validation_P1_nf1.csv"}],Prepend[(Lookup[#,cols]/._Missing->"NA")& /@ validationRows,cols],"CSV"];
 say["All 135 numerical points evaluated."];True
];

makePlots[] := Module[{rows,grp,tab,g4,g6,png,plot,r,maxRel},
 maxRel[v_List] := If[AllTrue[v,NumberQ],Max[v],Missing["UndefinedRelativeError"]];
 rows=Select[validationRows,# ["MassSquared"]===15625 && #["MuSquared"]===15625 && #["Nc"]===3&];
 grp=SortBy[GatherBy[rows,#["Beta"]&],First[#]["Beta"]&];
 tab=Table[<|"Beta"->First[r]["Beta"],"MaxSampledRelativeError4"->maxRel[Lookup[r,"RelativeError4"]],
  "MaxSampledRelativeError6"->maxRel[Lookup[r,"RelativeError6"]]|>,{r,grp}];
 write["P1_nf1_error_table.wl",<|"Data"->tab,"Angles"->{-2/3,-1/3,0,1/3,2/3},"Scope"->"Maximum over the sampled angles, m2=mu2=15625, Nc=3."|>];
 If[Environment["GGAA_NF1_NO_PLOTS"]==="1",say["Plots disabled; error table is saved."];Return[True]];
 g4={#["Beta"],#["MaxSampledRelativeError4"]}& /@ tab;g6={#["Beta"],#["MaxSampledRelativeError6"]}& /@ tab;
 png=CheckAbort[TimeConstrained[
   plot=ListLogLogPlot[{Select[g4,NumberQ[Last[#]]&&Last[#]>0&],Select[g6,NumberQ[Last[#]]&&Last[#]>0&]},Joined->True,Frame->True,
    FrameLabel->{"beta_A","Maximum sampled relative error"},PlotLabel->"P1: two-loop coefficient linear in Nf only",PlotLegends->{"through beta_A^4","through beta_A^6"},
    PlotStyle->{Directive[Thick],Directive[Thick,Dashed]},PlotRange->All,ImageSize->800];
   Export[FileNameJoin[{out,"P1_nf1_error.png"}],plot],300,$Aborted],$Aborted];
 If[!StringQ[png] || !FileExistsQ[png],AppendTo[warnings,"PNG export failed or timed out. The numerical error table is preserved."];say[Last[warnings]]];True
];

main[] := (
 require[StringQ[run] && DirectoryQ[out] && DirectoryQ[checkdir],"Run directory is not initialized by the launcher."];
 say["Mathematica "<>$Version];say[scope];
 stage["01_input_decomposition_and_depth",loadAndVerify[]];
 stage["02_fourth_order_series",series4=makeExpansion[4]];
 stage["03_sixth_order_series",series6=makeExpansion[6]];
 stage["04_symbolic_validation",checkAlgebra[]];
 stage["05_numerical_validation",checkNumbers[]];
 stage["06_error_tables_and_plot",makePlots[]];
 stage["07_combined_results",write["complete_P1_nf1_threshold.wl",<|"Series4"->series4,"Series6"->series6,"Validation"->validationRows,"Stages"->stages,"Warnings"->warnings,"Scope"->scope,"Kernel"->$Version|>]];
 True
);
outcome=CheckAbort[Catch[main[],"StopThreshold"],Failure["Aborted",<|"Stage"->current|>]];
If[outcome =!= True,say["FAILED: "<>ToString[outcome,InputForm]]];
write["run_outcome.wl",<|"Outcome"->outcome,"Stages"->stages,"Warnings"->warnings,"Scope"->scope|>];
Export[FileNameJoin[{run,"summary.txt"}],StringRiffle[{
 "gg -> AA: P1 two-loop coefficient linear in Nf, threshold expansion",
 "Execution: "<>If[outcome===True,"COMPLETED_P1_NF1_THRESHOLD_ONLY","FAILED_OR_INTERRUPTED"],
 "Kernel: "<>$Version,"Scope: "<>scope,
 "Earlier results were not modified.",
 "The series is compared with the unexpanded expression.",
 "Computed target orders: beta_A^4 and beta_A^6; angle, mass, scale and colour remain symbolic.",
 "Evaluated numerical points: "<>ToString[Length[validationRows]],
 "STAGES:",StringRiffle[(#["Stage"]<>" : "<>#["Status"])& /@ stages,"\n"],
 "WARNINGS:",If[warnings==={},"None recorded.",StringRiffle[warnings,"\n"]],
 "Remaining: coefficient of Nf in P2 and the Nf^0 coefficients.",
 "See results/complete_P1_nf1_threshold.wl, validation_P1_nf1.csv, symbolic_checks.wl and full_session.log."
},"\n"]<>"\n","Text"];
If[outcome===True,Export[FileNameJoin[{run,"COMPLETED_P1_NF1_THRESHOLD_ONLY.txt"}],scope<>"\n","Text"]];
Export[FileNameJoin[{run,"current_stage.txt"}],If[outcome===True,"COMPLETED_P1_NF1_THRESHOLD_ONLY", "FAILED: "<>current]<>"\n","Text"];
say["Final outcome: "<>ToString[outcome,InputForm]];
Exit[If[outcome===True,0,1]];
