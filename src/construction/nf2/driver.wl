(* Threshold expansion of the Nf^2 coefficients of P1 and P2. *)

Begin["AA2NF2`"];
$RecursionLimit=2048;
run=Environment["GGAA_NF2_RUN"];
If[!StringQ[run] || !DirectoryQ[run],Print["FATAL: no prepared run directory"];Exit[2]];
out=FileNameJoin[{run,"results"}]; If[!DirectoryQ[out],CreateDirectory[out]];
stages={}; warnings={}; stage="initializing";
statusText="Coefficient of Nf^2 of the two-loop finite remainders of P1 and P2.";
ass=0<beta<1 && -1<ct<1 && m2>0 && mu2>0 && Element[ell,Reals];
Clear[require,write,doStage,finiteQ,readExpression,numeric,pointRules,scaled];
require[ok_,msg_]:=If[!TrueQ[ok],Print["STOP: ",msg];Throw[Failure["NF2",<|"Stage"->stage,"Reason"->msg|>],"NF2Stop"]];
write[name_,val_]:=(Put[val,FileNameJoin[{out,name}]];val);
finiteQ[v_]:=NumberQ[v] && FreeQ[v,Indeterminate|_DirectedInfinity];
scaled[a_,b_]:=Max[Abs[Re[a-b]]/Max[1,Abs[Re[b]]],Abs[Im[a-b]]/Max[1,Abs[Im[b]]]];
SetAttributes[doStage,HoldRest];
doStage[name_,body_]:=Module[{t=AbsoluteTime[],ans},
 stage=name;Export[FileNameJoin[{run,"current_stage.txt"}],name,"Text"];
 Print[DateString[]," START ",name];
 ans=body;
 AppendTo[stages,<|"Stage"->name,"Seconds"->AbsoluteTime[]-t,"Status"->"COMPLETED"|>];
 write["stages.wl",stages];Print[DateString[]," COMPLETED ",name];ans];
(* Read the archived coefficient as exact input; no numerical literals are accepted here. *)
readExpression[name_String]:=Module[{txt,words,rules,h},
 txt=Import[FileNameJoin[{run,"payload","data",name}],"Text"];
 require[StringQ[txt],"Cannot read "<>name];
 words=DeleteDuplicates[StringCases[txt,RegularExpression["[A-Za-z][A-Za-z0-9]*"]]];
 require[Complement[words,{"M2","muR2","N","x","z","I","Pi","Log","beta","costh","L"}]==={},"Unexpected symbols in "<>name];
 require[StringDelete[txt,WhitespaceCharacter|DigitCharacter|LetterCharacter|"["|"]"|"("|")"|"+"|"-"|"*"|"/"|"^"|","]==="","Unexpected executable syntax in "<>name];
 rules=(RegularExpression["\\b"<>#[[1]]<>"\\b"]->"AA2NF2`"<>#[[2]])& /@
 {{"M2","m2"},{"muR2","mu2"},{"N","nc"},{"x","x"},{"z","z"},{"beta","beta"},{"costh","ct"},{"L","ell"}};
 h=Check[ToExpression[StringReplace[txt,rules],InputForm,HoldComplete],$Failed];
 require[Head[h]===HoldComplete,"Parser failed on "<>name];
 h=ReleaseHold[h];require[FreeQ[h,_Real],"Exact expression required in "<>name];h];
pointRules[bb_,cc_,mm_,uu_,col_,wp_]:=Module[{xx,zz},
 require[0<bb<1 && -1<cc<1 && mm>0 && uu>0 && col>0,"Point outside domain"];
 require[FreeQ[{bb,cc,mm,uu,col},_Real],"Use exact numerical inputs"];
 xx=(1-bb)/(1+bb);zz=(1+bb^2+2 bb cc)/(1-bb^2);
 Join[Thread[{x,z,m2,mu2,nc,beta,ct}->N[{xx,zz,mm,uu,col,bb,cc},wp]],
 {ell->N[Log[mm]-Log[uu],wp]}]];
numeric[e_,rules_,wp_]:=Module[{v},v=Block[{$MaxExtraPrecision=4 wp},N[e/.rules,wp]];require[finiteQ[v],"Nonfinite numerical value"];v];
(* The three stages below are kept separate so the symbolic and numerical checks can be inspected independently. *)
main[]:=Module[{raw,compact,pythonSeries,series4,series6,tests={},rows={},grids,pars,rr80,rr120,v80,v120,cv,sv4,sv6,i,p,k,ss,diff,jj,bs,cs,mm,uu,cc,bb,col,point,rf,ref,v,checks,results,plot,pe,table,fields,linearY,linearZ,explicit,constants,manifest},
 doStage["01_import_extracted_coefficients",
 raw=Table[readExpression["P"<>ToString[i]<>"_nf2.wl"],{i,2}];
 compact=Table[readExpression["P"<>ToString[i]<>"_nf2_factored.wl"],{i,2}];
 pythonSeries=Table[readExpression["P"<>ToString[i]<>"_nf2_series6.wl"],{i,2}];
 require[AllTrue[Join[raw,compact],FreeQ[#,Li22|PolyLog]&],"Unexpected special function in Nf^2 coefficient"];
 write["compact_nf2_formulas.wl",<|"P1"->compact[[1]],"P2"->compact[[2]],"Status"->statusText|>];
 manifest=<|"Kernel"->$Version,"System"->$SystemID,"Scope"->statusText,"P1SourceSHA256"->"0ad929948dad92936c5bf8a7ae00ef90b8956d0e4412e8fc74a68c77a8c315b1", "P2SourceSHA256"->"d7f76dbecc9c3aa4e41cf8c6d6b2d8191626c48d6e84a9e51e83ff5d342bd8e0"|>;
 write["manifest.wl",manifest];
 Print["Imported complete coefficients of Nf^2 for both projectors. No Y or Z definitions needed by these files."];
 ];
 doStage["02_independent_python_number_checks",
 rf=Import[FileNameJoin[{run,"payload","checks","python_nf2_checks.json"}],"RawJSON"];
 require[ListQ[rf]&&Length[rf]==8,"Missing Python benchmark records"];
 checks=Table[
 i=If[r["projector"]==="P1",1,2];bb=ToExpression[r["beta"]];cc=ToExpression[r["costh"]];mm=r["mass2"];uu=r["mu2"];
 require[StringMatchQ[r["real"],RegularExpression["-?[0-9]+\\.[0-9]+"]] && StringMatchQ[r["imag"],RegularExpression["-?[0-9]+\\.[0-9]+"]],"Unexpected numerical benchmark format"];
 ref=ToExpression[r["real"]<>"`70"]+I ToExpression[r["imag"]<>"`70"];
 v=numeric[compact[[i]],pointRules[bb,cc,mm,uu,3,100],100];
 <|"Projector"->i,"Beta"->bb,"CosTheta"->cc,"Difference"->scaled[v,ref],"Pass"->TrueQ[scaled[v,ref]<10^-40]|>,{r,rf}];
 write["python_number_checks.wl",checks];require[AllTrue[checks,TrueQ[#["Pass"]]&],"Cross-language numerical check failed"];
 Print["All 8 Python/Mathematica numerical comparisons passed."];
 ];
 doStage["03_threshold_series",
 linearY=(1+beta^2-2 beta ct)/(1-beta^2);linearZ=(1+beta^2+2 beta ct)/(1-beta^2);

 hfun[a_,b_]:=216 a^2+216 a b+72 b^2-(480+216 I Pi)a-(240+144 I Pi)b+200-69 Pi^2+240 I Pi;
 jfun[a_,b_,t_]:=36 t^2+(108 a+36 b-120-36 I Pi)t;
 blog=Log[4]-Log[1-beta^2];
 ry=m2*(linearY*linearZ-1)/linearY;rz=m2*(linearY*linearZ-1)/linearZ;
 explicit={64 I m2/(81(1-beta^2))*(hfun[ell,blog]+jfun[ell,blog,Log[linearY]]+jfun[ell,blog,Log[linearZ]]),
 16 I/81*((ry+rz)/2*hfun[ell,blog]+ry*jfun[ell,blog,Log[linearY]]+rz*jfun[ell,blog,Log[linearZ]])};
 makeSeries[e_,order_]:=Module[{w},w=TimeConstrained[Series[e,{beta,0,order},Assumptions->ass,Analytic->False],180,$Aborted];
 require[Head[w]===SeriesData && w[[5]]/w[[6]]>order,"Series unresolved or timed out"];
 w=MapAt[(FullSimplify[#,Assumptions->ass]& /@ #)&,w,3];w];
 series4=makeSeries[#,4]& /@ explicit;series6=makeSeries[#,6]& /@ explicit;
 Do[
 write["P"<>ToString[i]<>"_series4.wl",series4[[i]]];write["P"<>ToString[i]<>"_series6.wl",series6[[i]]];
 diff=TimeConstrained[FullSimplify[Normal[series6[[i]]]-pythonSeries[[i]],Assumptions->ass],120,$Aborted];
 AppendTo[tests,<|"Test"->"P"<>ToString[i]<>" Python symbolic coefficients","Pass"->TrueQ[diff===0],"Difference"->diff|>];
 diff=FullSimplify[Normal[Series[series6[[i]],{beta,0,4}]]-Normal[series4[[i]]],Assumptions->ass];
 AppendTo[tests,<|"Test"->"P"<>ToString[i]<>" fourth/sixth order agreement","Pass"->TrueQ[diff===0],"Difference"->diff|>];
 diff=FullSimplify[Normal[series6[[i]]]-(Normal[series6[[i]]]/.ct->-ct),Assumptions->ass];
 AppendTo[tests,<|"Test"->"P"<>ToString[i]<>" angular exchange","Pass"->TrueQ[diff===0],"Difference"->diff|>];
 AppendTo[tests,<|"Test"->"P"<>ToString[i]<>" odd powers vanish through six","Pass"->And@@(TrueQ[FullSimplify[SeriesCoefficient[series6[[i]],#],Assumptions->ass]===0]& /@ {1,3,5})|>];
 ,{i,2}];
 constants={SeriesCoefficient[series6[[1]],0],SeriesCoefficient[series6[[2]],0]};
 diff=FullSimplify[SeriesCoefficient[series6[[2]],2]-(1-ct^2)constants[[1]],Assumptions->ass];
 AppendTo[tests,<|"Test"->"P2 beta^2 coefficient equals (1-c^2) P1 constant","Pass"->TrueQ[diff===0],"Difference"->diff|>];
 write["symbolic_checks.wl",tests];require[AllTrue[tests,TrueQ[#["Pass"]]&],"Symbolic coefficient comparison failed"];
 write["threshold_leading_terms.wl",<|"P1Constant"->constants[[1]],"P2Constant"->constants[[2]],"P2Beta2"->SeriesCoefficient[series6[[2]],2],"Status"->statusText|>];
 write["coefficients_nf2.wl",Table[<|"Projector"->i,"Coefficients"->Table[<|"BetaPower"->k,"Coefficient"->SeriesCoefficient[series6[[i]],k]|>,{k,0,6}],"LogMassRatioSymbol"->ell,"Status"->statusText|>,{i,2}]];
 Print["All 9 symbolic checks passed. Computed orders 4 and 6 for the Nf^2 coefficient only."];
 ];
 doStage["04_unexpanded_coefficient_validation",
 pars={{15625,15625},{15625,62500},{22500,5625},{22500,22500}};
 bs={1/1000,1/100,1/20,1/10,1/5};cs={-2/3,0,1/3};
 Do[
 {mm,uu}=par;col=If[mm===22500,4,3];
 Print["Point: P",i," beta=",bb," c=",cc," M2=",mm," muR2=",uu," N=",col];
 rr80=pointRules[bb,cc,mm,uu,col,80];rr120=pointRules[bb,cc,mm,uu,col,120];
 v80=numeric[raw[[i]],rr80,80];v120=numeric[raw[[i]],rr120,120];cv=numeric[compact[[i]],rr120,120];
 sv4=numeric[Normal[series4[[i]]],rr120,120];sv6=numeric[Normal[series6[[i]]],rr120,120];
 require[scaled[v80,v120]<10^-30,"Unexpanded coefficient has not stabilized"];
 require[scaled[cv,v120]<10^-30,"Compact/unexpanded coefficient mismatch"];
 AppendTo[rows,<|"Projector"->i,"Beta"->bb,"CosTheta"->cc,"MassSquared"->mm,"MuSquared"->uu,"Nc"->col,
 "ReUnexpanded"->Re[v120],"ImUnexpanded"->Im[v120],"ReCompact"->Re[cv],"ImCompact"->Im[cv],
 "CompactScaledError"->scaled[cv,v120],"PrecisionChange80To120"->scaled[v80,v120],
 "P4AbsoluteError"->Abs[sv4-v120],"P6AbsoluteError"->Abs[sv6-v120],
 "P4RelativeError"->If[Abs[v120]>10^-50,Abs[(sv4-v120)/v120],Missing["NearZero"]],
 "P6RelativeError"->If[Abs[v120]>10^-50,Abs[(sv6-v120)/v120],Missing["NearZero"]]|>];
 write["validation_nf2.wl",<|"Rows"->rows,"Status"->statusText|>],
 {par,pars},{i,2},{cc,cs},{bb,bs}];
 fields=Keys[First[rows]];Export[FileNameJoin[{out,"validation_nf2.csv"}],Prepend[(Lookup[#,fields]/._Missing->"NA")& /@ rows,fields],"CSV"];
 Print[Length[rows]," point evaluations completed. Series errors are reported, not forced below a global tolerance."];
 ];
 doStage["05_comparison_tables_and_plots",
 Do[
 table=Table[jj=Select[rows,#["Projector"]===i && #["Beta"]===bb&];
 {bb,Max[Select[Lookup[jj,"P4RelativeError"],NumberQ]],Max[Select[Lookup[jj,"P6RelativeError"],NumberQ]]},{bb,bs}];
 write["P"<>ToString[i]<>"_error_table.wl",<|"Columns"->{"Beta","MaxSampledP4Error","MaxSampledP6Error"},"Data"->table,"Status"->"Maximum over the sampled angles and parameter sets."|>];
 If[Environment["GGAA_NF2_NO_PLOTS"]=!="1",
 pe=TimeConstrained[Check[
 plot=ListLogLogPlot[{Select[table[[All,{1,2}]],Last[#]>0&],Select[table[[All,{1,3}]],Last[#]>0&]},Joined->True,Frame->True,PlotRange->All,
 FrameLabel->{"beta_A","Maximum sampled relative error"},PlotLabel->"P"<>ToString[i]<>": two-loop Nf^2 coefficient only",PlotLegends->{"through beta_A^4","through beta_A^6"},PlotStyle->{Directive[Thick],Directive[Thick,Dashed]},ImageSize->640];
 Export[FileNameJoin[{out,"P"<>ToString[i]<>"_nf2_error.png"}],plot],$Failed],90,$Aborted];
 If[!StringQ[pe]||!FileExistsQ[pe],AppendTo[warnings,"P"<>ToString[i]<>" plot failed or timed out; numerical table saved."]]
 ];,{i,2}];
 write["complete_nf2_pilot.wl",<|"Compact"->compact,"Series4"->series4,"Series6"->series6,"SymbolicChecks"->tests,"Validation"->rows,"Warnings"->warnings,"Status"->statusText|>];
 ];True
];
result=CheckAbort[Catch[main[],"NF2Stop"],Failure["Aborted",<|"Stage"->stage|>]];
write["run_outcome.wl",<|"Outcome"->result,"Stages"->stages,"Warnings"->warnings,"Scope"->statusText|>];
Export[FileNameJoin[{run,"summary.txt"}],StringRiffle[Join[{
 "gg -> AA two-loop Nf^2 coefficient pilot",
 "Execution: "<>If[result===True,"COMPLETED_NF2_PILOT_ONLY","FAILED_OR_INTERRUPTED"],
 "Kernel: "<>$Version,
 "Scope: "<>statusText,
 "One-loop expansions were not recomputed.",
 "The remaining coefficients of Nf^0 and Nf^1 were NOT calculated in this pilot.",
 "Li22 and other higher-weight functions still require a verified method for those remaining contributions.",
 "The series are compared with the exact Nf^2 coefficients of the unexpanded expressions.",
 "See extraction_status.txt for whether native source re-extraction was repeated on this server.",
 "Truncation errors are listed for the validation grid.",
 "STAGES:"},(#["Stage"]<>" : "<>#["Status"])& /@ stages,{"WARNINGS:"},If[warnings==={}, {"None recorded"},warnings]],"\n"],"Text"];
If[result===True,Export[FileNameJoin[{run,"COMPLETED_NF2_PILOT_ONLY.txt"}],statusText,"Text"]];
Print["RESULT: ",If[result===True,"COMPLETED_NF2_PILOT_ONLY",result]];
Exit[If[result===True,0,1]];
