(* Definitions of the construction used to evaluate the unexpanded expressions. *)

Get[FileNameJoin[{DirectoryName[$InputFileName],"UnifiedMPL.wl"}]];
Begin["AAUnified`"];
componentLabels={"P2_nf1","P1_nf0","P2_nf0"};
knownLabels={"P1_nf1","P1_nf2","P2_nf2"};
scope="Coefficients P2 nf^1, P1 nf^0 and P2 nf^0, assembled with P1 nf^1 and the nf^2 coefficients.";
allBlocks=<||>;functions=<||>;computed=<||>;known=<||>;validationRows={};primitiveChecks={};algebraRows={};
makeDir[FileNameJoin[{run,"results"}]];makeDir[FileNameJoin[{run,"checkpoints"}]];

hashFile[p_]:=IntegerString[FileHash[p,"SHA256"],16,64];
(* Previously completed nf pieces are read back as exact SeriesData. *)
readKnown[name_String]:=Module[{e,syms,rules,allowed},
 e=Block[{$Context="AAUnified`",$ContextPath={"AAUnified`","System`"}},load[FileNameJoin[{run,"known",name}]]];
 If[AssociationQ[e],require[KeyExistsQ[e,"Series"],"Known record lacks Series: "<>name];e=e["Series"]];
 require[Head[e]===SeriesData,"Expected saved SeriesData in "<>name];
 syms=DeleteDuplicates[Cases[e,s_Symbol/;Context[s]=!="System`",Infinity]];
 allowed={"beta","ct","costh","m2","mass2","ell","mu2","nc","nf"};
 require[And@@(MemberQ[allowed,SymbolName[#]]&/@syms),"Unknown symbol in saved previous result: "<>name];
 rules=Map[Function[s,s->Switch[SymbolName[s],"beta",beta,"ct"|"costh",costh,"m2"|"mass2",m2,"ell",Log[m2]-Log[mu2],"mu2",mu2,"nc",nc,"nf",nf]],syms];
 e/.rules
];

homDegree[e_]:=Module[{vals},Which[e===0,Missing["Zero"],e===m2,1,FreeQ[e,m2],0,
 Head[e]===Times,vals=homDegree/@(List@@e);If[MemberQ[vals,_Missing],Missing["Zero"],Total[vals]],
 Head[e]===Power&&IntegerQ[e[[2]]],e[[2]]homDegree[e[[1]]],
 Head[e]===Plus,vals=DeleteDuplicates[DeleteCases[homDegree/@(List@@e),_Missing]];require[Length[vals]<=1,"Rational source has mixed mass degrees; do not normalize by guesswork."];If[Length[vals]==0,Missing["Zero"],First[vals]],
 True,require[False,"Unexpected mass dependence in rational coefficient."]]];

nativeParse[text_String] := Module[{words,residue,h,e},
 words=DeleteDuplicates[StringCases[text,RegularExpression["[A-Za-z][A-Za-z0-9]*"]]];
 require[Complement[words,{"costh","Order","LowerBound","Jet"}]==={},"Unexpected token in native rational output."];
 residue=StringDelete[text,WhitespaceCharacter|DigitCharacter|LetterCharacter|"\""|"<"|">"|"|"|"("|")"|"+"|"-"|"*"|"/"|"^"|","];
 require[residue==="","Unexpected native rational syntax."];
 h=Check[ToExpression[StringReplace[text,RegularExpression["\\bcosth\\b"]->"AAUnified`costh"],InputForm,HoldComplete],$Failed];
 require[Head[h]===HoldComplete,"Cannot parse native rational output."];
 e=ReleaseHold[h];require[FreeQ[e,_Real|Indeterminate|_DirectedInfinity],"Native rational output is not exact."];e
];

nativeExactZero[e_,label_String] := Module[{v},
 v=TimeConstrained[canon[e],180,$Aborted];
 If[v=!=0,save[FileNameJoin[{run,"results","native_rational_check_failure.wl"}],<|"Label"->label,"Difference"->v|>]];
 require[v===0,"Exact native rational coefficient test failed or timed out: "<>label];True
];

(* Reconstruct the rational input record used by the archived calculation. *)
loadAndReduce[p_String] := Module[{base,lines,specs,index,rows={},fields,id,k,pairs,src,cp,rec,old,j=0,
 native,audit,j4,j6,key,path,jhash,refs,rf,bb,rv,difference,sourcehash},
 base=FileNameJoin[{run,"collected",p}];
 lines=Import[FileNameJoin[{base,"functions.txt"}],"Lines"];
 require[ListQ[lines]&&Length[lines]>=4,"Cannot read primitive inventory."];
 specs=functionSpec[safeParse[#]]&/@lines;AssociateTo[functions,p->specs];
 index=Rest[Import[FileNameJoin[{base,"blocks.tsv"}],"Lines"]];
 Do[j++;fields=StringSplit[line,"\t"];require[Length[fields]>=2,"Malformed block index."];
 id=ToExpression[fields[[1]]];k=ToExpression[fields[[2]]];
 require[IntegerQ[id]&&IntegerQ[k]&&0<=k<=2,"Invalid block id or flavour order."];
 pairs=If[Length[fields]<3,{},(ToExpression/@StringSplit[# ,":"])&/@StringSplit[fields[[3]],","]];
 require[And@@(Length[#]==2&&IntegerQ[#[[1]]]&&IntegerQ[#[[2]]]&/@pairs),"Invalid monomial."];
 src=FileNameJoin[{base,"rational_"<>IntegerString[id,10,5]<>".wl"}];sourcehash=hashFile[src];
 old=FileNameJoin[{run,"checkpoints","reduced",p,"block_"<>IntegerString[id,10,5]<>".wl"}];
 If[FileExistsQ[old],
   rec=load[old];
   require[rec["Version"]===backendVersion&&rec["Monomial"]===pairs&&rec["NfPower"]===k&&rec["SourceBlockHash"]===sourcehash,"Rational checkpoint does not match its source block."];
   If[rec["RationalBeta"]=!=0,AppendTo[rows,Append[rec,"RationalMethod"->"ExistingExactClosedExpression"]]],
   native=FileNameJoin[{run,"checkpoints","native",p,"block_"<>IntegerString[id,10,5]}];
   audit=Import[FileNameJoin[{native,"audit.json"}],"RawJSON"];
   require[AssociationQ[audit]&&TrueQ[audit["passed"]]&&audit["method"]==="DirectRationalSeries-1.0"&&audit["source_sha256"]===sourcehash,"Native rational audit missing or mismatched: "<>src];
   Do[key="jet"<>ToString[order];path=FileNameJoin[{native,audit["files"][key]["name"]}];
     require[hashFile[path]===audit["files"][key]["sha256"],"Changed native rational series: "<>path];
     rf=nativeParse[Import[path,"Text"]];
     require[AssociationQ[rf]&&rf["Order"]===order&&IntegerQ[rf["LowerBound"]]&&AssociationQ[rf["Jet"]],"Invalid rational series record."];
     require[And@@(IntegerQ[#]&&rf["LowerBound"]<=#<=order&/@Keys[rf["Jet"]]),"Invalid rational exponent range."];
     require[FreeQ[Values[rf["Jet"]],beta|x|z|m2|mu2|nc|nf|_Log|_PolyLog|_MPL22],"Nonrational term in a native coefficient."];
     If[order===4,j4=rf["Jet"],j6=rf["Jet"]],{order,{4,6}}];
   Do[nativeExactZero[Lookup[j4,kk,0]-Lookup[j6,kk,0],p<>" block "<>ToString[id]<>" power "<>ToString[kk]],
      {kk,Union[Keys[j4],Select[Keys[j6],#<=4&]]}];
   refs=<||>;
   Do[path=FileNameJoin[{native,ref["File"]}];require[hashFile[path]===ref["SHA256"],"Changed exact reference function."];
      bb=ref["Numerator"]/ref["Denominator"];rf=nativeParse[Import[path,"Text"]];
      require[FreeQ[rf,beta|x|z|m2|mu2|nc|nf|_Log|_PolyLog|_MPL22],"Nonrational exact reference function."];
      AssociateTo[refs,bb->rf],{ref,audit["references"]}];
   require[Sort[Keys[refs]]===Sort[{1/1000,1/500,1/200,1/100,1/50,1/20,1/10,1/5,3/10}],"Exact references have the wrong beta grid."];
   rv=If[Length[j6]==0,7,Min[Keys[j6]]];
   rec=<|"Version"->backendVersion,"SourceBlockHash"->sourcehash,"Block"->id,"Projector"->p,"NfPower"->k,"Monomial"->pairs,
     "RationalMethod"->"DirectRationalSeries-1.0","Jet4"->j4,"Jet6"->j6,"Valuation"->rv,"ExactReferenceFunctions"->refs,
     "RationalBeta"->Missing["NotConstructed"],
     "Status"->"Empty series: zero through the stated order."|>;
   AppendTo[rows,rec]
 ];
 If[Mod[j,50]==0||j==Length[index],say[p<>" rational inputs "<>ToString[j]<>"/"<>ToString[Length[index]]<>"; retained "<>ToString[Length[rows]]]],{line,index}];
 AssociateTo[allBlocks,p->rows];
 save[FileNameJoin[{run,"results",p<>"_rational_inputs.wl"}],<|"Functions"->specs,"Blocks"->rows,"Version"->backendVersion,
 "Scope"->"Cached full rational expressions or new exact truncated rational series and separately computed unexpanded rational functions at the validation beta values. Not sampled fits."|>];rows
];

rationalSeriesFor[r_Association,order_Integer] := If[
 Lookup[r,"RationalMethod",""]==="DirectRationalSeries-1.0",r["Jet"<>ToString[order]],rationalJet[r["RationalBeta"],order]];

rationalReferenceFor[r_Association,bb_,cc_] := Module[{q,v},
 q=If[Lookup[r,"RationalMethod",""]==="DirectRationalSeries-1.0",
    require[KeyExistsQ[r["ExactReferenceFunctions"],bb],"No exact unexpanded reference at beta="<>ToString[bb,InputForm]];
    r["ExactReferenceFunctions"][bb],r["RationalBeta"]/.beta->bb];
 v=Quiet[Check[q/.costh->cc,$Failed],{Power::infy,Infinity::indet}];
 If[v===$Failed||!finiteQ[v],
   v=TimeConstrained[Limit[q,costh->cc],180,$Aborted];
   require[finiteQ[v]&&FreeQ[v,costh],"Unexpanded rational coefficient has no established finite angular limit, block "<>ToString[r["Block"]]]
 ];
 require[FreeQ[v,beta|costh|_Real]&&finiteQ[v],"An exact rational reference is not a finite exact number."];v
];

componentInfo[label_String]:=Module[{t=StringSplit[label,"_nf"]},{t[[1]],ToExpression[t[[2]]]}];
componentBlocks[label_String]:=With[{q=componentInfo[label]},Select[allBlocks[q[[1]]],#["NfPower"]===q[[2]]&]];

(* Primitive-library check retained from the construction run. *)
checkPrimitiveLibrary[]:=Module[{specs,entries={},sp,ij,aa,ww,h=6,jexpr,rv,sv,err,cp,key,cc,bb=1/1000,count=0},
 specs=DeleteDuplicates[Flatten[Values[functions]]];specs=Select[specs,#["Kind"]=!="Flavour"&&#["Kind"]=!="Constant"&];
 Do[count++;key=IntegerString[Hash[{backendVersion,sp,h},"SHA256"],16,64];cp=FileNameJoin[{run,"checkpoints","primitive_tests",key<>".wl"}];
 If[FileExistsQ[cp],entries=Join[entries,load[cp]],
 ij=primitiveJet[sp,h];jexpr=jetExpression[ij]/.Lb->Log[beta];rv={};
 Do[aa=primitiveNumber[sp,bb,cc,15625,15625,3,60];
 sv=finiteExpressionNumber[jexpr/.{beta->bb,costh->cc,m2->15625,mu2->15625,nc->3},60];
 require[finiteQ[sv],"Primitive jet has unevaluated numerical constants."];err=Abs[sv-aa]/Max[1,Abs[aa]];
 AppendTo[rv,<|"Specification"->sp,"Beta"->bb,"CosTheta"->cc,"Order"->h,"ScaledError"->err,"Pass"->TrueQ[err<10^-12]|>],{cc,{-1/3,1/3}}];
 save[cp,rv];entries=Join[entries,rv];require[And@@Lookup[rv,"Pass"],"A primitive expansion disagrees with its independent numerical evaluation."]];
 If[Mod[count,10]==0||count==Length[specs],say["Primitive tests "<>ToString[count]<>"/"<>ToString[Length[specs]]]],{sp,specs}];
 primitiveChecks=entries;save[FileNameJoin[{run,"results","primitive_numerical_checks.wl"}],entries];
 require[And@@Lookup[entries,"Pass"],"A cached primitive test is not successful."];True
];

(* Assemble one colour/flavour component from the saved rational and MPL pieces. *)
computeComponent[label_String,order_Integer]:=Module[{parts=componentBlocks[label],p=componentInfo[label][[1]],sum=<||>,cp,rec,j=0,v,h,product,rat,term,sp,k,e,full,coeff,versionKey},
 require[Length[parts]>0,"No nonzero blocks for "<>label];
 versionKey={backendVersion,"DirectRationalSeries-1.0",label,order};
 Do[j++;cp=FileNameJoin[{run,"checkpoints","assembly",label<>"_p"<>ToString[order],"term_"<>IntegerString[r["Block"],10,5]<>".wl"}];
 If[FileExistsQ[cp],rec=load[cp];require[rec["Key"]===versionKey,"Stale expansion checkpoint."];term=rec["Jet"],
 v=r["Valuation"];h=order-v;
 If[v>order,term=<||>,
 product=<|0->1|>;
 Do[sp=functions[p][[kv[[1]]+1]];k=kv[[2]];
 If[sp["Kind"]==="Constant",product=jetScale[product,sp["Value"]^k],require[k>=0,"Negative special-function power."];product=jetMul[product,jetPower[primitiveJet[sp,h],k,h],h]],{kv,r["Monomial"]}];
 require[And@@(IntegerQ[#]&&#>=0&/@Keys[product]),"Negative primitive power would require deeper rational coefficients."];
  term=jetScale[jetMul[rationalSeriesFor[r,order],product,order],m2];
 term=Map[coefficientCanonical,term]];
 require[FreeQ[Values[term],_Series|_Derivative|_MPL22],"Unevaluated primitive in coefficient block."];
 save[cp,<|"Key"->versionKey,"Block"->r["Block"],"Jet"->term|>]];
 sum=jetAdd[sum,term];
 If[Mod[j,50]==0||j==Length[parts],sum=Map[coefficientCanonical,sum];say[label<>" beta^"<>ToString[order]<>" blocks "<>ToString[j]<>"/"<>ToString[Length[parts]]]],{r,parts}];
 sum=Map[coefficientCanonical,sum];sum=KeySelect[sum,#<=order&];
 full=jetExpression[sum]/.Lb->Log[beta];
 require[FreeQ[full,Derivative|Series|MPL22],"Unresolved special-function derivative in assembled result."];
 rec=<|"Label"->label,"Order"->order,"Jet"->sum,"Expression"->full,"BoundaryConstants"->DeleteDuplicates[Cases[full,_MB,Infinity]],"Version"->backendVersion,
 "Status"->"Threshold series. MB[p,q,a,w] are the boundary integrals."|>;
 save[FileNameJoin[{run,"results",label<>"_series"<>ToString[order]<>".wl"}],rec];rec
];

checkComponent[label_String,r4_Association,r6_Association]:=Module[{diff,sym,rows,lo,logs,v,symaudit},
 diff=coefficientCanonical[r4["Expression"]-jetExpression[KeySelect[r6["Jet"],#<=4&]]/.Lb->Log[beta]];
 save[FileNameJoin[{run,"results",label<>"_common_order_difference.wl"}],diff];require[zeroQ[diff],"Common-order coefficient mismatch for "<>label];

 sym=TimeConstrained[FullSimplify[coefficientCanonical[r6["Expression"]-(r6["Expression"]/.costh->-costh)],Assumptions->assumptions],120,$Aborted];
 symaudit=If[zeroQ[sym],"PROVED",If[sym===$Aborted,"NOT_PROVED_TIMEOUT","NOT_REDUCED_TO_ZERO"]];
 If[symaudit=!="PROVED",warn[label<>" exact angular symmetry not shown symbolically; residual saved and checked numerically."]];
 lo=Select[Keys[r6["Jet"]],#<0&&!zeroQ[r6["Jet"][#]]&];logs=Max[0,Sequence@@(Exponent[#,Lb]&/@Values[r6["Jet"]])];
 rows=<|"Label"->label,"CommonOrderPass"->True,"AngularSymmetry"->symaudit,"AngularDifference"->sym,
 "PreservedNegativePowers"->lo,"MaxLogPower"->logs,"BoundaryConstants"->r6["BoundaryConstants"]|>;
 AppendTo[algebraRows,rows];save[FileNameJoin[{run,"results","component_symbolic_checks.wl"}],algebraRows];rows
];

unexpandedComponent[label_String,bb_,cc_,par_List,wp_Integer]:=Module[{p=componentInfo[label][[1]],bl=componentBlocks[label],sp,values,terms,q,mm=par[[1]],uu=par[[2]],nn=par[[3]],ids},
 ids=Union[Flatten[(#["Monomial"][[All,1]]&/@bl)]];
 values=Association[Table[id->primitiveNumber[functions[p][[id+1]],bb,cc,mm,uu,nn,wp],{id,ids}]];
 terms=Table[q=rationalReferenceFor[r,bb,cc];
 N[mm q,wp] Times@@(values[#[[1]]]^#[[2]]&/@r["Monomial"]),{r,bl}];
 q=Block[{$MaxExtraPrecision=5wp},Total[terms]];require[finiteQ[q],"Unexpanded coefficient is not finite."];q
];

stableReference[label_String,bb_,cc_,par_List]:=Module[{low,high,lowwp=80,err,wp},
 low=unexpandedComponent[label,bb,cc,par,lowwp];
 Do[high=unexpandedComponent[label,bb,cc,par,wp];err=scaled[low,high];
 If[TrueQ[err<10^-25],Return[<|"Value"->high,"LowerValue"->low,"Change"->err,"WorkingPrecision"->wp,"LowerPrecision"->lowwp|>]];
 say["Increasing reference precision for "<>label<>" at "<>ToString[{bb,cc},InputForm]];low=high;lowwp=wp,{wp,{120,160,220,300}}];
 require[False,"Reference did not stabilize through 300 working digits: "<>label]
];

seriesNumber[e_,bb_,cc_,par_List,wp_Integer]:=Module[{v},v=finiteExpressionNumber[e/.{beta->bb,costh->cc,m2->par[[1]],mu2->par[[2]],nc->par[[3]]},wp];require[finiteQ[v],"Series evaluation left unresolved constants."];v];
relativeError[a_,d_]:=If[Abs[d]>10^-40,Abs[a-d]/Abs[d],Missing["NearZeroReference"]];

(* Original validation grid, retained so the saved calculation can be reproduced. *)
validateAll[]:=Module[{bs={1/1000,1/500,1/200,1/100,1/50,1/20,1/10,1/5,3/10},cs={-2/3,-1/3,0,1/3,2/3},pars={{15625,15625,3},{15625,62500,3},{22500,5625,4}},par,bb,cc,key,cp,pointRows,r80,r120,stable,refs=<||>,refs80=<||>,refPrecisions=<||>,refInfo,full80,fullStability,fullWP,highFull,wp,s4,s6,labels=Join[componentLabels,knownLabels],expr4,expr6,fullref,full4,full6,rel,pass,guardRows={}},
 Do[key={backendVersion,par,bb,cc};cp=FileNameJoin[{run,"checkpoints","validation",IntegerString[Hash[key,"SHA256"],16,64]<>".wl"}];
 If[FileExistsQ[cp],pointRows=load[cp];validationRows=Join[validationRows,pointRows];Continue[]];
 say["Full validation "<>ToString[{par,bb,cc},InputForm]];refs=<||>;refs80=<||>;refPrecisions=<||>;pointRows={};
 Do[refInfo=stableReference[label,bb,cc,par];r80=refInfo["LowerValue"];r120=refInfo["Value"];stable=refInfo["Change"];
 AssociateTo[refs,label->r120];AssociateTo[refs80,label->r80];AssociateTo[refPrecisions,label->refInfo["WorkingPrecision"]];
 If[MemberQ[componentLabels,label],expr4=computed[label]["p4"]["Expression"];expr6=computed[label]["p6"]["Expression"],expr4=Normal[known[label]["p4"]];expr6=Normal[known[label]["p6"]]];
 s4=seriesNumber[expr4,bb,cc,par,Max[110,refInfo["WorkingPrecision"]]];s6=seriesNumber[expr6,bb,cc,par,Max[110,refInfo["WorkingPrecision"]]];rel=relativeError[s6,r120];
 AppendTo[pointRows,<|"Quantity"->label,"Parameters"->par,"Beta"->bb,"CosTheta"->cc,"ReUnexpanded"->Re[r120],"ImUnexpanded"->Im[r120],"ReSeries4"->Re[s4],"ImSeries4"->Im[s4],"ReSeries6"->Re[s6],"ImSeries6"->Im[s6],"AbsoluteError4"->Abs[s4-r120],"AbsoluteError6"->Abs[s6-r120],"RelativeError4"->relativeError[s4,r120],"RelativeError6"->rel,"ReferenceLastChange"->stable,"ReferenceWorkingPrecision"->refInfo["WorkingPrecision"]|>];

 If[bb===1/1000,require[TrueQ[Abs[s6-r120]/Max[1,Abs[r120]]<10^-8],"Small-beta sixth-order comparison fails for "<>label]],{label,labels}];
 Do[fullref=refs[p<>"_nf0"]+5refs[p<>"_nf1"]+25refs[p<>"_nf2"];full80=refs80[p<>"_nf0"]+5refs80[p<>"_nf1"]+25refs80[p<>"_nf2"];fullStability=scaled[full80,fullref];fullWP=Max[Lookup[refPrecisions,{p<>"_nf0",p<>"_nf1",p<>"_nf2"}]];
 If[!TrueQ[fullStability<10^-25],Do[If[wp<=fullWP,Continue[]];highFull=unexpandedComponent[p<>"_nf0",bb,cc,par,wp]+5unexpandedComponent[p<>"_nf1",bb,cc,par,wp]+25unexpandedComponent[p<>"_nf2",bb,cc,par,wp];fullStability=scaled[fullref,highFull];fullref=highFull;fullWP=wp;If[TrueQ[fullStability<10^-25],Break[]],{wp,{160,220,300}}]];
 require[TrueQ[fullStability<10^-25],"The assembled reference is not precision-stable through the allowed precision cap: "<>p];
 full4=seriesNumber[assembled[p]["p4"],bb,cc,par,Max[110,fullWP]];full6=seriesNumber[assembled[p]["p6"],bb,cc,par,Max[110,fullWP]];
 AppendTo[pointRows,<|"Quantity"->p<>"_FULL_nf5","Parameters"->par,"Beta"->bb,"CosTheta"->cc,"ReUnexpanded"->Re[fullref],"ImUnexpanded"->Im[fullref],"ReSeries4"->Re[full4],"ImSeries4"->Im[full4],"ReSeries6"->Re[full6],"ImSeries6"->Im[full6],"AbsoluteError4"->Abs[full4-fullref],"AbsoluteError6"->Abs[full6-fullref],"RelativeError4"->relativeError[full4,fullref],"RelativeError6"->relativeError[full6,fullref],"ReferenceLastChange"->fullStability,"ReferenceWorkingPrecision"->fullWP|>];
 If[bb===1/1000,require[TrueQ[Abs[full6-fullref]/Max[1,Abs[fullref]]<10^-8],"Small-beta assembled comparison fails for "<>p]],{p,{"P1","P2"}}];
 save[cp,pointRows];validationRows=Join[validationRows,pointRows];save[FileNameJoin[{run,"results","validation_all.wl"}],validationRows],{par,pars},{cc,cs},{bb,bs}];
 True
];

checkNumericalSymmetry[]:=Module[{groups,rows={},g,pos,neg,vp,vn,err},
 groups=GatherBy[validationRows,{#["Quantity"],#["Parameters"],#["Beta"]}&];
 Do[Do[neg=SelectFirst[g,#["CosTheta"]===-pos["CosTheta"]&,Missing["Absent"]];
 require[AssociationQ[neg],"An angular-symmetry partner is missing from the numerical grid."];
 vp=pos["ReUnexpanded"]+I pos["ImUnexpanded"];vn=neg["ReUnexpanded"]+I neg["ImUnexpanded"];err=scaled[vp,vn];
 AppendTo[rows,<|"Quantity"->pos["Quantity"],"Parameters"->pos["Parameters"],"Beta"->pos["Beta"],"PositiveCosTheta"->pos["CosTheta"],"ScaledDifference"->err,"Pass"->TrueQ[err<10^-25]|>],
 {pos,Select[g,#["CosTheta"]>0&]}],{g,groups}];
 save[FileNameJoin[{run,"results","numerical_angular_symmetry.wl"}],rows];require[And@@Lookup[rows,"Pass"],"Numerical angular symmetry failed. Inspect the saved differences."];rows
];
maximumRelative[values_List]:=With[{v=Select[values,NumberQ]},If[Length[v]==0,Missing["NoDefinedRelativeError"],Max[v]]];

plotAndExport[]:=Module[{qs,rs,groups,table,plot,g,cols},
 cols=Keys[First[validationRows]];
 Export[FileNameJoin[{run,"results","validation_all.csv"}],Prepend[(Lookup[#,cols]/._Missing->"NA")&/@validationRows,cols],"CSV"];
 qs=DeleteDuplicates[Lookup[validationRows,"Quantity"]];
 Do[rs=Select[validationRows,#["Quantity"]===q&];groups=GatherBy[rs,#["Beta"]&];
 table=SortBy[Table[{First[g]["Beta"],maximumRelative[Lookup[g,"RelativeError4"]],maximumRelative[Lookup[g,"RelativeError6"]]},{g,groups}],First];
 save[FileNameJoin[{run,"results",q<>"_error_table.wl"}],<|"Data"->table,"Scope"->"Maxima over the sampled angles and parameter sets."|>];
 If[Environment["GGAA_UNIFIED_NO_PLOTS"]=!="1",plot=ListLogLogPlot[{Select[table[[All,{1,2}]],NumberQ[Last[#]]&&Last[#]>0&],Select[table[[All,{1,3}]],NumberQ[Last[#]]&&Last[#]>0&]},Joined->True,Frame->True,FrameLabel->{"beta_A","Maximum sampled relative error"},PlotLabel->q,PlotLegends->{"through beta_A^4","through beta_A^6"},ImageSize->700];
 If[TimeConstrained[Check[Export[FileNameJoin[{run,"results",q<>"_errors.png"}],plot],$Failed],300,$Failed]===$Failed,warn["Plot export failed: "<>q<>". Numerical table preserved."]]],{q,qs}];True
];

End[];
