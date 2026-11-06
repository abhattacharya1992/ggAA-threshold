(* Numerical validation of the saved series against the unexpanded expressions. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "ArchivedDefinitions.wl"}]];
Begin["AAUnified`"];

validationVersion = "GGAA-Saved-Validation-1.0";
validationStart = DateString[{"ISODate", "T", "Time"}];
config = Import[FileNameJoin[{run, "validation_config.json"}], "RawJSON"];
inputFingerprint = config["input_fingerprint"];
sourceRun = config["source_run"];
referenceTolerance = 10^-25;
minimumScaledDigits = 30;
smallBetaTolerance = 10^-8;
validationBetas = {1/1000, 1/500, 1/200, 1/100, 1/50, 1/20, 1/10, 1/5, 3/10};
validationAngles = {-2/3, -1/3, 0, 1/3, 2/3};
validationParameters = {{15625,15625,3}, {15625,62500,3}, {22500,5625,4}};
validationLabels = Join[componentLabels, knownLabels];
validationRows = {}; referenceTraces = {}; symmetryRows = {};
inheritedChecks = {}; inheritedWarnings = {}; restoreChecks = {};
scope = "Numerical validation of the saved series through beta^4 and beta^6 against the unexpanded expressions; full sums at nf=5.";

Clear[save,load];
save[path_String,value_] := Module[{tmp=path<>".tmp."<>ToString[$ProcessID], ok},
 makeDir[DirectoryName[path]];
 ok=Check[Put[value,tmp]; True,False];
 require[TrueQ[ok] && FileExistsQ[tmp],"Cannot write "<>path];
 ok=Check[RenameFile[tmp,path,OverwriteTarget->True];True,False];
 require[TrueQ[ok] && FileExistsQ[path],"Cannot commit "<>path]; path
];
load[path_String] := Module[{v},
 require[FileExistsQ[path],"Missing saved file: "<>path];
 v=Block[{$Context="AAUnified`",$ContextPath={"AAUnified`","System`"}},Check[Get[path],$Failed]];
 require[v=!=$Failed,"Cannot read "<>path];v
];

(* Validation-only mode: construction routines are deliberately unavailable. *)
Clear[loadAndReduce,computeComponent,checkPrimitiveLibrary];
loadAndReduce[___] := require[False,"Disabled in validation-only mode: rational construction."];
computeComponent[___] := require[False,"Disabled in validation-only mode: symbolic expansion."];
checkPrimitiveLibrary[___] := require[False,"Primitive-series construction is disabled; archived tests are retained."];

hashKey[key_] := IntegerString[Hash[key,"SHA256"],16,64];
cachePath[group_String,key_] := FileNameJoin[{run,"checkpoints",group,hashKey[key]<>".wl"}];
checkpointKey[args___] := {validationVersion,inputFingerprint,args};

scaledDigits[v_] := If[!finiteQ[v],-Infinity,Min[
 Accuracy[Re[v]]+Log[10,Max[1,Abs[Re[v]]]],
 Accuracy[Im[v]]+Log[10,Max[1,Abs[Im[v]]]]]];

(* Accept a numerical reference only after two successive precision levels agree. *)
precisionLoop[evaluator_,precisions_List,tol_,observer_] := Module[
 {low,high,lowwp,wp,err,record,history={},answer=Missing["NotAccepted"],pass},
 require[Length[precisions]>=2 && And@@(IntegerQ[#]&&#>0&/@precisions),"Bad precision sequence."];
 lowwp=First[precisions];low=evaluator[lowwp];
 If[!finiteQ[low],
   answer=Failure["NonfiniteInitialValue",<|"Precision"->lowwp,"Value"->low|>],
   Do[
    high=evaluator[wp];
    If[!finiteQ[high],
     record=<|"LowerPrecision"->lowwp,"WorkingPrecision"->wp,"LowerValue"->low,
       "Value"->high,"Pass"->False,"Reason"->"NonfiniteValue"|>;
     AppendTo[history,record];observer[record];
     answer=Failure["NonfiniteValue",record];Break[]
    ];
    err=scaled[low,high];
    pass=TrueQ[err<tol] && TrueQ[scaledDigits[low]>=minimumScaledDigits] &&
         TrueQ[scaledDigits[high]>=minimumScaledDigits];
    record=<|"LowerPrecision"->lowwp,"WorkingPrecision"->wp,"LowerValue"->low,
      "Value"->high,"Change"->err,"LowerScaledDigits"->scaledDigits[low],
      "HigherScaledDigits"->scaledDigits[high],"Tolerance"->tol,"Pass"->pass|>;
    AppendTo[history,record];observer[record];
    If[pass,answer=KeyDrop[record,{"Pass"}];Break[]];
    low=high;lowwp=wp,
    {wp,Rest[precisions]}]
 ];
 If[MissingQ[answer],answer=Failure["PrecisionSequenceExhausted",<|"Precisions"->precisions,"History"->history|>]];
 <|"Accepted"->AssociationQ[answer],"Result"->answer,"History"->history|>
];

controlFlowTests[] := Module[{rows={},r,calls={},ev,obs=Function[v,Null]},
 ev=Function[w,AppendTo[calls,w]; 2+3 I];
 r=precisionLoop[ev,{80,120,160},referenceTolerance,obs];
 AppendTo[rows,<|"Test"->"Accept at first stable pair, no fall-through", "Pass"->TrueQ[r["Accepted"] && calls==={80,120} && r["Result"]["WorkingPrecision"]===120]|>];
 calls={};ev=Function[w,AppendTo[calls,w]; If[w===80,1+1/100,1]];
 r=precisionLoop[ev,{80,120,160,220},referenceTolerance,obs];
 AppendTo[rows,<|"Test"->"Escalate once then return accepted pair", "Pass"->TrueQ[r["Accepted"] && calls==={80,120,160} && r["Result"]["WorkingPrecision"]===160]|>];
 calls={};ev=Function[w,AppendTo[calls,w]; 1+1/w];
 r=precisionLoop[ev,{80,120,160},referenceTolerance,obs];
 AppendTo[rows,<|"Test"->"Exhausted sequence fails", "Pass"->TrueQ[!r["Accepted"] && FailureQ[r["Result"]] && calls==={80,120,160}]|>];
 ev=Function[w,If[w===80,0,referenceTolerance]];
 r=precisionLoop[ev,{80,120},referenceTolerance,obs];
 AppendTo[rows,<|"Test"->"Tolerance equality does not pass strict comparison", "Pass"->TrueQ[!r["Accepted"]]|>];
 ev=Function[w,If[w===80,1,Indeterminate]];
 r=precisionLoop[ev,{80,120,160},referenceTolerance,obs];
 AppendTo[rows,<|"Test"->"Nonfinite evaluation is rejected", "Pass"->TrueQ[!r["Accepted"] && FailureQ[r["Result"]]]|>];
 ev=Function[w,N[1,10]];
 r=precisionLoop[ev,{80,120},referenceTolerance,obs];
 AppendTo[rows,<|"Test"->"Identical low-precision values are not certification", "Pass"->TrueQ[!r["Accepted"]]|>];
 save[FileNameJoin[{run,"results","reference_control_flow_tests.wl"}],rows];
 require[And@@Lookup[rows,"Pass"],"A reference-control-flow regression test failed."];
 say["All six precision-control regression tests passed."];rows
];

(* Restore the symbolic series and the checks already completed in the construction run. *)
restoreSaved[] := Module[{r,full,expected,checks={},diff,labels=validationLabels,sp},
 require[AssociationQ[config]&&config["version"]===validationVersion,"Wrong continuation configuration."];
 require[backendVersion==="GGAA-Unified-RealMPL22-0.1","Unsupported mathematical engine version."];
 Do[
  r=load[FileNameJoin[{run,"input","saved",p<>"_rational_inputs.wl"}]];
  require[AssociationQ[r]&&r["Version"]===backendVersion&&ListQ[r["Blocks"]]&&ListQ[r["Functions"]],"Bad saved rational input for "<>p];
  AssociateTo[allBlocks,p->r["Blocks"]];AssociateTo[functions,p->r["Functions"]];
  require[And@@(AssociationQ[#]&&MemberQ[{0,1,2},#["NfPower"]]&&#["Projector"]===p&/@allBlocks[p]),"Invalid rational block metadata for "<>p];
  say["Restored "<>p<>": "<>ToString[Length[allBlocks[p]]]<>" rational blocks; no reductions run."],{p,{"P1","P2"}}];
 Do[AssociateTo[known,label-><|"p4"->readKnown[label<>"_p4.wl"],"p6"->readKnown[label<>"_p6.wl"]|>],{label,knownLabels}];
 Do[
  AssociateTo[computed,label->Association[Table[
    r=load[FileNameJoin[{run,"input","saved",label<>"_series"<>ToString[ord]<>".wl"}]];
    require[AssociationQ[r]&&r["Label"]===label&&r["Order"]===ord&&r["Version"]===backendVersion,"Saved component record mismatch: "<>label];
    "p"<>ToString[ord]->r,{ord,{4,6}}]]],{label,componentLabels}];
 assembledSymbolic=<||>;
 Do[
  full=load[FileNameJoin[{run,"input","saved",p<>"_full_series.wl"}]];
  require[AssociationQ[full]&&KeyExistsQ[full,"ThroughBeta4"]&&KeyExistsQ[full,"ThroughBeta6"],"Missing saved full form factor "<>p];
  AssociateTo[assembledSymbolic,p-><|"p4"->full["ThroughBeta4"],"p6"->full["ThroughBeta6"]|>];
  Do[
   expected=computed[p<>"_nf0"]["p"<>ToString[ord]]["Expression"]+
    nf If[p==="P1",Normal[known["P1_nf1"]["p"<>ToString[ord]]],computed["P2_nf1"]["p"<>ToString[ord]]["Expression"]]+
    nf^2 Normal[known[p<>"_nf2"]["p"<>ToString[ord]]];
   diff=assembledSymbolic[p]["p"<>ToString[ord]]-expected;
   If[diff=!=0,diff=TimeConstrained[coefficientCanonical[diff],180,$Aborted]];
   AppendTo[checks,<|"Projector"->p,"Order"->ord,"Test"->"Saved full expression equals saved component sum","Difference"->diff,"Pass"->TrueQ[diff===0]|>],
   {ord,{4,6}}],{p,{"P1","P2"}}];
 save[FileNameJoin[{run,"results","restored_assembly_checks.wl"}],checks];
 require[And@@Lookup[checks,"Pass"],"Saved assembly consistency failed; no coefficients were changed."];
 assembled=assembledSymbolic/.nf->5;
 inheritedChecks=load[FileNameJoin[{run,"input","saved","component_symbolic_checks.wl"}]];
 primitiveChecks=load[FileNameJoin[{run,"input","saved","primitive_numerical_checks.wl"}]];
 require[ListQ[primitiveChecks]&&Length[primitiveChecks]>0&&And@@Lookup[primitiveChecks,"Pass"],"Archived primitive checks are not all successful."];
 Do[If[Lookup[r,"AngularSymmetry",""]=!="PROVED",AppendTo[inheritedWarnings,r["Label"]<>": exact angular identity remains "<>ToString[r["AngularSymmetry"]]]],{r,inheritedChecks}];
 save[FileNameJoin[{run,"results","inherited_analytical_status.wl"}],<|"Checks"->inheritedChecks,"Warnings"->inheritedWarnings,
   "Meaning"->"Status carried over from the construction run."|>];
 require[And@@(Length[componentBlocks[#]]>0&/@labels),"A flavour coefficient lacks reference blocks."];
 restoreChecks=checks;
 say["Restored all component and full series. Four saved-assembly checks passed. Numerical validation only."];
 True
];

cachedBoundary[p_Integer,q_Integer,a_,w_,wp_Integer] := Module[{key,path,r,v},
 key={backendVersion,"MBNumber",p,q,a,w,wp};path=cachePath["boundary_numbers",key];
 If[FileExistsQ[path],r=load[path];require[r["Key"]===key&&finiteQ[r["Value"]],"Invalid boundary cache."];v=r["Value"],
  v=mpNumber[p,q,a,w,wp];require[finiteQ[v],"Nonfinite boundary constant."];
  save[path,<|"Key"->key,"Value"->v|>]];v
];
Clear[finiteExpressionNumber];
finiteExpressionNumber[e_,wp_Integer] := Module[{atoms,rules,v},
 atoms=DeleteDuplicates[Cases[e,_MB,Infinity]];
 rules=(#->cachedBoundary[#[[1]],#[[2]],#[[3]],#[[4]],wp+10])&/@atoms;
 v=Block[{$MaxExtraPrecision=5wp},N[e/.rules,wp]];
 require[finiteQ[v],"Expression retains unresolved numerical functions."];v
];

componentValue[label_String,bb_,cc_,par_List,wp_Integer] := Module[{key,path,r,v,t},
 key=checkpointKey["Unexpanded",label,bb,cc,par,wp];path=cachePath["component_values",key];
 If[FileExistsQ[path],r=load[path];require[r["Key"]===key&&finiteQ[r["Value"]],"Invalid unexpanded-value cache."];v=r["Value"],
  t=AbsoluteTime[];say["Unexpanded "<>label<>" at "<>ToString[{par,bb,cc},InputForm]<>"; digits="<>ToString[wp]];
  v=unexpandedComponent[label,bb,cc,par,wp];require[finiteQ[v],"Nonfinite unexpanded coefficient."];
  save[path,<|"Key"->key,"Value"->v,"Seconds"->AbsoluteTime[]-t|>]];v
];

(* Common precision-escalation wrapper for component and full-form-factor comparisons. *)
acceptedComparison[kind_String,id_,bb_,cc_,par_List,precisions_List,evaluator_] := Module[
 {key,path,tracePath,r,report,history={},observer,result},
 key=checkpointKey[kind,id,bb,cc,par,precisions,referenceTolerance,minimumScaledDigits];
 path=cachePath["accepted_precision_checks",key];tracePath=cachePath["precision_traces",key];
 If[FileExistsQ[path],
  r=load[path];require[r["Key"]===key&&AssociationQ[r["Result"]],"Mismatched accepted comparison."];
  result=r["Result"];
  require[finiteQ[result["Value"]]&&finiteQ[result["LowerValue"]]&&
    TrueQ[scaled[result["LowerValue"],result["Value"]]<referenceTolerance]&&
    TrueQ[Min[scaledDigits[result["LowerValue"]],scaledDigits[result["Value"]]]>=minimumScaledDigits],
    "Cached precision acceptance no longer meets its criterion."],
  observer=Function[rec,
    AppendTo[history,rec];save[tracePath,<|"Key"->key,"History"->history|>];
    say["Precision comparison "<>kind<>" "<>ToString[id,InputForm]<>" "<>
      ToString[{bb,cc},InputForm]<>" at "<>ToString[Lookup[rec,"LowerPrecision"]]<>"/"<>
      ToString[Lookup[rec,"WorkingPrecision"]]<>" digits; accepted="<>ToString[Lookup[rec,"Pass",False]]]
  ];
  report=precisionLoop[evaluator,precisions,referenceTolerance,observer];
  save[tracePath,Append[report,"Key"->key]];
  require[TrueQ[report["Accepted"]],"Precision acceptance failed for "<>kind<>" "<>ToString[id,InputForm]<>". Inspect "<>tracePath];
  result=report["Result"];save[path,<|"Key"->key,"Result"->result,"History"->report["History"]|>]
 ];result
];

Clear[stableReference];
stableReference[label_String,bb_,cc_,par_List] := acceptedComparison[
 "UnexpandedComponent",label,bb,cc,par,{80,120,160,220,300},
 Function[w,componentValue[label,bb,cc,par,w]]
];
fullValue[p_String,bb_,cc_,par_List,wp_Integer] :=
 Block[{$MaxExtraPrecision=5wp},componentValue[p<>"_nf0",bb,cc,par,wp]+
   5 componentValue[p<>"_nf1",bb,cc,par,wp]+25 componentValue[p<>"_nf2",bb,cc,par,wp]];

seriesExpression[label_String,ord_Integer] := Which[
 MemberQ[componentLabels,label],computed[label]["p"<>ToString[ord]]["Expression"],
 MemberQ[knownLabels,label],Normal[known[label]["p"<>ToString[ord]]],
 label==="P1_FULL_nf5",assembled["P1"]["p"<>ToString[ord]],
 label==="P2_FULL_nf5",assembled["P2"]["p"<>ToString[ord]],
 True,require[False,"Unknown series label."]
];
seriesValue[label_String,ord_Integer,bb_,cc_,par_List,wp_Integer] := Module[{key,path,r,v,t},
 key=checkpointKey["SavedSeries",label,ord,bb,cc,par,wp];path=cachePath["series_values",key];
 If[FileExistsQ[path],r=load[path];require[r["Key"]===key&&finiteQ[r["Value"]],"Invalid saved-series value cache."];v=r["Value"],
  t=AbsoluteTime[];
  v=seriesNumber[seriesExpression[label,ord],bb,cc,par,wp];
  save[path,<|"Key"->key,"Value"->v,"Seconds"->AbsoluteTime[]-t|>]];v
];
seriesChecked[label_String,ord_Integer,bb_,cc_,par_List,refwp_Integer] := Module[{precisions,base},
 base=Max[110,refwp];precisions=DeleteDuplicates[Join[{base,base+40},Select[{160,220,300,340},#>base+40&]]];
 acceptedComparison["SavedSeries",{label,ord},bb,cc,par,precisions,
   Function[w,seriesValue[label,ord,bb,cc,par,w]]]
];

saveRows[] := (
 save[FileNameJoin[{run,"results","validation_all.wl"}],validationRows];
 Export[FileNameJoin[{run,"validation_progress.txt"}],"Numerical rows saved: "<>ToString[Length[validationRows]]<>" / 1080\n","Text"];
);

rowFor[label_String,bb_,cc_,par_List] := Module[{key,path,r,ref,s4,s6,v,a4,a6,row,isFull,p,smallpass},
 key=checkpointKey["ValidationRow",label,bb,cc,par];path=cachePath["validation_rows",key];
 If[FileExistsQ[path],r=load[path];require[r["Key"]===key&&AssociationQ[r["Row"]],"Mismatched numerical row checkpoint."];row=r["Row"],
  isFull=MemberQ[{"P1_FULL_nf5","P2_FULL_nf5"},label];
  ref=If[isFull,p=StringTake[label,2];
    acceptedComparison["UnexpandedFull",p,bb,cc,par,{80,120,160,220,300},Function[w,fullValue[p,bb,cc,par,w]]],
    stableReference[label,bb,cc,par]];
  s4=seriesChecked[label,4,bb,cc,par,ref["WorkingPrecision"]];
  s6=seriesChecked[label,6,bb,cc,par,ref["WorkingPrecision"]];
  v=ref["Value"];a4=s4["Value"];a6=s6["Value"];
  smallpass=If[bb===1/1000,TrueQ[Abs[a6-v]/Max[1,Abs[v]]<smallBetaTolerance],True];
  row=<|"Quantity"->label,"Parameters"->par,"Beta"->bb,"CosTheta"->cc,
   "ReUnexpanded"->Re[v],"ImUnexpanded"->Im[v],"ReSeries4"->Re[a4],"ImSeries4"->Im[a4],
   "ReSeries6"->Re[a6],"ImSeries6"->Im[a6],"AbsoluteError4"->Abs[a4-v],"AbsoluteError6"->Abs[a6-v],
   "AbsoluteRealError4"->Abs[Re[a4-v]],"AbsoluteImagError4"->Abs[Im[a4-v]],
   "AbsoluteRealError6"->Abs[Re[a6-v]],"AbsoluteImagError6"->Abs[Im[a6-v]],
   "RelativeError4"->relativeError[a4,v],"RelativeError6"->relativeError[a6,v],
   "ReferenceLastChange"->ref["Change"],"ReferenceWorkingPrecision"->ref["WorkingPrecision"],
   "ReferenceLowerPrecision"->ref["LowerPrecision"],"Series4LastChange"->s4["Change"],
   "Series6LastChange"->s6["Change"],"Series4WorkingPrecision"->s4["WorkingPrecision"],
   "Series6WorkingPrecision"->s6["WorkingPrecision"],"SmallBetaGateApplied"->(bb===1/1000),"SmallBetaGatePass"->smallpass|>;
  save[path,<|"Key"->key,"Row"->row|>]
 ];
 AppendTo[validationRows,row];saveRows[];
 say["Saved row "<>ToString[Length[validationRows]]<>"/1080: "<>label<>" "<>ToString[{par,bb,cc},InputForm]<>"; p6 relative error="<>ToString[N[row["RelativeError6"],6],InputForm]];
 require[TrueQ[row["SmallBetaGatePass"]],"Small-beta comparison fails for "<>label<>". The failed row has been saved; do not bypass this gate."];
 row
];

runNumericalGrid[] := Module[{labels=Join[validationLabels,{"P1_FULL_nf5","P2_FULL_nf5"}]},
 validationRows={};
 Do[
  say["POINT "<>ToString[{par,bb,cc},InputForm]];
  Do[rowFor[label,bb,cc,par],{label,labels}],
  {par,validationParameters},{cc,validationAngles},{bb,validationBetas}];
 require[Length[validationRows]===1080,"The 1080 numerical records were not completed."];
 True
];

runAngularChecks[] := Module[{groups,g,pos,neg,vr1,vr2,vs1,vs2,err,rows={},seriesrows={},e},
 groups=GatherBy[validationRows,{#["Quantity"],#["Parameters"],#["Beta"]}&];
 Do[
  Do[
   neg=SelectFirst[g,#["CosTheta"]===-pos["CosTheta"]&,Missing["Absent"]];
   require[AssociationQ[neg],"Missing angular partner in numerical data."];
   vr1=pos["ReUnexpanded"]+I pos["ImUnexpanded"];vr2=neg["ReUnexpanded"]+I neg["ImUnexpanded"];err=scaled[vr1,vr2];
   AppendTo[rows,<|"Quantity"->pos["Quantity"],"Parameters"->pos["Parameters"],"Beta"->pos["Beta"],"PositiveCosTheta"->pos["CosTheta"],"ScaledDifference"->err,"Pass"->TrueQ[err<10^-25]|>];
   Do[
    vs1=pos["ReSeries"<>ToString[ord]]+I pos["ImSeries"<>ToString[ord]];
    vs2=neg["ReSeries"<>ToString[ord]]+I neg["ImSeries"<>ToString[ord]];e=scaled[vs1,vs2];
    AppendTo[seriesrows,<|"Quantity"->pos["Quantity"],"Order"->ord,"Parameters"->pos["Parameters"],"Beta"->pos["Beta"],"PositiveCosTheta"->pos["CosTheta"],"ScaledDifference"->e,"Pass"->TrueQ[e<10^-25]|>],
    {ord,{4,6}}],{pos,Select[g,#["CosTheta"]>0&]}],{g,groups}];
 save[FileNameJoin[{run,"results","numerical_angular_symmetry.wl"}],rows];
 save[FileNameJoin[{run,"results","numerical_series_angular_symmetry.wl"}],seriesrows];
 require[And@@Lookup[rows,"Pass"],"Unexpanded numerical angular symmetry failed. Inspect saved results."];
 require[And@@Lookup[seriesrows,"Pass"],"Saved-series numerical angular symmetry failed. Inspect saved results."];
 symmetryRows=rows;True
];

exportTables[] := Module[{fields,qs,rs,groups,table,plot,png,g,nmissing},
 require[Length[validationRows]>0,"No validation records to export."];
 fields=Keys[First[validationRows]];
 Export[FileNameJoin[{run,"results","validation_all.csv"}],Prepend[(Lookup[#,fields]/._Missing->"NA")&/@validationRows,fields],"CSV"];
 qs=DeleteDuplicates[Lookup[validationRows,"Quantity"]];
 Do[
  rs=Select[validationRows,#["Quantity"]===q&];groups=GatherBy[rs,#["Beta"]&];
  table=SortBy[Table[{First[g]["Beta"],maximumRelative[Lookup[g,"RelativeError4"]],maximumRelative[Lookup[g,"RelativeError6"]],Length[g],Count[Lookup[g,"RelativeError6"],_Missing]},{g,groups}],First];
  save[FileNameJoin[{run,"results",q<>"_error_table.wl"}],<|"Columns"->{"Beta","MaxRelativeError4","MaxRelativeError6","SampleCount","UndefinedRelativeErrors"},"Data"->table,
    "Scope"->"Maxima over sampled angles and parameter choices only, not uniform bounds. Absolute errors and undefined ratios retained in validation_all.wl."|>];
  If[Environment["GGAA_UNIFIED_NO_PLOTS"]=!="1",
   plot=ListLogLogPlot[{Select[table[[All,{1,2}]],NumberQ[Last[#]]&&Last[#]>0&],Select[table[[All,{1,3}]],NumberQ[Last[#]]&&Last[#]>0&]},Joined->True,Frame->True,
     FrameLabel->{"beta_A","Maximum sampled relative error"},PlotLabel->q,PlotLegends->{"through beta_A^4","through beta_A^6"},ImageSize->700];
   png=TimeConstrained[Check[Export[FileNameJoin[{run,"results",q<>"_errors.png"}],plot],$Failed],300,$Failed];
   If[png===$Failed,warn["Plot export failed for "<>q<>"; numerical data retained."]]
  ],{q,qs}];True
];

validationMain[] := (
 stage["01_reference_control_flow_regression",controlFlowTests[]];
 stage["02_restore_saved_amplitudes_and_reference_blocks",restoreSaved[]];
 stage["03_component_and_full_form_factor_validation",runNumericalGrid[]];
 stage["04_numerical_angular_symmetry",runAngularChecks[]];
 stage["05_error_tables_and_plots",exportTables[]];
 stage["06_combined_validation_record",
  save[FileNameJoin[{run,"results","complete_unified_validation.wl"}],<|
    "SavedComponents"->computed,"PreviouslyCompleted"->known,"FullFormFactors"->assembledSymbolic,
    "RestoredAssemblyChecks"->restoreChecks,"InheritedAnalyticalChecks"->inheritedChecks,
    "ArchivedPrimitiveChecks"->primitiveChecks,"NumericalValidation"->validationRows,"NumericalAngularChecks"->symmetryRows,
    "Warnings"->warnings,"InheritedWarnings"->inheritedWarnings,"Scope"->scope,"ValidationVersion"->validationVersion,
    "EngineVersion"->backendVersion,"InputFingerprint"->inputFingerprint|>];True];
 True
);

outcome=CheckAbort[Catch[validationMain[],"UnifiedStop"],Failure["Interrupted",<|"Stage"->currentStage|>]];
If[outcome=!=True,say["FAILED: "<>ToString[outcome,InputForm]]];

finalization=CheckAbort[Catch[
 save[FileNameJoin[{run,"results","run_outcome.wl"}],<|"Outcome"->outcome,"Stages"->stages,"RowsSaved"->Length[validationRows],"Warnings"->warnings,"InheritedWarnings"->inheritedWarnings,"Scope"->scope|>];
 If[Length[validationRows]>0,saveRows[]];
 Export[FileNameJoin[{run,"summary.txt"}],StringRiffle[Join[{
  "gg -> AA validation-only continuation",
  "Execution: "<>If[outcome===True,"COMPLETED_SAVED_FORM_FACTOR_VALIDATION","FAILED_OR_INTERRUPTED"],
  "Kernel: "<>$Version,"Source recovery run: "<>sourceRun,"Numerical rows saved: "<>ToString[Length[validationRows]],
  "No rational reduction, symbolic expansion or native C++ compilation was rerun.",
  "Each reference value is accepted after its precision loop.",
  "Reference and series values must agree to 30 digits between working precisions.",
  "Results are written to this run directory.",
  "Saved series compared with the unexpanded expressions.",
  "Relative errors are measured on a finite grid; execution success does not certify a global accuracy target.",
  "Exact angular identities and boundary-constant relations are checked separately.",
  "STAGES:"},(#["Stage"]<>" : "<>#["Status"]&/@stages),{"NEW WARNINGS:"},If[warnings==={},{"None recorded."},warnings],
  {"INHERITED ANALYTICAL WARNINGS:"},If[inheritedWarnings==={},{"See the copied original report."},inheritedWarnings],
  {"Scope: "<>scope,"Started: "<>validationStart,"Finished: "<>DateString[{"ISODate","T","Time"}]}],"\n"]<>"\n","Text"];
 If[outcome===True,require[Length[validationRows]===1080,"Completion marker requires 1080 rows."];
  Export[FileNameJoin[{run,"COMPLETED_SAVED_VALIDATION.txt"}],"COMPLETED_SAVED_FORM_FACTOR_VALIDATION\nNumerical rows: 1080\n","Text"]];
 True,"UnifiedStop"],$Failed];
If[finalization=!=True,Print["FATAL: report finalization failed; existing checkpoint files remain."];Exit[3]];
say["Final outcome: "<>ToString[outcome,InputForm]];
Exit[If[outcome===True,0,1]];
End[];
