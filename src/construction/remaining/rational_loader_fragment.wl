(* Reads the exact rational series prepared outside Mathematica. *)

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
