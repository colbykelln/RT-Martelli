LoadPackage("kbmag");

# Use a local autcos wrapper that raises gpmakefsa's otherwise hard-coded
# word-difference ceiling.  All binaries besides that wrapper are symlinks to
# the installed KBMAG executables.
_KBExtDir := Directory(Filename(DirectoryCurrent(), "kbmag_custom"));;

# Coxeter presentation of the integral Lorentzian reflection group Gamma^5.
F := FreeGroup("s1","s2","s3","s4","s5","s6");;
s := GeneratorsOfGroup(F);;
m := [
  [1,3,2,2,2,2],
  [3,1,3,2,2,2],
  [2,3,1,3,3,2],
  [2,2,3,1,2,2],
  [2,2,3,2,1,4],
  [2,2,2,2,4,1]
];;
rels := List(s, x -> x^2);;
for i in [1..6] do
  for j in [i+1..6] do Add(rels, (s[i]*s[j])^m[i][j]); od;
od;
W := letters -> Product(letters, i -> s[i]);;

# The 16 RT side pairings, recovered and matrix-verified by
# coxeter_words_for_T.sage.  One pairing is the identity and is omitted.
hwords := [
 [2,3,2,1,4,3,2,5,3,2,1,6,5,3,2,1],
 [2,1,3,2,1,4,3,5,3,2,1,4,3,2,6,5,3,1],
 [2,3,2,1,4,3,5,3,2,4,3,5,6,5],
 [2,3,5,3,2,1,4,6,5,3,2,1,4,3,2,6],
 [2,3,2,1,5,3,2,1,4,3,5,6,5,3,4,6,2,3,5,3,4,2,3,1,2,1],
 [2,3,4,3,2,1,5,3,2,4,3,5,6,5,3,4,1,2,3,5,4,3,2,1],
 [3,2,1,5,3,2,1,6,5,3,2,1,6,5,3,4,2,3],
 [2,1,3,5,3,2,1,6,5,3,2,1,4,3,5,6,2,3,5,3,4,3,2,1],
 [2,3,2,4,5,3,6,5,3,2,1,6,5,3,2,4,3,5,2,3],
 [2,1,3,2,4,5,3,4,6,5,3,2,6,5,3,4,3,2],
 [2,3,5,6,5,3,2,4,3,5,6,5,3,4],
 [2,4,3,5,3,4,6,5,3,2,1,6,5,3,2,6],
 [2,1,3,5,3,2,4,3,6,5,3,2,4,3,6,5,6,5,4,3,1,2],
 [2,3,4,5,3,2,1,6,5,3,2,4,3,5,6,5,3,4,1,2,3,5,6,2,3,5,1,2,3,2],
 [2,3,5,3,2,1,4,3,5,6,5,3,2,4,3,6,5,6,4,2,3,5,2,3,1,2]
];;

G := F/rels;;
gs := GeneratorsOfGroup(G);;
WG := letters -> Product(letters, i -> gs[i]);;
hgens := List(hwords, WG);;

if IsBound(ONLY_KB_SUBGROUP) and ONLY_KB_SUBGROUP = true then
  Rsub := KBMAGRewritingSystem(G);;
  rsubopts := OptionsRecordOfKBMAGRewritingSystem(Rsub);;
  rsubopts.maxeqns := 2000000;; rsubopts.tidyint := 500;;
  Ssub := SubgroupOfKBMAGRewritingSystem(Rsub, hgens);;
  SetInfoLevel(InfoRWS, 1);
  subok := KnuthBendixOnCosetsWithSubgroupRewritingSystem(Rsub,Ssub);;
  Print("SUBGROUP_REWRITING_SYSTEM_SUCCESS: ", subok, "\n");
  if subok then
    SS := RewritingSystemOfSubgroupOfKBMAGRewritingSystem(Rsub,Ssub);;
    Print("SUBGROUP_RULES: ", Length(Rules(SS)), "\n");
    PrintTo("T_subgroup_rewriting_system_summary.txt",
      "Alphabet: ", SS!.alphabet, "\n",
      "Number of rules: ", Length(Rules(SS)), "\n");
    if IsBound(CONTINUE_SUBGROUP_KB) and CONTINUE_SUBGROUP_KB = true then
      ssopts := OptionsRecordOfKBMAGRewritingSystem(SS);;
      ssopts.maxeqns := 2000000;; ssopts.maxstates := 2000000;;
      ssopts.tidyint := 500;;
      ssconf := KnuthBendix(SS);;
      Print("T_SUBGROUP_CONFLUENT: ", ssconf, "\n");
      Print("T_SUBGROUP_FINAL_RULES: ", Length(Rules(SS)), "\n");
    fi;
    if IsBound(WTLEX_SUBGROUP) and WTLEX_SUBGROUP = true then
      ResetRWS(SS);
      weights := Concatenation(List(hwords,w->[Length(w),Length(w)]));;
      SetOrderingOfKBMAGRewritingSystem(SS,"wtlex",weights);
      wtok := AutomaticStructure(SS);;
      Print("T_WTLEX_AUTOMATIC_STRUCTURE_SUCCESS: ",wtok,"\n");
      if wtok then
        Print("T_WTLEX_WORD_ACCEPTOR_STATES: ",
          NumberOfStatesFSA(WordAcceptor(SS)),"\n");
        Print("T_WTLEX_GENERAL_MULTIPLIER_STATES: ",
          NumberOfStatesFSA(GeneralMultiplier(SS)),"\n");
        LogTo("T_wtlex_side_pairing_automatic_structure.txt");
        Print("word_acceptor :=\n"); WriteFSA(WordAcceptor(SS));
        Print("general_multiplier :=\n"); WriteFSA(GeneralMultiplier(SS));
        Print("first_word_difference_automaton :=\n");
        WriteFSA(FirstWordDifferenceAutomaton(SS));
        Print("second_word_difference_automaton :=\n");
        WriteFSA(SecondWordDifferenceAutomaton(SS));
        LogTo();
      fi;
    fi;
    if IsBound(AUTO_AFTER_SUBGROUP) and AUTO_AFTER_SUBGROUP = true then
      # AutomaticStructure must begin from the extracted defining equations,
      # rather than from the preceding (incomplete) KB run.
      ResetRWS(SS);
      ssopts := OptionsRecordOfKBMAGRewritingSystem(SS);;
      ssopts.maxeqns := 2000000;; ssopts.maxwdiffs := 500000;;
      ssopts.tidyint := 500;;
      autsubok := AutomaticStructure(SS, false, false, false);;
      Print("T_AUTOMATIC_STRUCTURE_SUCCESS: ", autsubok, "\n");
      if autsubok then
        Print("T_WORD_ACCEPTOR_STATES: ", NumberOfStatesFSA(WordAcceptor(SS)), "\n");
        Print("T_GENERAL_MULTIPLIER_STATES: ",
          NumberOfStatesFSA(GeneralMultiplier(SS)), "\n");
        LogTo("T_word_acceptor.gap"); WriteFSA(WordAcceptor(SS)); LogTo();
        LogTo("T_general_multiplier.gap"); WriteFSA(GeneralMultiplier(SS)); LogTo();
      fi;
    fi;
  fi;
elif IsBound(ONLY_RS) and ONLY_RS = true then
  H := Subgroup(G, hgens);;
  Print("STANDARD_GAP_INDEX: ", Index(G,H), "\n");
  isoH := IsomorphismFpGroup(H);;
  HP := Image(isoH);;
  Print("RAW_RS_GENERATORS: ", Length(GeneratorsOfGroup(HP)), "\n");
  Print("RAW_RS_RELATORS: ", Length(RelatorsOfFpGroup(HP)), "\n");
  HPS := SimplifiedFpGroup(HP);;
  PrintTo("T_standard_RS_presentation.txt",
    "Group: ", HP, "\nGenerators: ", GeneratorsOfGroup(HP),
    "\nRelators: ", RelatorsOfFpGroup(HP),
    "\nGenerators as Coxeter words: ",
    List(GeneratorsOfGroup(HP), x -> PreImagesRepresentative(isoH, x)),
    "\nSide-pairing generators in this presentation: ",
    List(hgens, x -> Image(isoH, x)), "\n");
  Print("SIMPLIFIED_RS_GENERATORS: ", Length(GeneratorsOfGroup(HPS)), "\n");
  Print("SIMPLIFIED_RS_RELATORS: ", Length(RelatorsOfFpGroup(HPS)), "\n");
else
R := KBMAGRewritingSystem(G);;
opts := OptionsRecordOfKBMAGRewritingSystem(R);;
opts.maxeqns := 2000000;; opts.maxwdiffs := 500000;; opts.tidyint := 500;;
S := SubgroupOfKBMAGRewritingSystem(R, hgens);;
SetInfoLevel(InfoRWS, 1);
Print("Computing the index-3840 RT subgroup inside Gamma^5...\n");
ok := AutomaticStructureOnCosetsWithSubgroupPresentation(
  R, S, false, false, false
);;
Print("COSET_AUTOMATIC_STRUCTURE_SUCCESS: ", ok, "\n");
if ok then
  Print("INDEX: ", Index(R,S), "\n");
  HP := PresentationOfSubgroupOfKBMAGRewritingSystem(R,S);;
  PrintTo("T_coxeter_subgroup_presentation.txt", HP, "\n");
  Print("SUBGROUP_PRESENTATION_GENERATORS: ", Length(GeneratorsOfGroup(HP)), "\n");
  Print("SUBGROUP_PRESENTATION_RELATORS: ", Length(RelatorsOfFpGroup(HP)), "\n");
fi;
fi;
QUIT;
