LoadPackage("kbmag");
Read("hom15_gap_data.txt");

FT := FreeGroup("a","b");; tg := GeneratorsOfGroup(FT);; a:=tg[1];; b:=tg[2];;
if IsBound(USE_FIVE_RELATOR_PRESENTATION) and USE_FIVE_RELATOR_PRESENTATION then
trels := [
 a^-3*b^-2*a^3*b^-1*a^-1*b^2*a*b,
 a^-4*b*(b*a)^2*a*b^-1*a^-1*b^-1*a^2*b^-1,
 a*b^-1*a^-1*b^-2*a*(b^2*a^-1)^2*a^-1*b*a*b*a^-3*b,
 (b^-1*a^2)^2*b^-1*a^-1*b^-2*a^-3*(b^2*a)^2*b*a^-2,
 a^2*(a*b^-1)^2*b^-1*a^3*b*a^-2*b^-2*a*(a*b^-2)^2
];;
certificate_file := "hom15_inverse_kbmag_5rel_certificate.txt";;
else
trels := [
 a^-3*b^-2*a^3*b^-1*a^-1*b^2*a*b,
 a^-4*b^2*a*b*a^2*b^-1*a^-1*b^-1*a^2*b^-1,
 a*b^-1*a^-1*b^-2*a*b^2*a^-1*b^2*a^-2*b*a*b*a^-3*b,
 a^3*b^-1*a*b^-2*a^3*b*a^-2*b^-2*a^2*b^-2*a*b^-2,
 a^3*b^-1*a*b^-1*a^-3*b^2*a^-1*b^2*a*b^-1*a^2*b^-1*a^2*b^-1*a^-1*b^-3,
 a^4*b^-1*a*b^-1*a^-3*b^2*a^-2*b^2*a*b*a^-1*b*a^-2*b^-2*a^2*b^-2,
 a^-3*b^2*a*b^2*a*b*a^-1*b*a^-2*b^-2*a^-1*b^2*a*b*a^-1*b^-2*a^-1*b^2*a*b*a^-1*b*a^-2*b^-1
];;
certificate_file := "hom15_inverse_kbmag_certificate.txt";;
fi;
T := FT/trels;; Tgens:=GeneratorsOfGroup(T);;

WordFromLetters := function(list, gens)
  local value, letter;
  value := One(Parent(gens[1]));
  for letter in list do
    if letter > 0 then value := value*gens[letter];
    else value := value*gens[-letter]^-1; fi;
  od;
  return value;
end;;

invF1 := [-1,-1,-1,2,2,-1,-1,2,1,1,-2,-1,-1,-1,2,1];;
invF2 := [2,-1,-1,-2,-2,-1,2,2,1,2,1,1,-2,1,-2];;
invimgs := [WordFromLetters(invF1,Tgens),WordFromLetters(invF2,Tgens)];;
tests := List(P_RELATORS, r -> WordFromLetters(r,invimgs));;
Add(tests, WordFromLetters(IMAGE_A,invimgs)*Tgens[1]^-1);
Add(tests, WordFromLetters(IMAGE_B,invimgs)*Tgens[2]^-1);

LogTo(certificate_file);
R := KBMAGRewritingSystem(T);;
opts := OptionsRecordOfKBMAGRewritingSystem(R);;
opts.maxeqns := 100000;; opts.tidyint := 200;;
SetInfoLevel(InfoRWS,1);
ok := KnuthBendix(R);;
Print("KB_COMPLETE=",ok," RULES=",Length(Rules(R)),"\n");
for i in [1..Length(tests)] do
  red := ReducedForm(R,UnderlyingElement(tests[i]));
  Print("TEST_",i,"_IDENTITY=",red=One(FT)," REDUCED_LENGTH=",Length(red),"\n");
od;
LogTo();
QUIT;
