F := FreeGroup("s1","s2","s3","s4","s5","s6");;
s := GeneratorsOfGroup(F);;
m := [[1,3,2,2,2,2],[3,1,3,2,2,2],[2,3,1,3,3,2],
      [2,2,3,1,2,2],[2,2,3,2,1,4],[2,2,2,2,4,1]];;
rels := List(s, x -> x^2);;
for i in [1..6] do for j in [i+1..6] do Add(rels,(s[i]*s[j])^m[i][j]); od; od;
G := F/rels;; gs := GeneratorsOfGroup(G);;
W := letters -> Product(letters, i -> gs[i]);;
a := W([2,3,4,5,3,2,1,6,5,3,2,4,3,5,6,5,3,4,1,2,3,5,6,2,3,5,1,2,3,2]);;
b := W([2,3,5,3,2,1,4,3,5,6,5,3,2,4,3,6,5,6,4,2,3,5,2,3,1,2]);;
H := Subgroup(G,[a,b]);;
Print("INDEX=",Index(G,H),"\n");
iso := IsomorphismFpGroup(H);;
P := Image(iso);;
Print("P_GENERATORS=",GeneratorsOfGroup(P),"\n");
Print("P_RELATORS=",RelatorsOfFpGroup(P),"\n");
Print("IMAGE_A=",Image(iso,a),"\n");
Print("IMAGE_B=",Image(iso,b),"\n");
Print("IMAGE_A_LETTERS=",LetterRepAssocWord(UnderlyingElement(Image(iso,a))),"\n");
Print("IMAGE_B_LETTERS=",LetterRepAssocWord(UnderlyingElement(Image(iso,b))),"\n");
Print("PREIMAGES_P_GENERATORS=",
      List(GeneratorsOfGroup(P),x->PreImagesRepresentative(iso,x)),"\n");
Print("PREIMAGE_LETTERS=",
      List(GeneratorsOfGroup(P),x->LetterRepAssocWord(
        UnderlyingElement(PreImagesRepresentative(iso,x)))),"\n");
PrintTo("hom15_gap_data.txt",
  "P_RELATORS:=",List(RelatorsOfFpGroup(P),x->LetterRepAssocWord(UnderlyingElement(x))),";\n",
  "IMAGE_A:=",LetterRepAssocWord(UnderlyingElement(Image(iso,a))),";\n",
  "IMAGE_B:=",LetterRepAssocWord(UnderlyingElement(Image(iso,b))),";\n");
QUIT;
