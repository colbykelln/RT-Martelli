# Isomorphism certificate for homomorphism 15

Let `T` be the two-generator, seven-relator group in `all_in_one.txt`, with
generators `x,y`.  Let `P=<F1,F2 | R1,...,R6>` be the GAP
Reidemeister--Schreier presentation of the matrix subgroup `HMat` inside the
Coxeter group `Gamma^5`.  The presentation and the images of `x,y` under
homomorphism 15 are generated in `hom15_gap_data.txt` (where the exported
variables retain the internal names `IMAGE_A` and `IMAGE_B`).

After identifying `P` with `HMat`, an inverse `P -> T` to the transported
homomorphism is:

```
F1 |-> x^-1 x^-1 x^-1 y y x^-1 x^-1 y x x y^-1 x^-1 x^-1 x^-1 y x
F2 |-> y x^-1 x^-1 y^-1 y^-1 x^-1 y y x y x x y^-1 x y^-1
```

Equivalently, in exponent notation,

```
F1 |-> x^-3 y^2 x^-2 y x^2 y^-1 x^-3 y x
F2 |-> y x^-2 y^-2 x^-1 y^2 x y x^2 y^-1 x y^-1.
```

There are three independent parts to the certificate.

1. `verify_hom15.g` constructs `P` directly from the index-3840 subgroup of
   the Coxeter group and exports its six relators and the forward images.
2. `find_inverse_words_hom15.sage` finds the two inverse words above and
   words for all 16 natural side-pairing matrices using the corrected
   meet-in-the-middle reconstruction from `M5.py`; it verifies every equality
   by exact multiplication in `GL(6,QQ)`.
3. `check_inverse_relators_kbmag.g` substitutes the proposed inverse into all
   six relators of `P`, and also forms
   `inverse(hom15(x))*x^-1` and `inverse(hom15(y))*y^-1`.  All eight words
   reduce to the empty word using oriented consequences of the seven defining
   relators of `T`.

The KBMAG system used for the last step need not be confluent: reduction of a
word to the empty word is still a valid derivation from the defining
relations.  Confluence would only be needed to infer nontriviality from a
nonempty irreducible word.

Thus the displayed assignment defines a homomorphism `P -> T`, and its
composite after homomorphism 15 fixes `x` and `y`.  Since homomorphism 15 was
already shown to be surjective, it is injective and hence an isomorphism.

Reproduce the checks with:

```
gap -q verify_hom15.g
DOT_SAGE=/tmp/codex_sage_hom15 sage find_inverse_words_hom15.sage
gap -q check_inverse_relators_kbmag.g
```
