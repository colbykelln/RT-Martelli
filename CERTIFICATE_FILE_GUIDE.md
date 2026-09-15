# Certificate file guide for the write-up

The paper uses the two-generator, seven-relator presentation
`T=<x,y | r1,...,r7>`.  Some scripts use `a,b` internally; in the paper these
are to be read as `x,y`.

## HNN extension and the seven-relator presentation

- Source construction: `mapping_torus_pi1.py` in the user's `RT_stuff`
  directory.
- Raw Regina 7.3 transcript: `T_onto_H.txt` in the user's
  `unfiled terminal output` directory.
- The transcript records the HNN-extension presentation, the Regina
  `intelligentSimplify()` isomorphism, and the resulting seven-relator
  presentation.

Reference these in the paragraph where `pi_1(M)` is replaced by `T`.

## Matrix generators and the surjection `T -> HMat`

- Search and exact verifier: `find_inverse_words_hom15.sage`.
- Machine-readable certificate: `hom15_inverse_matrix_certificate.txt`.
- Copy-ready word identities: `hom15_surjectivity_words.tex`.
- Full copy-ready matrix list and summary of checks:
  `HMat_matrices_and_checks.tex`.
- Generator for the full TeX file: `make_hmat_tex.py`.

Reference `find_inverse_words_hom15.sage` and
`hom15_inverse_matrix_certificate.txt` in the proof of surjectivity.  Insert
`hom15_surjectivity_words.tex` into the paper (or an appendix) for the actual
identities `h_i=W_i(h15,h16)`.

The first sixteen words in the certificate match the corrected depth-10 run
in `T_onto_H.txt`.  All are checked by exact arithmetic over `QQ`.

## Construction of `P` and the identification `P ~= HMat`

- Coxeter factorizations: `coxeter_words_for_T.sage`.
- Index computation and Reidemeister--Schreier construction:
  `automatic_T_via_coxeter.g`.
- Resulting presentation and rewriting data:
  `T_standard_RS_presentation.txt`.
- GAP construction specialized to homomorphism 15: `verify_hom15.g`.
- Exported relators and forward images: `hom15_gap_data.txt`.
- Exact matrix check of `F1`, `F2`, and all six relators:
  `presentation_matrices.sage`.

Reference these where the paper introduces
`P=<F1,F2 | R1,...,R6>` and explains why `P` is isomorphic to `HMat`.

## Explicit inverse and final isomorphism certificate

- Overview and displayed inverse words: `HOM15_ISOMORPHISM_CERTIFICATE.md`.
- GAP/KBMAG verification script: `check_inverse_relators_kbmag.g`.
- Output certificate: `hom15_inverse_kbmag_certificate.txt`.

Reference the GAP script and its output where the paper verifies that the
proposed map `P -> T` respects all six relators and that the composite fixes
`x` and `y`.  The refreshed output has eight successful identity checks.

## Files not to use as primary certificates

- Do not use the files with `5rel` in their names; the paper uses the
  seven-relator presentation from the Regina transcript.
- Do not use `verify_hom15_inverse.py` or
  `hom15_inverse_verification.txt`; that Regina attempt was inconclusive.
- Do not quote the depth-9 run in `T_onto_H.txt`; it did not find `h7` within
  that search bound.  Use the complete depth-10 run and the regenerated exact
  certificate instead.
