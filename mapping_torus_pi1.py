"""
Helper function to get HNN extension written by Claude, some version from early 2026
Fundamental group of the mapping torus of φ with fiber F.

From Theorem 9, π₁(F) has generators a_i, b_i (i = 1,...,6, mod 6) and relations:
  R1: a_{i+2} = a_i a_{i+1}
  R2: b_{i+2} = b_i b_{i+1}
  R3: a_i⁻¹ b_{i+1} a_{i+2} = b_i⁻¹ a_{i+1} b_{i+2}

The automorphism φ_* acts as:
  a1 -> a3 b5 a4⁻¹ b3⁻¹ a1⁻¹
  a2 -> a3 b3⁻¹ a1⁻¹
  a3 -> a1 b3 a2⁻¹ b1⁻¹ a5⁻¹
  a4 -> a1 b1⁻¹ a5⁻¹
  a5 -> a5 b1 a6⁻¹ b5⁻¹ a3⁻¹
  a6 -> a5 b5⁻¹ a3⁻¹
  b1 -> a3⁻¹
  b2 -> a3 a1⁻¹
  b3 -> a1⁻¹
  b4 -> a1 a5⁻¹
  b5 -> a5⁻¹
  b6 -> a5 a3⁻¹

The mapping torus M_φ has:
  π₁(M_φ) = π₁(F) ⋊_φ* ℤ
           = < a_i, b_i, t | (relations of π₁(F)),
                              t a_i t⁻¹ = φ_*(a_i),
                              t b_i t⁻¹ = φ_*(b_i) >

Generator index convention (0-based for Regina):
  a1=0, a2=1, a3=2, a4=3, a5=4, a6=5,
  b1=6, b2=7, b3=8, b4=9, b5=10, b6=11,
  t=12
"""

import regina


def idx(name):
    """Return 0-based generator index by name."""
    table = {
        'a1': 0, 'a2': 1, 'a3': 2, 'a4': 3, 'a5': 4, 'a6': 5,
        'b1': 6, 'b2': 7, 'b3': 8, 'b4': 9, 'b5': 10, 'b6': 11,
        't': 12,
    }
    return table[name]


def word(*letters):
    """
    Build a Regina GroupExpression from a sequence of (generator_index, exponent) pairs.
    e.g. word((0, 1), (2, -1)) represents a1 * a3⁻¹
    """
    expr = regina.GroupExpression()
    for (gen, exp) in letters:
        expr.addTermFirst(gen, exp)  # addTermFirst prepends; build in reverse
    # Actually addTermFirst adds at front, so build list then add in order
    return expr


def make_word(pairs):
    """
    Build a GroupExpression from a list of (gen_index, exponent) pairs,
    in left-to-right order.
    """
    expr = regina.GroupExpression()
    for (gen, exp) in pairs:
        expr.addTermLast(gen, exp)
    return expr


def build_pi1_mapping_torus():
    """
    Construct and return a Regina GroupPresentation for π₁(M_φ).

    Generators (13 total):
      indices 0-5:  a1,...,a6
      indices 6-11: b1,...,b6
      index  12:    t  (stable letter of the HNN extension)

    Relations:
      (A) Fiber relations from π₁(F):
          R1_i: a_{i+2} a_{i+1}⁻¹ a_i⁻¹ = 1   (a_{i+2} = a_i a_{i+1})
          R2_i: b_{i+2} b_{i+1}⁻¹ b_i⁻¹ = 1   (b_{i+2} = b_i b_{i+1})
          R3_i: a_i⁻¹ b_{i+1} a_{i+2} b_{i+2}⁻¹ a_{i+1}⁻¹ b_i = 1
                (from a_i⁻¹ b_{i+1} a_{i+2} = b_i⁻¹ a_{i+1} b_{i+2})
      (B) HNN relations:
          t a_i t⁻¹ φ_*(a_i)⁻¹ = 1  for i=1,...,6
          t b_i t⁻¹ φ_*(b_i)⁻¹ = 1  for i=1,...,6
    """

    G = regina.GroupPresentation(13)  # 13 generators

    # ------------------------------------------------------------------
    # Helper: indices for a_i and b_i (1-indexed, mod 6 with base 1)
    # ------------------------------------------------------------------
    def a(i):
        return (i - 1) % 6          # 0-based index

    def b(i):
        return 6 + (i - 1) % 6     # 0-based index

    T = 12  # index for t

    # ------------------------------------------------------------------
    # (A) Fiber relations, i = 1,...,6 (mod 6)
    # ------------------------------------------------------------------
    for i in range(1, 7):
        i1 = i + 1  # i+1 mod 6 (use a/b helpers which do mod)
        i2 = i + 2

        # R1: a_{i+2} = a_i a_{i+1}  =>  a_{i+2} a_{i+1}^{-1} a_i^{-1} = 1
        G.addRelation(make_word([(a(i2), 1), (a(i1), -1), (a(i), -1)]))

        # R2: b_{i+2} = b_i b_{i+1}  =>  b_{i+2} b_{i+1}^{-1} b_i^{-1} = 1
        G.addRelation(make_word([(b(i2), 1), (b(i1), -1), (b(i), -1)]))

        # R3: a_i^{-1} b_{i+1} a_{i+2} = b_i^{-1} a_{i+1} b_{i+2}
        #  => a_i^{-1} b_{i+1} a_{i+2} b_{i+2}^{-1} a_{i+1}^{-1} b_i = 1
        G.addRelation(make_word([
            (a(i),  -1),
            (b(i1),  1),
            (a(i2),  1),
            (b(i2), -1),
            (a(i1), -1),
            (b(i),   1),
        ]))

    # ------------------------------------------------------------------
    # (B) HNN relations:  t x t^{-1} = φ_*(x)
    #     written as:  t x t^{-1} φ_*(x)^{-1} = 1
    #
    # φ_* images (as lists of (gen_index, exp) pairs):
    # ------------------------------------------------------------------
    phi = {
        # a_i
        a(1): [(a(3), 1), (b(5), 1), (a(4), -1), (b(3), -1), (a(1), -1)],
        a(2): [(a(3), 1), (b(3), -1), (a(1), -1)],
        a(3): [(a(1), 1), (b(3), 1), (a(2), -1), (b(1), -1), (a(5), -1)],
        a(4): [(a(1), 1), (b(1), -1), (a(5), -1)],
        a(5): [(a(5), 1), (b(1), 1), (a(6), -1), (b(5), -1), (a(3), -1)],
        a(6): [(a(5), 1), (b(5), -1), (a(3), -1)],
        # b_i
        b(1): [(a(3), -1)],
        b(2): [(a(3), 1), (a(1), -1)],
        b(3): [(a(1), -1)],
        b(4): [(a(1), 1), (a(5), -1)],
        b(5): [(a(5), -1)],
        b(6): [(a(5), 1), (a(3), -1)],
    }

    for gen_idx, image_pairs in phi.items():
        # t gen t^{-1} image^{-1} = 1
        rel = [(T, 1), (gen_idx, 1), (T, -1)]
        # append inverse of image (reverse list, negate exponents)
        rel += [(g, -e) for (g, e) in reversed(image_pairs)]
        G.addRelation(make_word(rel))

    return G


def main():
    G = build_pi1_mapping_torus()

    print("=" * 60)
    print("π₁(mapping torus of φ) — raw presentation")
    print("=" * 60)
    print(f"Number of generators: {G.countGenerators()}")
    print(f"Number of relations:  {G.countRelations()}")
    print()
    print("Generator names (0-based index → name):")
    names = ['a1','a2','a3','a4','a5','a6',
             'b1','b2','b3','b4','b5','b6','t']
    for i, n in enumerate(names):
        print(f"  g{i} = {n}")
    print()
    print("Relations (as printed by Regina):")
    for i in range(G.countRelations()):
        print(f"  {G.relation(i)}")

    print()
    print("=" * 60)
    print("Simplifying with Regina...")
    print("=" * 60)
    G.intelligentSimplify()
    print(f"After simplification:")
    print(f"  Generators: {G.countGenerators()}")
    print(f"  Relations:  {G.countRelations()}")
    print()
    for i in range(G.countRelations()):
        print(f"  {G.relation(i)}")

    return G


if __name__ == "__main__":
    G = main()
