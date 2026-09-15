"""Print exact matrices for F1,F2 and the six relators of P."""

import ast
import re
from pathlib import Path
from sage.all import QQ, identity_matrix, matrix

I = identity_matrix(QQ, 6)

def swap(i, j):
    value = matrix(I)
    value.swap_rows(i, j)
    return value

s = [
    swap(0, 1),
    swap(1, 2),
    swap(2, 3),
    matrix(QQ, [
        [ 0,-1,-1, 0, 0, 1],
        [-1, 0,-1, 0, 0, 1],
        [-1,-1, 0, 0, 0, 1],
        [ 0, 0, 0, 1, 0, 0],
        [ 0, 0, 0, 0, 1, 0],
        [-1,-1,-1, 0, 0, 2],
    ]),
    swap(3, 4),
    matrix.diagonal(QQ, [1,1,1,1,-1,1]),
]

# Preimages in Gamma^5 returned by GAP's IsomorphismFpGroup.
p_words = [
    [1,3,4,6,5,6,-2,-3,-5,-3,5,6,5,3,-6,-5,-4,-3,-2,3,4,6,5,6,
     -2,-3,-5,-3,-4,-1,4,3,2,6,5,3,-6,-5,-6,-4,-3,-2,4,3,2,6,5,3,
     -6,-5,-6,-4,-3,-2],
    [2,3,4,6,5,6,-3,-5,-6,-2,-3,-4,2,3,4,6,5,6,-3,-5,-6,-2,-3,-4,
     1,4,3,5,3,2,-6,-5,-6,-4,-3,2,3,4,5,6,-3,-5,-6,-5,-4,-1,4,3,
     2,6,5,3,-6,-5,-6,-4,-3,-2],
]

def coxeter_value(word):
    value = I
    for letter in word:
        # All six Coxeter generators are involutions.
        value *= s[abs(letter)-1]
    return value

F = [coxeter_value(word) for word in p_words]

data = Path("hom15_gap_data.txt").read_text()
relator_text = re.search(r"P_RELATORS:=(.*?);\s*\nIMAGE_A", data, re.S).group(1)
relators = ast.literal_eval(relator_text)

def p_value(word):
    value = I
    for letter in word:
        generator = F[abs(letter)-1]
        value *= generator if letter > 0 else generator.inverse()
    return value

for i, value in enumerate(F, 1):
    print(f"F{i}=")
    print(value)

for i, word in enumerate(relators, 1):
    value = p_value(word)
    print(f"R{i}_letters={word}")
    print(f"R{i}_matrix=")
    print(value)
    assert value == I
