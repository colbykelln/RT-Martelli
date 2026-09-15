"""Find and exactly verify words in hom15(a),hom15(b) for the RT matrices."""

import json
from pathlib import Path
from sage.all import QQ, identity_matrix, matrix

with open("/Users/Colby/Documents/Math/code_4_math_stuff/RT_stuff/ratcliffe-tschantz.json") as stream:
    rt = json.load(stream)

targets = [matrix(QQ, value) for value in rt["pasting_matrices"].values()]

# Coxeter words for the two generators F1,F2 of the certified
# Reidemeister--Schreier presentation of the same subgroup.
p_words = [
    [1,3,4,6,5,6,-2,-3,-5,-3,5,6,5,3,-6,-5,-4,-3,-2,3,4,6,5,6,
     -2,-3,-5,-3,-4,-1,4,3,2,6,5,3,-6,-5,-6,-4,-3,-2,4,3,2,6,5,3,
     -6,-5,-6,-4,-3,-2],
    [2,3,4,6,5,6,-3,-5,-6,-2,-3,-4,2,3,4,6,5,6,-3,-5,-6,-2,-3,-4,
     1,4,3,5,3,2,-6,-5,-6,-4,-3,2,3,4,5,6,-3,-5,-6,-5,-4,-1,4,3,
     2,6,5,3,-6,-5,-6,-4,-3,-2],
]

# Faithful geometric representation of Gamma^5 used in the other scripts.
def swap(i, j):
    value = identity_matrix(QQ, 6)
    value.swap_rows(i, j)
    return value
s1 = swap(0, 1)
s2 = swap(1, 2)
s3 = swap(2, 3)
s4 = matrix(QQ, [[0,-1,-1,0,0,1],[-1,0,-1,0,0,1],
                 [-1,-1,0,0,0,1],[0,0,0,1,0,0],[0,0,0,0,1,0],
                 [-1,-1,-1,0,0,2]])
s5 = swap(3, 4)
s6 = matrix.diagonal(QQ, [1,1,1,1,-1,1])
coxeter = [s1,s2,s3,s4,s5,s6]
for word in p_words:
    value = identity_matrix(QQ, 6)
    for letter in word:
        value *= coxeter[abs(letter)-1]
    targets.append(value)
a = matrix(QQ, [
    [-1,2,-1,-2,0,3], [0,-2,0,1,0,-2], [0,2,-1,-2,-1,3],
    [-1,2,0,-2,-1,3], [0,1,0,-2,0,2], [-1,4,-1,-4,-1,6]])
b = matrix(QQ, [
    [0,1,-2,0,0,2], [-1,0,0,0,0,0], [0,0,-1,-1,0,1],
    [0,0,-1,0,-1,1], [0,2,-2,-1,-1,3], [0,2,-3,-1,-1,4]])

identity = identity_matrix(QQ, 6)
letters = (1, 2, -1, -2)
matrices = {1: a, 2: b, -1: a.inverse(), -2: b.inverse()}

def key(m):
    return tuple(m.list())

def evaluate(word):
    value = identity
    for letter in word:
        value *= matrices[letter]
    return value

def table_to_depth(depth):
    table = {key(identity): (identity, ())}
    frontier = [(identity, ())]
    for length in range(1, depth + 1):
        following = []
        for value, word in frontier:
            for letter in letters:
                if word and letter == -word[-1]:
                    continue
                new_value = value * matrices[letter]
                new_word = word + (letter,)
                k = key(new_value)
                if k not in table:
                    table[k] = (new_value, new_word)
                    following.append((new_value, new_word))
        frontier = following
        print("DEPTH", length, "TABLE_SIZE", len(table))
    return table

SEARCH_DEPTH = 10
table = table_to_depth(SEARCH_DEPTH)

def no_cancel(word, letter):
    return not word or letter != -word[-1]

def meet_in_the_middle(target):
    """Match the corrected search in M5.py: target = left * right."""
    direct = table.get(key(target))
    if direct is not None:
        return direct[1]

    back_level = [(identity, ())]
    back_seen = {key(identity)}
    for _ in range(SEARCH_DEPTH):
        next_back = []
        for right, right_word in back_level:
            for letter in letters:
                if not no_cancel(right_word, letter):
                    continue
                new_right = right * matrices[letter]
                new_right_word = right_word + (letter,)
                right_key = key(new_right)
                if right_key in back_seen:
                    continue
                back_seen.add(right_key)

                # We want left * new_right = target, so
                # left = target * new_right^{-1}.
                found = table.get(key(target * new_right.inverse()))
                if found is not None:
                    return found[1] + new_right_word
                next_back.append((new_right, new_right_word))
        back_level = next_back
    return None

answers = []
for number, target in enumerate(targets, start=1):
    answer = meet_in_the_middle(target)
    print("TARGET", number, "WORD", answer)
    assert answer is not None
    assert evaluate(answer) == target
    answers.append(answer)

Path("hom15_inverse_matrix_certificate.txt").write_text(
    "ALL_18_MATRIX_EQUALITIES=true\n"
    "SEARCH_ALGORITHM=corrected_M5_meet_in_the_middle\n"
    "SEARCH_HALF_DEPTH=10\n"
    "WORD_ENCODING: 1=x, 2=y, -1=x^-1, -2=y^-1\n"
    "The first 16 targets are the natural HMat generators; targets 17 and 18 "
    "are F1 and F2 of the Reidemeister--Schreier presentation.\n"
    + "".join(f"word(h_{i})={answers[i-1]}\n" for i in range(1,17))
    + f"inverse_image(F1)={answers[16]}\n"
    + f"inverse_image(F2)={answers[17]}\n"
)

def tex_word(word):
    if not word:
        return "I_6"
    pieces = []
    start = 0
    while start < len(word):
        stop = start + 1
        while stop < len(word) and word[stop] == word[start]:
            stop += 1
        letter = word[start]
        exponent = (stop - start) * (1 if letter > 0 else -1)
        symbol = "X" if abs(letter) == 1 else "Y"
        pieces.append(symbol if exponent == 1 else f"{symbol}^{{{exponent}}}")
        start = stop
    return " ".join(pieces)

tex_lines = []
for i, answer in enumerate(answers[:16], start=1):
    punctuation = "." if i == 16 else ","
    tex_lines.append(f"h_{{{i}}} &= {tex_word(answer)}{punctuation}\\\\")

Path("hom15_surjectivity_words.tex").write_text(
    "% Generated and exactly verified by find_inverse_words_hom15.sage.\n"
    "% Here X=h_{15}=theta(x) and Y=h_{16}=theta(y).\n"
    "\\begin{align*}\n"
    + "\n".join(tex_lines)
    + "\n\\end{align*}\n"
)
