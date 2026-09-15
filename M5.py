import itertools, regina
from tessellation import Tessellated_manifold, Tessellated_manifold_isometry, Tessellated_manifold_isometry_group
from regina import Perm4
from P5 import P5, P5_iso, I, J, K, rt
from sage.all import QuaternionAlgebra, Matrix, MatrixGroup, QQ
from collections import deque


def get_S5():
    """
    Returns a 5-dimensional sphere, using two copies of P5 glued via identity.
    """
    def pasting_map_sphere(facet):
        p = facet.pol
        manifold = p.manifold
        id = P5_iso()
        return (manifold.polytopes[1-p.index], id)

    return Tessellated_manifold(P5, [0,1], pasting_map_sphere)

def get_M5():
    """
    Returns the manifold obtained with 2^8 copies of P5, induced by the coloring of P5 with 8 colors.
    """
    def facet_state(facet):
        state = True
        if facet.label in {-1, -I, -J, -K}:
            label = -facet.label
            state = not state
        else:
            label = facet.label

        i = [1, I, J, K].index(label)
        if (facet.pol.index[i] + facet.pol.index[i+4]) % 2 == 0:
            return state
        else:
            return not state

    def pasting_map_M5(facet):
        p = facet.pol
        manifold = p.manifold
        x = list(p.index)
        id = P5_iso()
        x[facet.color] = 1 - x[facet.color]
        return (manifold.polytopes[tuple(x)], id)

    return Tessellated_manifold(P5, [x for x in itertools.product(*[[0,1]]*8)], pasting_map_M5, facet_state)

# Action of Q_8 on P5
q8_iso = {1: P5_iso(), -1: P5_iso([-1, -1, -1, -1]), I: P5_iso([1, -1, 1, -1], Perm4(3, 2, 1, 0)), J: P5_iso([1, -1, -1, 1], Perm4(1, 0, 3, 2))}
q8_iso[K] = q8_iso[I] * q8_iso[J]
for i in [I, J, K]:
    q8_iso[-i] = q8_iso[i].inverse()

# Isometry which inverts the small cusps of P5
quotienting_iso = P5_iso(refl=[1, -1, -1, -1], perm=Perm4(0,2,1,3))

def get_minimal_quotient():
    """
    Constructs a small manifold gluing two copies of P5, using Q8.
    The result is isomorphic to the manifold constructed by Ratcliffe and Tschantz.
    """
    def facet_state(facet):
        if facet.pol.index[0] % 2 == 0:
            return facet.label in {1, I, J, K}
        else:
            return facet.label in {-1, I, J, K}

    def pasting_map_minimal_quotient(facet):
        p = facet.pol
        manifold = p.manifold
        level, rotation = p.index

        if level % 2 == 0:
            if level==2 and facet.state:
                return None
            iso = q8_iso[{1: 1, -1: 1, I: J, -I: J, J: K, -J: K, K: I, -K: I}[facet.label]]
            new_p = manifold.polytopes[(level + (1 if facet.state else -1), rotation * (I if facet.index[4] == 1 else -I))]
            return (new_p, iso)

        else:
            if level==-1 and not facet.state:
                return None
            iso = q8_iso[{1: 1, -1: 1, I: -K, -I: -K, J: -I, -J: -I, K: -J, -K: -J}[facet.label]]
            new_p = manifold.polytopes[(level + (1 if facet.state else -1), rotation * (-I if facet.index[4] == 1 else I))]
            return (new_p, iso)

    mnf = Tessellated_manifold(P5, [(-1, I), (-1, -I), (0, 1), (0, -1), (1, I), (1, -I), (2, 1), (2, -1)], pasting_map_minimal_quotient, facet_state)
    isom = Tessellated_manifold_isometry(mnf, mnf.polytopes[(0,1)], mnf.polytopes[(1, I)], quotienting_iso)
    return mnf.get_quotient(Tessellated_manifold_isometry_group(isom))


def get_M5_cyclic_covering():
    """
    Returns the cyclic covering of M5, given by the choice of a state on
    each facet of P5.
    """
    def facet_state(facet):
        # Sets the facet's state
        res = facet.label in {1, I, J, K}
        if facet.pol.index[facet.color] == 1:
            res = not res
        if facet.pol.index[(facet.color + 4) % 8] == 1:
            res = not res
        return res

    def pasting_map_leveled_M5(facet):
        p = facet.pol
        manifold = p.manifold
        x = list(p.index[0:8])
        level = p.index[8]

        new_level = level + (1 if facet.state else -1)

        if new_level not in {0, 1, 2}:
            return None

        id = P5_iso()
        x[facet.color] = 1 - x[facet.color]
        return (manifold.polytopes[tuple(x) + (new_level, )], id)

    return Tessellated_manifold(P5, [x for x in itertools.product(*[[0,1]]*8 + [[1, 0, 2]]) if sum(x) % 2 == 0], pasting_map_leveled_M5, facet_state)

def get_M5_two_quotient():
    """
    Returns the quotient of the cyclic covering of M5 by the action of Q8 and the
    additional isometry which inverts the two small cusps.
    The result is isomorphic to the manifold constructed by Ratcliffe and Tschantz.
    """
    leveled_M5 = get_M5_cyclic_covering()
    isometries = [
        Tessellated_manifold_isometry(leveled_M5, leveled_M5.polytopes[(0,)*9], leveled_M5.polytopes[x], y) for x, y in [
            ((1,1,0,0,1,1,0,0,0), P5_iso()),
            ((1,0,1,0,1,0,1,0,0), P5_iso()),
            ((1,0,0,1,1,0,0,1,0), P5_iso()),
            ((1,0,1,0,0,0,0,0,0), q8_iso[I]),
            ((1,0,0,1,0,0,0,0,0), q8_iso[J]),
            ((1,0,0,0,0,0,0,0,1), quotienting_iso)
        ]
    ]
    G = Tessellated_manifold_isometry_group(*isometries, iterations=3)
    return leveled_M5.get_quotient(G)

def get_M5_other_two_quotient():
    """
    Returns another quotient of the cyclic covering of M5 by the action of Q8 and the
    additional isometry which inverts the two small cusps.
    We take a different index two subgroup.
    The result is NOT isomorphic to the manifold constructed by Ratcliffe and Tschantz.
    """
    leveled_M5 = get_M5_cyclic_covering()
    isometries = [
        Tessellated_manifold_isometry(leveled_M5, leveled_M5.polytopes[(0,)*9], leveled_M5.polytopes[x], y) for x, y in [
            ((1,1,0,0,1,1,0,0,0), P5_iso()),
            ((1,0,1,0,1,0,1,0,0), P5_iso()),
            ((1,0,0,1,1,0,0,1,0), P5_iso()),
            ((1,0,1,0,0,0,0,0,0), q8_iso[I]),
            ((1,0,0,0,0,1,0,0,0), q8_iso[K]),
            ((1,0,0,0,0,0,0,0,1), quotienting_iso)
        ]
    ]
    G = Tessellated_manifold_isometry_group(*isometries, iterations=3)
    return leveled_M5.get_quotient(G)



def get_rt_minimal_manifold():
    """
    Constructs the manifold described by Ratcliffe and Tschantz.
    """
    def pasting_map_rt(facet):
        p = facet.pol
        number = facet.number + 16*p.index
        target_number = rt["mappings"][number]
        def get_matrix(i,j):
            standard_reflection = p.facets[(-1, 1, -1, -1, -1)].reflection_matrix
            return p.facet_from_number(j % 16).reflection_matrix * \
                (standard_reflection if j >= 16 else Matrix.identity(6)) * \
                Matrix(rt["pasting_matrices"][str(i)]) * \
                (standard_reflection if i >= 16 else Matrix.identity(6))

        m = get_matrix(number, target_number) if str(number) in rt["pasting_matrices"] else get_matrix(target_number, number).inverse()
        try:
            iso = P5_iso.from_lorentzian_matrix(m)
        except:
            raise Exception("Cannot paste {} to {}".format(number, target_number))
        return (p.manifold.polytopes[1 if target_number >= 16 else 0], iso)

    return Tessellated_manifold(P5, [0,1], pasting_map_rt)


## OTHER APPROACH
## in sage in terminal inline this function is equivalent to GM = MatrixGroup([Matrix(m) for m in rt["pasting_matrices"].values()]); GM
## having run M5.py in sage -python 
def get_pi1_rt_mat_gp():
    gens = []
    for m in rt["pasting_matrices"].values():
        # a check that the matrix is indeed in SO(5,1)
        assert (Matrix(m).transpose() * Matrix.diagonal([1, 1, 1, 1, 1, -1]) * Matrix(m)) == Matrix.diagonal([1, 1, 1, 1, 1, -1]), m
        # matrices are actually integer matrices
        gens = gens + [Matrix(QQ,m)]
    G = MatrixGroup(gens)
    return G

identity = Matrix([[1, 0, 0, 0, 0, 0],[0, 1, 0, 0, 0, 0],[0, 0, 1, 0, 0, 0],[0, 0, 0, 1, 0, 0],[0, 0, 0, 0, 1, 0],[0, 0, 0, 0, 0, 1]])

### help from CHATGPT for below functions 


'''
def generate_distinct_matrix_words(M_gens, max_len):
    """
    Generate all *distinct* matrix elements in the subgroup of GL(n, R)
    generated by M_gens, for words up to length max_len.
    """
    seen = {str(identity)}
    current_level = [identity]

    gen_list = list(M_gens) + [g.inverse() for g in M_gens]

    for _ in range(max_len):
        next_level = []
        for m in current_level:
            for g in gen_list:
                new_mat = m * g
                if str(new_mat) not in seen:
                    seen.add(str(new_mat))
                    next_level.append(new_mat)
        current_level = next_level

    return list(seen)

def generate_distinct_matrix_words(M_gens, max_len):
    """
    Generate all *distinct* matrix elements in the subgroup of GL(n, R)
    generated by M_gens, for words up to length max_len.
    """
    seen = [identity]
    seen_set = set()  # for fast str-based deduplication
    seen_set.add(str(identity))

    current_level = [identity]
    gen_list = list(M_gens) + [g.inverse() for g in M_gens]

    for _ in range(max_len):
        next_level = []
        for m in current_level:
            for g in gen_list:
                new_mat = m * g
                key = str(new_mat)
                if key not in seen_set:
                    seen_set.add(key)
                    seen.append(new_mat)
                    next_level.append(new_mat)
        current_level = next_level

    return seen
'''
def generate_distinct_matrix_words_keep_track(M_gens, max_len):
    """
    Returns a list of (matrix, word) pairs where:
    - matrix is the resulting matrix product
    - word is a tuple of generator indices (+ for gens, - for inverse gens)
    """
    tup_ID = tuple(tuple(QQ(x) for x in row) for row in identity.rows())
    seen = {tup_ID: (identity, ())}
    current_level = [(identity, ())]

    gen_list = list(M_gens) + [g.inverse() for g in M_gens]
    gen_indices = list(range(len(M_gens))) + list(range(-len(M_gens), 0))


    yield (identity, ())

    for _ in range(max_len):
        next_level = []
        for mat, word in current_level:
            for g, idx in zip(gen_list, gen_indices):
                new_mat = mat * g
                new_word = word + (idx,)
                # continuing if we get identity matrix with nontrivial (unreduced) word
                if new_mat == identity and new_word != ():
                    continue
                key = tuple(tuple(QQ(x) for x in row) for row in new_mat.rows())
                if key not in seen:
                    seen[key] = (new_mat, new_word)
                    yield (new_mat, new_word)
                    next_level.append((new_mat, new_word))
        current_level = next_level
'''
### Trying 'yeild' to make all_homs look at each word one at a time 
### rather than having generate_distinct finish and return the whole list
    for _ in range(max_len):
        next_level = []
        for mat, word in current_level:
            for g, idx in zip(gen_list, gen_indices):
                new_mat = mat * g
                if new_mat == identity and word != ():
                    continue
                key = tuple(tuple(QQ(entry) for entry in row) for row in new_mat.rows())
                if key not in seen:
                    new_word = word + (idx,)
                    seen[key] = (new_mat, new_word)
                    next_level.append((new_mat, new_word))
        current_level = next_level

    return list(seen.values())
'''
def all_homomorphisms_to_matrix_group(G_pres, Mat, max_word_len):
    """
    Compute all homomorphisms from G_pres to a matrix group generated by M_gens,
    allowing images to be words of length <= max_word_len in M_gens.
    """
    # casting G_pres as a sage group. was regina before if gotten from triangulation.fundamentalgroup e.g. M = get_minimal_quotient(); X = 
    # M.triangulate(); G = X.fundamentalGroup();
    G_pres = G_pres.sage()
    
    gens = list(G_pres.gens())
    rels = list(G_pres.relations())

    gens_str = [str(gen) for gen in gens]

    print("generator key = ",{**{i: gens_str[i] for i in range(len(gens_str))}, **{i - len(gens_str): gens_str[i] + '^-1' for i in range(len(gens_str))}})

    M_gens = Mat.gens()

    # Enumerate all candidate matrix group elements of word length up to max_word_len
    M_word_pairs = list(generate_distinct_matrix_words_keep_track(M_gens, max_word_len))

    ## Maybe I want to modify this so that not only do I track what the products of the matrices are from generate_distinct_matrix_words but also 
    ## what the words (in the generators of the matrix gp M) actually are. 
    homs = []
    count = 0

    for assignment in itertools.product(M_word_pairs, repeat=len(gens)):
        gen_to_mat = dict(zip(gens_str, assignment))  # maps "x0" → (matrix, word)
        count += 1

        if all(check_relation(r, gen_to_mat, identity) for r in rels):
            homs.append(gen_to_mat)
            print("***Found one! Count:", count)
        elif count % 100000 == 0:
            print("Checked:", count)

    return homs

## be careful about parent type of M. need it to be a matrix group, living inside of GL(6, ZZ), and for sage to recognize it as a group. 
## in order to use .parent().identity, and .is_one(). 
## currently only have it as a matrix group in libgap in sage

## CHeCK rewrite below: rel is a relator in G_pres.sage(). 

def check_relation(rel, gen_to_mat, id_mat):
    """
    Check whether the matrix product corresponding to a relation evaluates to the identity.
    rel is an element of a finitely presented group (not iterable).
    """
    mat = id_mat
    word = rel.syllables() #https://doc.sagemath.org/html/en/reference/groups/sage/groups/free_group.html#sage.groups.free_group.FreeGroupElement.syllables 
    # looks like ((a, 2), (b, -1), (a, 3)) if word is a^2 * b^-1 * a^3.
    
    for gen, exp in word:
        def matrix_from_string(s):
            """
            Convert a multiline string of rows into a Sage matrix.
            The string should have rows separated by newlines, and entries separated by spaces.
            Brackets around each row are optional.
        
            Example:
            '[1 0 0]\n[0 1 0]\n[0 0 1]'
            or
            '1 0 0\n0 1 0\n0 0 1'
            """
            rows = []
            for line in s.strip().split('\n'):
                line = line.strip().strip('[]')  # remove whitespace and brackets
                row = [QQ(x) for x in line.split()]
                rows.append(row)
            return Matrix(rows)
        #m = matrix_from_string(gen_to_mat[str(gen)])
        #m = matrix_from_string(gen_to_mat[str(gen)][0])
        m = gen_to_mat[str(gen)][0]
        '''
        #some formatting fixes 
        m = m.replace('\n', ',')
        m = '['+m+']'
        print("m = ", m)
        mat = Matrix(m)
        '''
        if exp < 0:
            m = m.inverse()
            exp = -exp
        mat *= m**exp

    return mat.is_one()

# Checks if matrix group G2 is generated by matrix group G1
def check_generation(G1, G2):
    G1 = MatrixGroup(G1)
    G2 = MatrixGroup(G2)
    set2 = G2.gens()
    # 'in' does not work for infinite groups 
    return all((print(f"in group: {m}") or True) if m in G1 else (print(f"not in group: {m}") or False) for m in set2)



def bad_is_generated_by(target_gp, generating_gp, max_word_len=10):
    """
    For each matrix in target_gp's generators, attempt to express it as a word in 
    generating_gp's generators, avoiding immediate generator-inverse cancellations 
    and not attmepting to write generate the identity matrix 
    
    Returns:
        Dict: {target_matrix: word_as_indices or None}
    """
    print("Checking for generation...")
    generators = generating_gp.gens()
    targets = target_gp.gens()
    
    gen_list = list(generators) + [g.inverse() for g in generators]
    gen_indices = list(range(len(generators))) + list(range(-len(generators), 0))

    results = {}

    for target in targets:
        seen = {str(identity)}
        queue = deque([(identity, (), None)])  # track last index

        found = False
        for _ in range(max_word_len):
            next_queue = deque()
            while queue:
                mat, word, last_idx = queue.popleft()
                if target == identity:
                    results[target]= word
                    found = true
                    break 
                
                for g, idx in zip(gen_list, gen_indices):
                    # Prune immediate cancellations like (i, -i)
                    if last_idx is not None and idx == -last_idx:
                        continue

                    new_mat = mat * g
                    new_word = word + (idx,)
                    key = str(new_mat)
                    if key in seen:
                        continue

                    seen.add(key)
                    if new_mat == target:
                        results[target] = new_word
                        found = True
                        print(f"Found target:\n{target}\nas word: {new_word}")
                        break

                    next_queue.append((new_mat, new_word, idx))
                if found:
                    break
            if found:
                break
            queue = next_queue

        if not found:
            results[target] = None

    return results



def is_generated_by(target_gp, generating_gp, max_word_len=10):
    """
    For each matrix in `targets`, attempt to express it as a word in `generators`.
    
    Args:
        targets: list of matrices (e.g. from G2.gens())
        generators: list of matrices (e.g. from G1.gens())
        max_word_len: max length of words to search

    Returns:
        A dict: {target_matrix: word_as_indices or None}
    """
    print("Checking for generation...")
    generators = generating_gp.gens()
    targets = target_gp.gens()
    
    gen_list = list(generators) + [g.inverse() for g in generators]
    gen_indices = list(range(len(generators))) + list(range(-len(generators), 0))

    gen_dict = dict(zip(gen_indices, gen_list))

    results = {}

    for target in targets:
        seen = {str(identity)}
        queue = deque([(identity, ())])

        found = False
        for _ in range(max_word_len):
            next_queue = deque()
            if found:
                break
            while queue:
                mat, word = queue.popleft()

                if target == identity:
                    results[target]= ()
                    found = True
                    break 
                    
                for g, idx in zip(gen_list, gen_indices):
                    # Avoid immediately canceling moves like g * g⁻¹
                    # if word is not the empty word and the most recently multiplied index
                    # is that of the inverse of the current matrix about to be multiplied (g), 
                    # then don't multiply by g because it will just be g⁻¹* g
                    if idx >= 0:
                        if word and word[-1] == idx - len(generators):
                            continue
                    if idx < 0:
                        if word and word[-1] == idx + len(generators):
                            continue
                    new_mat = mat * g
                    new_word = word + (idx,)
                    key = str(new_mat)

                    if key in seen:
                        continue
                    seen.add(key)

                    if new_mat == target:
                        results[target] = new_word
                        found = True
                        break

                    next_queue.append((new_mat, new_word))
                if found:
                    break
            if found:
                break
            queue = next_queue

        if not found:
            results[target] = None

    return results
    
def extract_image_matrix_groups(homomorphisms):
    """
    Given a list of homomorphisms (as returned by all_homomorphisms_to_matrix_group),
    extract the matrix groups generated by the image of each homomorphism.

    Each homomorphism is a dict mapping generator strings to (matrix, word) tuples.
    This returns a list of MatrixGroup instances.
    """
    image_groups = []
    for hom in homomorphisms:
        image_matrices = [hom[g][0] for g in hom]  # extract only matrices
        G = MatrixGroup(image_matrices)
        image_groups.append(G)
    return image_groups

def word_to_matrix(word, generators):
    """
    Given a word as a tuple of generator indices (from is_generated_by),
    and a list of generators, compute the resulting matrix product.
    
    Positive index i means generators[i], negative -i means generators[-i].inverse().
    Note about python lists: generators[-i]=generators[i-1]
    
    Note: generators e.g.= [g0, g1, g2]
    gen_list = [g0, g1, g2, g0⁻¹, g1⁻¹, g2⁻¹]
    gen_indices = [0, 1, 2, -3, -2, -1]
    """
    gen_list = list(generators) + [g.inverse() for g in generators]
    gen_indices = list(range(len(generators))) + list(range(-len(generators), 0))

    gen_dict = dict(zip(gen_indices, gen_list))
    mat = identity
    for idx in word:
        mat *= gen_dict[idx]
    return mat

def trace_of_powers(mat, max_pow):
    mat = Matrix(mat)
    powers = mat.powers(max_pow)
    trace_of_powers = [m.trace() for m in powers]
    return trace_of_powers
    
  #from claude, giving it my is_gen_by function and asking for a 'meet in the middle' analogue 

def mitm_is_generated_by(target_gp, generating_gp, max_half_len=10):
    """
    Meet-in-the-middle analogue of is_generated_by.

    For each generator T of target_gp, searches for a word W in the generators
    of generating_gp such that W = T, where W has length at most 2 * max_half_len.

    Strategy:
        - Forward pass: enumerate all products of generators of length 0..max_half_len,
          storing them in a dict  forward = { matrix_key: (matrix, word) }.
        - For each target T, backward pass: enumerate all products P of generators
          of length 0..max_half_len, and check whether T * P^{-1} is in forward.
          If so, we have found  (T * P^{-1}) * P = T, i.e. word = left_word + right_word.

    Args:
        target_gp:     MatrixGroup whose generators we want to express as words.
        generating_gp: MatrixGroup whose generators we use to build words.
        max_half_len:  Maximum length of each half of the search (total depth = 2x this).

    Returns:
        dict: { target_matrix: word_as_tuple_of_indices, or None if not found }
              Word indices follow the same convention as is_generated_by:
              0  -> generator 0
              1  -> generator 1
              -2 -> inverse of generator 0
              -1 -> inverse of generator 1
    """
    generators = [Matrix(QQ, g) for g in generating_gp.gens()]
    targets     = [Matrix(QQ, g) for g in target_gp.gens()]

    n = len(generators)
    gen_list    = generators + [g.inverse() for g in generators]
    gen_indices = list(range(n)) + list(range(-n, 0))

    def mat_key(m):
        # Hashable, exact key for a QQ matrix
        return tuple(tuple(QQ(e) for e in row) for row in m.rows())

    def no_cancel(word, idx):
        """Return True if appending idx to word does NOT immediately cancel."""
        if not word:
            return True
        last = word[-1]
        # idx and last cancel iff one is the inverse index of the other
        if idx >= 0 and last == idx - n:
            return False
        if idx < 0  and last == idx + n:
            return False
        return True

    # ?? Forward pass: build table of all words of length 0..max_half_len ??????
    print(f"Building forward table (half-depth {max_half_len})...")
    forward = {}                            # mat_key -> (matrix, word_tuple)
    forward[mat_key(identity)] = (identity, ())
    current_level = [(identity, ())]

    for depth in range(max_half_len):
        next_level = []
        for mat, word in current_level:
            for g, idx in zip(gen_list, gen_indices):
                if not no_cancel(word, idx):
                    continue
                new_mat  = mat * g
                new_word = word + (idx,)
                key = mat_key(new_mat)
                if key not in forward:
                    forward[key] = (new_mat, new_word)
                    next_level.append((new_mat, new_word))
        current_level = next_level
        print(f"  depth {depth+1}: {len(forward)} entries so far")

    print(f"Forward table complete: {len(forward)} distinct matrices.\n")

    # ?? For each target, do the backward pass ?????????????????????????????????
    results = {}

    for target in targets:
        print(f"Searching for target:\n{target}")

        # Check length-0 case
        key0 = mat_key(target)
        if key0 in forward:
            left_word = forward[key0][1]
            results[target] = left_word
            print(f"  Found at depth {len(left_word)}: word = {left_word}\n")
            continue

        found      = False
        back_level = [(identity, ())]
        back_seen  = {mat_key(identity)}

        for depth in range(max_half_len):
            next_back = []
            for right_mat, right_word in back_level:
                for g, idx in zip(gen_list, gen_indices):
                    if not no_cancel(right_word, idx):
                        continue
                    new_right_mat  = right_mat * g
                    new_right_word = right_word + (idx,)
                    rkey = mat_key(new_right_mat)
                    if rkey in back_seen:
                        continue
                    back_seen.add(rkey)

                    # We want: left_mat * new_right_mat = target
                    # => left_mat = target * new_right_mat^{-1}
                    candidate_left = target * new_right_mat.inverse()
                    lkey = mat_key(candidate_left)

                    if lkey in forward:
                        left_word  = forward[lkey][1]
                        full_word = left_word + new_right_word
                        results[target] = full_word
                        found = True
                        total_len = len(left_word) + len(new_right_word)
                        print(f"  Found at total depth {total_len}: word = {full_word}")
                        break

                    next_back.append((new_right_mat, new_right_word))

                if found:
                    break
            if found:
                break
            back_level = next_back

        if not found:
            results[target] = None
            print(f"  Not found within depth {2 * max_half_len}.")
        print()

    return results