# Modified for SpectroIBIS by Brodie Bulcock 31/10/24
# See https://github.com/berhane/arbalign/blob/master/ArbAlign-scipy.py for pre-modded source code


from numpy import array, dot, linalg, sqrt, transpose
from scipy.optimize import linear_sum_assignment
from collections import Counter
from PrinCoords_modded import main as prin_coords


def kabsch(A, B):
    """
    Kabsch Algorithm as implemented by Jimmy Charnley Kromann

    Calculate RMSD between two XYZ files

    by: Jimmy Charnley Kromann <jimmy@charnley.dk> and
    Lars Andersen Bratholm <larsbratholm@gmail.com>
    project: https://github.com/charnley/rmsd
    license: https://github.com/charnley/rmsd/blob/master/LICENSE

    A - set of coordinates
    B - set of coordinates

    Performs the kabsch algorithm to calculate the RMSD between A and B

    Returns an RMSD
    """
    A_new = array(A)
    A_new = A_new - sum(A_new) / len(A_new)
    A = A_new
    B_new = array(B)
    B_new = B_new - sum(B_new) / len(B_new)
    B = B_new

    # Compute covariance matrix

    C = dot(transpose(A), B)

    # Compute singular value decomposition (SVD)

    V, S, W = linalg.svd(C)
    d = (linalg.det(V) * linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    # Compute rotation matrix

    U = dot(V, W)

    # Rotate A

    A = dot(A, U)

    return rmsd(A, B)


def rmsd(V, W):
    """
    V - set of coordinates
    W - set of coordinates

    Returns root-mean-square deviation from two sets of vectors V and W.
    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i] - w[i]) ** 2.0 for i in range(D)])
    return sqrt(rmsd / N)


def coords_sorter(coords):
    """Returns a tuple (a, b) sorted by atom labels and coordinates
     where a is a list of coordinate labels and b is a set of principal coordinates
    (i.e) a = ["O", "H", "H"], b = [[x0,y0,z0],[x1,y1,z1],[x2,y2,z2]]
    Sorts the file by atom labels first and coordinates second
    such that atoms of the same label/type are grouped together. Returns aligned coords too
    """

    sortedlabels = []
    sortedcoords = []
    sortedorder = []
    sortedlines = []
    atomcount = 0

    # First convert coords to principal coordinates

    aligned_coords = prin_coords(coords)

    # Rearrange coords to the format accepted by ArbAlign [[element, x, y, z, atomcount], ...]

    for i in range(len(coords[0])):
        sortedlines.append([coords[0][i], aligned_coords[i][0], aligned_coords[i][1], aligned_coords[i][2], atomcount])
        atomcount += 1
    # sort by element followed by first coordinate (x) then second coordinate (y)

    sortedlines.sort(key=lambda x: (x[0], x[1], x[2]))
    for i in range(atomcount):
        sortedlabels.append(sortedlines[i][0])
        sortedcoords.append([float(sortedlines[i][1]), float(sortedlines[i][2]), float(sortedlines[i][3])])
        sortedorder.append(sortedlines[i][4])
    NA = len(sortedlabels)
    return sortedlabels, sortedcoords, NA, sortedorder, aligned_coords


def parse_for_atom(labels, coords, atom):
    """
    labels - a list of coordinate labels
    coords - a set of coordinates
    atom - the atom that is to be parsed

    Returns a set of coordinates corresponding to parsed atom
    """
    atom_coords = []
    for i in range(len(labels)):
        if labels[i] == atom:
            atom_coords.append(coords[i])
    return atom_coords


def transform_coords(coords, swap, reflect):
    """
    coords - a set of coordinates
    swap - the swap transformation (i.e. (0, 2, 1) --> (x, z, y))
    reflect - the reflection transformation (i.e. (1, -1, 1) --> (x, -y, z))

    Returns the transformed coordinates
    """
    new_coords = []
    for i in range(len(coords)):
        new_coords.append(
            [coords[i][swap[0]] * reflect[0], coords[i][swap[1]] * reflect[1], coords[i][swap[2]] * reflect[2]]
        )
    return new_coords


def permute_atoms(coords, permutation, atom_indices):
    """
    coords - a set of coordinates
    permuation - a permutation of atoms
    atom_indices - indices of all desired atoms in [coords]

    Returns the coordinates after permuting just the specified atom
    """
    new_coords = coords[:]
    for i in range(len(permutation)):
        j = atom_indices[permutation[i]]
        k = atom_indices[i]
        new_coords[k] = coords[j]
    return new_coords


def permute_all_atoms(labels, coords, permutation):
    """
    labels - atom labels
    coords - a set of coordinates
    permuation - a permutation of atoms

    Returns the permuted labels and coordinates
    """
    new_coords = coords[:]
    new_labels = labels[:]
    for i in range(len(permutation)):
        new_coords[permutation[i]] = coords[i]
        new_labels[permutation[i]] = labels[i]
    return new_labels, new_coords


def get_atom_indices(labels, atom):
    """
    labels - a list of coordinate labels ("Elements")
    atom - the atom whose indices in labels are sought
    Returns a list of all location of [atom] in [labels]
    """
    indices = []
    for i in range(len(labels)):
        if labels[i] == atom:
            indices.append(i)
    return indices


def main(coords_a, coords_b):

    description = """
This code uses the Kuhn-Munkres or Hungarian algorithm to optimally align two
arbitrarily ordered isomers. Given two isomers A and B whose Cartesian
coordinates are given in XYZ format, it will optimally align B on A to minimize
the Kabsch root-mean-square deviation (RMSD) between structure A and B after

1) a Kuhn-Munkres assignment/reordering (quick) 
2) a Kuhn-Munkres assignment/reordering factoring in axes swaps and reflections (~48x slower)

We recommend the second method although the first one would still be better
than RMSD calculations without atom reorderings.  

A web server with this implementation is available at http://www.arbalign.org

While this script is kept as minimal as possible in order to ensure ease of use
and portability, it does require these two Python packages beyond what's
included in standard python installations.

1) Python Numpy module 
2) Python Hungarian module by Harold Cooper 
   (Hungarian: Munkres' Algorithm for the Linear Assignment Problem in Python.
   https://github.com/Hrldcpr/Hungarian) 
   This is a wrapper to a fast C++ implementation of the Kuhn-Munkres algorithm.
   The installation instructions are described at https://github.com/Hrldcpr/Hungarian

Other optional tools are:

1) PrinCoords.py - using principal coordinates generally yields better
   alignment (lower RMSDs).  A Python script to convert molecules from arbitrary
   to principal coordinate system is included.

2) In cases where one wants to use atom types including connectivity and
   hybridization information, it is necessary to use OpenBabel to convert the
   Cartesian coordinates to SYBYL Mol2 (sy2) and MNA (mna) formats.  

The best way to take advantage of these two optional tools is probably to use
the attached driver script (ArbAlign-driver.py) The syntax looks like

   Usage: ArbAlign-driver.py -<flag> <filename_1.xyz> <filename_2.xyz>"
        : where the <flag> is "
        : -l   match by atom or element label "
        : -t   match by SYBYL atom type"
        : -c   match by NMA atom connectivity type"
        "
     Eg.: ArbAlign-driver.py -b -N cluster1.xyz cluster2.xyz"
        : ArbAlign-driver.py -T cluster1.xyz cluster2.xyz"
        : ArbAlign-driver.py -C cluster1.xyz cluster2.xyz"
        "
   This matches the Cartesian coordinates of the file1 and file2 using the \
         Kuhn-Munkres algorithm based on atom labels (-l), type (-t) or \
         connectivity (-t). "
   It produces s-file1.xyz and s-file2-matched.xyz which are the sorted and \
         matched file1 and file2.xyz, respectively."

"""

    epilog = """
The code will provide the following:
1) The initial Kabsch RMSD
2) The final Kabsch RMSD after the application of the Kuhn-Munkres algorithm 
3) The coordinates corresponding to the best alignment of B on A to a file called B-aligned_to-A.xyz

If you find this script useful for any publishable work, please cite the corresponding paper:
  Berhane Temelso, Joel M. Mabey, Toshiro Kubota, Nana Appiah-padi, George C. Shields
  J. Chem. Info. Model. 2017, 57(5), 1045-1054 
"""

    """
   Read in the original coordinates and labels of xyz1 and xyz2, 
   and sort them by atom labels so that atoms of the same label/name are grouped together

   Then, count how many types of atoms, and determine their numerical frequency
   """
    a_labels, a_coords, NA_a, order, a_prin_coords = coords_sorter(coords_a)
    Uniq_a = list(set(a_labels))
    list.sort(Uniq_a)
    N_uniq_a = len(Uniq_a)
    Atom_freq_a = dict(Counter(a_labels))

    b_labels, b_coords, NA_b, junk, b_prin_coords = coords_sorter(coords_b)
    Uniq_b = list(set(b_labels))
    list.sort(Uniq_b)
    N_uniq_b = len(Uniq_b)
    Atom_freq_b = dict(Counter(b_labels))

    """
   If the number and type of atoms in the two structures are not equal, exit with 
   an error message
   """
    if (NA_a == NA_b) & (Uniq_a == Uniq_b) & (Atom_freq_a == Atom_freq_b):
        num_atoms = NA_a
        num_uniq = N_uniq_a
        Uniq = Uniq_a  # list(set(a_labels))
        Atom_freq = Atom_freq_a
    A_all = array(a_coords)
    A_all = A_all - sum(A_all) / len(A_all)
    B_all = array(b_coords)
    B_all = B_all - sum(B_all) / len(B_all)

    """
   Dynamically generate hashes of coordinates and atom indices for every atom type
   """
    a_Coords = {}
    a_Indices = {}
    b_Coords = {}
    b_Indices = {}
    Perm = {}
    for i in range(len(Uniq)):
        a_Coords[Uniq[i]] = "a_" + Uniq[i] + "coords"
        b_Coords[Uniq[i]] = "b_" + Uniq[i] + "coords"
        a_Indices[Uniq[i]] = "a_" + Uniq[i] + "indices"
        b_Indices[Uniq[i]] = "b_" + Uniq[i] + "indices"
        Perm[Uniq[i]] = "perm_" + Uniq[i]
        vars()[Perm[Uniq[i]]] = []
        vars()[a_Coords[Uniq[i]]] = parse_for_atom(a_labels, a_coords, str(Uniq[i]))
        vars()[a_Indices[Uniq[i]]] = get_atom_indices(a_labels, str(Uniq[i]))
        vars()[b_Coords[Uniq[i]]] = parse_for_atom(b_labels, b_coords, str(Uniq[i]))
        vars()[b_Indices[Uniq[i]]] = get_atom_indices(b_labels, str(Uniq[i]))
    l = 0
    A = array(vars()[a_Coords[Uniq[l]]])
    A = A - sum(A) / len(A)
    B = array(vars()[b_Coords[Uniq[l]]])
    B = B - sum(B) / len(B)

    """
   For each atom type, we can do a Kuhn-Munkres assignment in the initial 
   coordinates or the many swaps and reflections thereof

   If a single Kuhn-Munkres assignment is requested with a -s or --simple flag,
   no swaps and reflections are considered. Otherwise, the default is to perform 
   a combination of 6 axes swaps and 8 reflections and do Kuhn-Munkres assignment 
   on all 48 combinations. 
   """
    swaps = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
    reflects = [(1, 1, 1), (-1, 1, 1), (1, -1, 1), (1, 1, -1), (-1, -1, 1), (-1, 1, -1), (1, -1, -1), (-1, -1, -1)]
    B_t = []
    for i in swaps:
        for j in reflects:
            B_t.append([transform_coords(B, i, j), i, j])
    rmsds = []
    # Performs the munkres algorithm on each set of transformed coordinates

    for i in range(len(B_t)):
        l = 0
        cost_matrix = array([[linalg.norm(a - b) for b in B_t[i][0]] for a in A])
        LAP0, LAP1 = linear_sum_assignment(cost_matrix)
        vars()[Perm[Uniq[l]]] = []
        for j in range(len(LAP1)):
            vars()[Perm[Uniq[l]]] += [(j, LAP1[j])]
        vars()[Perm[Uniq[l]]] = sorted(vars()[Perm[Uniq[l]]], key=lambda x: x[0])
        vars()[Perm[Uniq[l]]] = [x[1] for x in vars()[Perm[Uniq[l]]]]

        # If there's more than one atom type, loop through each unique atom type

        if num_uniq == 1:
            b_perm = permute_atoms(b_coords, vars()[Perm[Uniq[l]]], vars()[b_Indices[Uniq[l]]])
            b_final = transform_coords(b_perm, B_t[i][1], B_t[i][2])
            rmsds.append([kabsch(a_coords, b_final), B_t[i][1], B_t[i][2], b_final, vars()[Perm[Uniq[l]]]])
            rmsds = sorted(rmsds, key=lambda x: x[0])
        else:
            b_perm = permute_atoms(b_coords, vars()[Perm[Uniq[l]]], vars()[b_Indices[Uniq[l]]])
            b_trans = transform_coords(b_perm, B_t[i][1], B_t[i][2])
            while l < num_uniq:
                if l > 0:
                    vars()[b_Coords[Uniq[l]]] = parse_for_atom(b_labels, b_final, Uniq[l])
                else:
                    vars()[b_Coords[Uniq[l]]] = parse_for_atom(b_labels, b_trans, Uniq[l])
                x = array(vars()[b_Coords[Uniq[l]]])
                y = array(vars()[a_Coords[Uniq[l]]])
                cost_matrix = array([[linalg.norm(a - b) for b in x] for a in y])
                LAP0, LAP1 = linear_sum_assignment(cost_matrix)
                vars()[Perm[Uniq[l]]] = []
                for k in range(len(LAP1)):
                    vars()[Perm[Uniq[l]]] += [(k, LAP1[k])]
                vars()[Perm[Uniq[l]]] = sorted(vars()[Perm[Uniq[l]]], key=lambda x: x[0])
                vars()[Perm[Uniq[l]]] = [x[1] for x in vars()[Perm[Uniq[l]]]]
                b_final = permute_atoms(b_trans, vars()[Perm[Uniq[l]]], vars()[b_Indices[Uniq[l]]])
                b_trans = b_final
                l += 1
                q = l - 1
                rmsds.append([kabsch(a_coords, b_final), B_t[i][1], B_t[i][2], b_final])
                rmsds = sorted(rmsds, key=lambda x: x[0])
    FinalRMSD = float(rmsds[0][0])

    b_final_labels, b_final_coords = permute_all_atoms(b_labels, rmsds[0][3], order)
    return FinalRMSD, a_prin_coords, b_final_coords
