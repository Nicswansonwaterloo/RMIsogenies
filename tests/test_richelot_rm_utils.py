from sage.all import vector
from richelot_rm_utils import *
from richelot_products import get_arbitrary_square_example

def test_golden_ratio_endomorphism():
    p = 2**11 * 3 - 1
    E1, E2 = get_arbitrary_square_example(p)
    P2_1, Q2_1 = E1.torsion_basis(2)
    P2_2, Q2_2 = E2(P2_1), E2(Q2_1)
    def vector_to_point(v, basis):
        return v[0] * basis[0] + v[1] * basis[1] + v[2] * basis[2] + v[3] * basis[3]

    basis = [CouplePoint(P2_1, E2(0)), CouplePoint(E1(0), P2_2), CouplePoint(Q2_1, E2(0)), CouplePoint(E1(0), Q2_2)]
    M_phi = golden_ratio_action_on_symplectic_torsion(2, 11)
    Zle = Integers(2**11)
    for i in range(4):
        v = vector(Zle, [0, 0, 0, 0])
        v[i] = 1
        assert vector_to_point(M_phi * v, basis) == golden_ratio_endomorphism(vector_to_point(v, basis))


    P2e_1, Q2e_1 = E1.torsion_basis(2**11)
    P2e_2, Q2e_2 = E2(P2e_1), E2(Q2e_1)
    basis = [CouplePoint(P2e_1, E2(0)), CouplePoint(E1(0), P2e_2), CouplePoint(Q2e_1, E2(0)), CouplePoint(E1(0), Q2e_2)]
    M_phi = golden_ratio_action_on_symplectic_torsion(2, 11)
    Zle = Integers(2**11)
    for i in range(4):
        v = vector(Zle, [0, 0, 0, 0])
        v[i] = 1
        assert vector_to_point(M_phi * v, basis) == golden_ratio_endomorphism(vector_to_point(v, basis))

test_golden_ratio_endomorphism()