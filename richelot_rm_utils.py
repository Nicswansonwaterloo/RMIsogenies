from sage.all import Matrix, Integers

from sage.schemes.hyperelliptic_curves.invariants import absolute_igusa_invariants_kohel
from dependencies.Castryck_Decru_SageMath.richelot_aux import FromJacToJac
from dependencies.Theta_SageMath.theta_structures.couple_point import CouplePoint


def golden_ratio_endomorphism(cp_pt):
    """ The endomorphism of E x E given by (P, Q) |-> (Q, P + Q). """
    P, Q = cp_pt
    return CouplePoint(Q, P + Q)

def golden_ratio_action_on_symplectic_torsion(ell=2, e=1):
    """ Returns the matrix of the action of the golden ratio endomorphism on the ell^e-torsion,
        with respect to the basis [ (P2_1, 0), (0, P2_2), (Q2_1, 0), (0, Q2_2) ].
    """
    Zle = Integers(ell**e)
    return Matrix(Zle, [[0, 1, 0, 0],
                        [1, 1, 0, 0],
                        [0, 0, 0, 1],
                        [0, 0, 1, 1]])

def get_isogeny_from_jacobian_two_kernel(kernel_generators):
    # TODO: Detect when codomain is jacobian vs product
    return FromJacToJac(kernel_generators[0], kernel_generators[1])

class RMVertex:
    """
        A vertex in the (2,2)-isogeny graph of supersingular abelian surfaces with RM.
        
        - The vertex can either be a product of two elliptic curves or the Jacobian of a genus 2 curve. This is stored in the attribute 'variety'.
        - The RM is represented by its action on 2^r-torsion, for some r. This integer is stored in the attribute 'r'.
        - The action on the 2^r-torsion is stored in the attribute 'action', which is a 4x4 matrix with coefficients in Z/(2^r)Z.
    """
    def __init__(self, variety, r, action, torsion_basis):
        self.variety = variety
        self.r = r
        self.action = action
        self.torsion_basis = torsion_basis

    def is_product(self):
        return isinstance(self.variety, tuple)

    def is_jacobian(self):
        return not self.is_product()

    def __repr__(self):
        if isinstance(self.variety, tuple):
            return f"Product: ({self.variety[0].j_invariant()}, {self.variety[1].j_invariant()})"
        return f"Jacobian: {absolute_igusa_invariants_kohel(self.variety)}"
    
    def __eq__(self, other):
        return repr(self.variety) == repr(other.variety)