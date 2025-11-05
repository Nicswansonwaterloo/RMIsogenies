from sage.all import Matrix, GF, vector, HyperellipticCurve, EllipticCurve
from richelot_rm.genus_two_structures import GenusTwoJacobianStructure, GenusTwoProductStructure
from richelot_rm.jacobian_point import JacobianPoint
from richelot_rm.product_point import ProductPoint


def is_2_kernel_jac(kernel):
    if not (isinstance(kernel, (tuple, list)) and len(kernel) == 2):
        raise ValueError("Kernel must be a tuple of length 2.")
    if not all(isinstance(P, JacobianPoint) for P in kernel):
        raise TypeError("Kernel generators must be Jacobian points.")
    
    gen1, gen2 = kernel
    h = gen1.parent().h
    G1, _ = gen1
    G2, _ = gen2
    G3, r3 = h.quo_rem(G1 * G2)
    return r3 == 0


def is_2_kernel_jac_split(kernel):
    if not is_2_kernel_jac(kernel):
        raise ValueError("The given kernel does not define a 2-isogeny on a jacobian of genus 2")
    
    gen1, gen2 = kernel
    g2_structure = gen1.parent()
    h = g2_structure.h

    G1, _ = gen1
    G2, _ = gen2
    G3 = h // (G1 * G2)

    delta = Matrix(G.padded_list(3) for G in (G1, G2, G3))
    if delta.determinant():
        # Determinant is non-zero, no splitting
        return False

    return True


def jacobian_to_product_2_isogeny(kernel):
    if not is_2_kernel_jac_split(kernel):
        raise ValueError("The given kernel does not define a split 2-isogeny on a jacobian of genus 2")
    
    gen1, gen2 = kernel
    g2_structure = gen1.parent()
    h = g2_structure.h
    x = g2_structure.x

    G1, _ = gen1
    G2, _ = gen2
    G3 = h // (G1 * G2)
    M = Matrix(G.padded_list(3) for G in (G1,G2,G3))

    homography_needed = not M.column(1).is_zero()
    if not homography_needed:
        # No homography needed
        H11, H21, H31 = M.column(2)
        H10, H20, H30 = M.column(0)
    else:
        # Find homography
        u, v, w = M.right_kernel().gen()
        d = u/2
        (ad, _), (b, _) = (x**2 - v*x + w*d/2).roots()
        a = ad/d

        # Apply transform G(x) -> G((a*x+b)/(x+d))*(x+d)^2
        # The coefficients of x^2 are M * (1, a, a^2)
        # The coefficients of 1 are M * (d^2, b*d, b^2)
        H11, H21, H31 = M * vector([1, a, a*a])
        H10, H20, H30 = M * vector([d*d, b*d, b*b])
        assert G1((a*x+b)/(x+d))*(x+d)**2 == H11*x**2+H10

    h2 = (H11*x**2+H10)*(H21*x**2+H20)*(H31*x**2+H30)
    H2 = HyperellipticCurve(h2)

    p1 = (H11*x+H10)*(H21*x+H20)*(H31*x+H30)
    p2 = (H11+H10*x)*(H21+H20*x)*(H31+H30*x)
    # We will need to map to actual elliptic curve
    p1norm = (x + H10*H21*H31)*(x + H20*H11*H31)*(x + H30*H11*H21)
    p2norm = (x + H11*H20*H30)*(x + H21*H10*H30)*(x + H31*H10*H20)
    E1 = EllipticCurve([0, p1norm[2], 0, p1norm[1], p1norm[0]])
    E2 = EllipticCurve([0, p2norm[2], 0, p2norm[1], p2norm[0]])
    codomain = GenusTwoProductStructure(E1, E2)

    def morphE1(x, y):
        # from y^2=p1 to y^2=p1norm
        return (H11*H21*H31*x, H11*H21*H31*y)
    def morphE2(x, y):
        # from y^2=p1 to y^2=p2norm
        return (H10*H20*H30*x, H10*H20*H30*y)
    
    def isogeny(D: JacobianPoint):
        U, V = D # Lets call the roots of U: x_a, x_b
        if homography_needed:
            # apply homography
            # y = v1 x + v0 =>
            U_ = U[0] * (x+d)**2 + U[1]*(a*x+b)*(x+d) + U[2]*(a*x+b)**2
            V_ = V[0] * (x+d)**3 + V[1]*(a*x+b)*(x+d)**2
            V_ = V_ % U_
        else:
            U_, V_ = U, V

        # Note that y_a = v1 x_a + v0, y_b = v1 x_b + v0 are the y-coordinates corresponding to x_a, x_b
        v1, v0 = V_[1], V_[0]
        # prepare symmetric functions
        s = - U_[1] / U_[2] # SUM of roots of U_: x_a + x_b
        p = U_[0] / U_[2] # PRODUCT of roots of U_: x_a * x_b
        assert p != 0 # nether x_a not x_b can be zero. Perhaps a change of coordinates is needed.

        # For E1: x -> x^2 := z, and y -> y := w
        # We need to compute the divisor on E1, corresponding to the sum of the image of (x_a, y_a) and (x_b, y_b), then sum it to get the image on E1.
        if s.is_zero():
            # This is the case where the divisor on E1 is the sum of two points with the same x-coordinate: x_a^2 = (-x_b)^2
            if v0.is_zero():
                # Then v0 = 0 <=> y_a = - y_b
                E1_image = E1(0)
            else:
                # Here y_a = y_b. This gives v1 x_a + v0 = v1 (-x_a) + v0 => 2 v1 x_a = 0 => v_1 = 0
                # So we double the point (x_a^2, y_a) = (-p, v_0) on E1.
                E1_image = 2 * E1(morphE1(-p, v0)) #How do figure out the sign of y_a? 
        else:
            U1 = x**2 - (s*s - 2*p)*x + p**2 # roots are x_a^2, x_b^2
            if not v0.is_zero():
                # The general case
                V1 = (p1 - v1**2 * x + v0**2) / (2*v0) # V1(x_a^2) = w_a, V1(x_b^2) = w_b where w_a^2 = p1(x_a^2)
                V1 = V1 % U1 # Reduce to Mumford coordinates
                U1red = (p1 - V1**2) // U1 # V1^2 - p1 have the roots x_a^2, x_b^2 BUT ALSO the extra root x_c. which simulateniously lives on E1, but also lies on the line through (x_a, y_a) and (x_b, y_b). Hence x_c is the x-coordinate of P_{x_a^2} + P_{x_b^2} on E1.
                xP1 = -U1red[0] / U1red[1] # Here U1red is linear, and we are solving for the extra root x_c
                yP1 = V1(xP1)
            else:
                # Special case when v0 = 0
                # Then p1 - v1^2 x = 0 at x_a^2, x_b^2 and the third root is immediately calculable
                U1red = (p1 - v1**2 * x) // U1
                xP1 = -U1red[0] / U1red[1]
                # Now its a matter of finding V1(z) the line through (x_a^2, w_a) and (x_b^2, w_b) to recover the y-coordinate.
                # The slope is (w_b - w_a) / (x_b^2 - x_a^2) = (v1 x_b - v1 x_a) / (x_b^2 - x_a^2) = v1 / (x_a + x_b) = v1 / s
                yP1 = (v1 / s) * (xP1 + p)
            
            assert yP1**2 == p1(xP1)
            E1_image = E1(morphE1(xP1, yP1))

        # For E2: x -> 1/x^2 := z AND y -> y/x^3 := w
        if s.is_zero():
            # Points will still be mapped to the same z-coordinate
            # However, w_a = (v1 x_a + v0)/x_a^3 and w_b = (v1 x_a - v0)/(x_a)^3 in general will have w_a != - w_b
            # w_a = -w_b  <=> v1 = 0
            if v1.is_zero():
                E2_image = E2(0)
            # w_a = w_b <=> v0 = 0
            elif v0.is_zero():
                # In this case w_a = w_b = v1 x_a / x_a^3 = v1 / x_a^2
                E2_image = 2 * E2(morphE2(1/-p, - v1/p)) # This may be incorrect?
            else:
                raise ValueError("Cannot have s=0 with v1, v0 both non-zero. This is mathematically impossible unless x_a = 0...") # this is because y_a^2 = y_b^2 implies (v1 x_a + v0)^2 = (v1 (-x_a) + v0)^2 => 4 v1 v0 x_a = 0
        else:
            C = v1*(s*s-2*p) + v0*s # This is (x_ay_a + x_by_b)
            U2 = x**2 - (s*s-2*p)/p**2*x + 1/p**2 # The roots are 1/x_a^2 := z_a, 1/x_b^2 := z_b
            if C == 0:
                # We have sv0 = v1 (s^2 - 2p) with s != 0. This means x_a y_a = - x_b y_b
                # This means we can directly compute the slope between (z_a, w_a) and (z_b, w_b) on E2
                slope = v1 + (v0*(s**2 - 2*p) / (p*s))
                intercept = v0 / p
                V2 = slope * x + intercept
                V2 = V2 % U2
                U2red = (p2 - V2**2) // U2
                xP2 = -U2red[0] / U2red[1]
                yP2 = V2(xP2)
                assert yP2**2 == p2(xP2)
                E2_image = E2(morphE2(xP2, yP2))
            else:
                # General case
                V21 = x**2 * C # this is z^2 * (x_ay_a + x_by_b)
                V20 = p2 + x**4 * (p*(v1**2*p + v1*v0*s + v0**2)) # w^2 + z^4(x_a x_b y_ay_b)

                # V21 * w = V20 modulo U2 (for (z_a, w_a) and (z_b, w_b))
                _d, V21inv, _ = V21.xgcd(U2)
                assert _d.is_one(), f"GCD not 1: {_d}\n s: {s}, p: {p}\n V21: {V21}, U2: {U2}"
                V2 = (V21inv * V20) % U2
                assert V2**2 % U2 == p2 % U2

                # Reduce coordinates
                U2red = (p2 - V2**2) // U2
                xP2 = -U2red[0] / U2red[1]
                yP2 = V2(xP2)
                assert yP2**2 == p2(xP2)
                E2_image = E2(morphE2(xP2, yP2))

        return ProductPoint(E1_image, E2_image)
    
    return codomain, isogeny
    

def jacobian_to_jacobian_2_isogeny(kernel):
    gen1, gen2 = kernel
    g2_structure = gen1.parent()
    J = g2_structure.jac
    h = g2_structure.h
    x = g2_structure.x

    G1, _ = gen1
    G2, _ = gen2
    if (G1[2] != 1 and G1[2] != 0):
        G1 = G1 / G1[2]
    if (G2[2] != 1 and G2[2] != 0):
        G2 = G2 / G2[2]
    
    G3 = h // (G1 * G2)
    M = Matrix(G.padded_list(3) for G in (G1,G2,G3))
    delta = M.inverse()

    H1 = -delta[0][0]*x**2 + 2*delta[1][0]*x - delta[2][0]
    H2 = -delta[0][1]*x**2 + 2*delta[1][1]*x - delta[2][1]
    H3 = -delta[0][2]*x**2 + 2*delta[1][2]*x - delta[2][2]

    # This is the so-called Richelot Correspondance from Ben Smith's thesis:
    h_codomain = H1 * H2 * H3
    codomain = GenusTwoJacobianStructure(h_codomain)
    
    def isogeny(D: JacobianPoint):
        U, V = D
        
        # Make monic
        if not U[2].is_one():
            U = U / U[2]
        V = V % U
        if U == G1 or U == G2:
            return JacobianPoint(J(0))

        # Sum and product of (xa, xb)
        s, p = -U[1], U[0]
        # Compute X coordinates (non reduced, degree 4)
        g1red = G1 - U
        g2red = G2 - U
        # Make sure everything is monic!
        assert g1red[2].is_zero() and g2red[2].is_zero()
        g11, g10 = g1red[1], g1red[0]
        g21, g20 = g2red[1], g2red[0]
        # see above
        Px = (g11*g11*p + g11*g10*s + g10*g10) * (H1 * H1) \
           + (2*g11*g21*p + (g11*g20+g21*g10)*s + 2*g10*g20) * (H1 * H2) \
           + (g21*g21*p + g21*g20*s + g20*g20) * (H2 * H2)

        # Compute Y coordinates (non reduced, degree 3)
        assert V[2].is_zero()
        v1, v0 = V[1], V[0]
        # coefficient of y^2 is V(xa)V(xb)
        Py2 = v1*v1*p + v1*v0*s + v0*v0
        # coefficient of y is h1(x) * (V(xa) Gred1(xb) (x-xb) + V(xb) Gred1(xa) (x-xa))
        # so we need to symmetrize:
        # V(xa) Gred1(xb) (x-xb)
        # = (v1*xa+v0)*(g11*xb+g10)*(x-xb)
        # = (v1*g11*p + v1*g10*xa + v0*g11*xb + v0*g10)*x
        # - xb*(v1*g11*p + v1*g10*xa + v0*g11*xb + v0*g10)
        # Symmetrizing xb^2 gives u1^2-2*u0
        Py1 = (2*v1*g11*p + v1*g10*s + v0*g11*s + 2*v0*g10)*x \
          - (v1*g11*s*p + 2*v1*g10*p + v0*g11*(s*s-2*p) + v0*g10*s)
        Py1 *= H1
        # coefficient of 1 is Gred1(xa) Gred1(xb) h1(x)^2 U(x)
        Py0 = H1 * H1 * U * (g11*g11*p + g11*g10*s + g10*g10)

        # Now reduce the divisor, and compute Cantor reduction.
        # Py2 * y^2 + Py1 * y + Py0 = 0
        # y = - (Py2 * hnew + Py0) / Py1
        _, Py1inv, _ = Py1.xgcd(Px)
        Py = (- Py1inv * (Py2 * h_codomain + Py0)) % Px
        assert Px.degree() == 4
        assert Py.degree() <= 3

        Dx = ((h_codomain - Py ** 2) // Px)
        Dy = (-Py) % Dx
        assert (h_codomain - Dy**2) % Dx == 0, f"not div {(h_codomain - Dy**2) % Dx}"
        jac_divisor = codomain.jac([Dx, Dy])
        return JacobianPoint(jac_divisor)
    
    return codomain, isogeny

def get_symplectic_two_torsion_jac(jac_structure: GenusTwoJacobianStructure):
    J = jac_structure.jac
    Rx = jac_structure.Rx
    x = jac_structure.x
    h = jac_structure.h
    roots = h.roots(multiplicities=False)
    if h.degree() == 6:
        assert len(roots) == 6 
    elif h.degree() == 5:
        assert len(roots) == 5 
        
    T1x = (x - roots[0]) * (x - roots[1])
    T3x = (x - roots[0]) * (x - roots[2])  # must share a factor with T1, none others
    T2x = (x - roots[3]) * (x - roots[4])  # must share a factor with T4, none others
    if h.degree() == 5:
        T4x = (x - roots[3])  # must share a factor with T2, none others
    else:
        T4x = (x - roots[3]) * (x - roots[5])  # must share a factor with T2, none others

    T1 = JacobianPoint(J([Rx(T1x), Rx(0)]))
    T2 = JacobianPoint(J([Rx(T2x), Rx(0)]))
    T3 = JacobianPoint(J([Rx(T3x), Rx(0)]))
    T4 = JacobianPoint(J([Rx(T4x), Rx(0)]))

    return [T1, T2, T3, T4]


def compute_2_isogeny_from_jacobian(kernel):
    if is_2_kernel_jac_split(kernel):
        return jacobian_to_product_2_isogeny(kernel)
    return jacobian_to_jacobian_2_isogeny(kernel)
