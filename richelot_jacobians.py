from sage.all import Matrix, GF, PolynomialRing, EllipticCurve, HyperellipticCurve, vector, ZZ

from couple_point import CouplePoint

def is_jac_kernel_split(h, kernel_generators):
    """
    Checks if the kernel of a Jacobian isogeny splits.

    :param h: The polynomial defining the hyperelliptic curve.
    :param kernel_generators: A list of two kernel generators.
    :return: A tuple (is_split, G1, G2, G3) where is_split is a boolean.
    """
    D11, D12 = kernel_generators[0]
    D21, D22 = kernel_generators[1]
    G1 = D11
    G2 = D21
    G3, r3 = h.quo_rem(G1 * G2)
    assert r3 == 0, f"h roots: {h.roots()}, \n G1 roots: {G1.roots()}, \n G2 roots: {G2.roots()}\n G3 roots: {G3.roots()}"

    delta = Matrix(G.padded_list(3) for G in (G1,G2,G3))
    if delta.determinant():
        # Determinant is non-zero, no splitting
        return False, G1, G2, G3
    
    return True, G1, G2, G3

def get_isogeny_from_jacobian_two_kernel(kernel_generators, h):
    """
    Computes an isogeny from a Jacobian with a given 2-torsion kernel.

    :param kernel_generators: The kernel generators of the isogeny.
    :param h: The polynomial defining the hyperelliptic curve.
    :return: A tuple (codomain, isogeny_map).
    """
    is_split, G1, G2, G3 = is_jac_kernel_split(h, kernel_generators)
    if is_split:
        return FromJacToProd(G1, G2, G3)
    else:
        return FromJacToJac(G1, G2, G3)


def FromProdToJac(P, Q):
    """
    Constructs an isogeny from a product of two elliptic curves to a Jacobian of a hyperelliptic curve. Given two 2-torsion points generating the kernel on the product.

    :param P: A 2-torsion point on C x E. (CouplePoint)
    :param Q: A 2-torsion point on C x E. (CouplePoint)
    :return: A tuple (h, D11, D12, D21, D22, isogeny) where h is the new curve polynomial,
             (D11, D12) and (D21, D22) are image points, and isogeny is the map.
    """

    assert P.order() == 2 and Q.order() == 2, "Points must be of order 2. Got orders: {}, {}".format(P.order(), Q.order())

    P_c, P_e = P
    Q_c, Q_e = Q
    C = P_c.curve()
    Fp2 = C.base()
    Rx = PolynomialRing(Fp2, name="x")
    x = Rx.gen()

    a1, a2, a3 = P_c[0], Q_c[0], (P_c + Q_c)[0]
    b1, b2, b3 = P_e[0], Q_e[0], (P_e + Q_e)[0]
    # print(f"P_c: {P_c}, P_e: {P_e}")
    # print(f"Q_c: {Q_c}, Q_e: {Q_e}")
    # print(f"pairing: {P.weil_pairing(Q, 2)}")
    # print(f"a1: {a1}, a2: {a2}, a3: {a3}\nb1: {b1}, b2: {b2}, b3: {b3}")
    # Compute coefficients
    M = Matrix(Fp2, [
        [a1*b1, a1, b1],
        [a2*b2, a2, b2],
        [a3*b3, a3, b3]])
    R, S, T = M.inverse() * vector(Fp2, [1,1,1])
    RD = R * M.determinant()
    # print(f"R: {R}, S: {S}, T: {T}, RD: {RD}")
    da = (a1 - a2)*(a2 - a3)*(a3 - a1)
    db = (b1 - b2)*(b2 - b3)*(b3 - b1)

    s1, t1 = - da / RD, db / RD
    s2, t2 = -T/R, -S/R

    a1_t = (a1 - s2) / s1
    a2_t = (a2 - s2) / s1
    a3_t = (a3 - s2) / s1
    h = s1 * (x**2 - a1_t) * (x**2 - a2_t) * (x**2 - a3_t)

    H = HyperellipticCurve(h)
    J = H.jacobian()

    def isogeny(pair):
        # Argument members may be None to indicate the zero point.

        # The projection maps are:
        # H->C: (xC = s1/x²+s2, yC = s1 y)
        # so we compute Mumford coordinates of the divisor f^-1(P_c): a(x), y-b(x)
        Pc, P = pair
        if Pc:
            xPc, yPc = Pc.xy()
            JPc = J([s1 * x**2 + s2 - xPc, Rx(yPc / s1)])
        # Same for E
        # H->E: (xE = t1 x² + t2, yE = t1 y/x^3)
        if P:
            xP, yP = P.xy()
            JP = J([(xP - t2) * x**2 - t1, yP * x**3 / t1])
        if Pc and P:
            return JPc + JP
        if Pc:
            return JPc
        if P:
            return JP

    return h, isogeny

class RichelotCorr:
    """
    The Richelot correspondance between hyperelliptic
    curves y²=g1*g2*g3 and y²=h1*h2*h3=hnew(x)

    It is defined by equations:
        g1(x1) h1(x2) + g2(x1) h2(x2) = 0
    and y1 y2 = g1(x1) h1(x2) (x1 - x2)

    Given a divisor D in Mumford coordinates:
        U(x) = x^2 + u1 x + u0 = 0
        y = V(x) = v1 x + v0
    Let xa and xb be the symbolic roots of U.
    Let s, p by the sum and product (s=-u1, p=u0)

    Then on x-coordinates, the image of D is defined by equation:
        (g1(xa) h1(x) + g2(xa) h2(x))
      * (g1(xb) h1(x) + g2(xb) h2(x))
    which is a symmetric function of xa and xb.
    This is a non-reduced polynomial of degree 4.

    Write gred = g-U = g1*x + g0
    then gred(xa) gred(xb) = g1^2*p + g1*g0*s + g0^2
    and  g1red(xa) g2red(xb) + g1red(xb) g2red(xa)
       = 2 g11 g21 p + (g11*g20+g10*g21) s + 2 g10*g20

    On y-coordinates, the image of D is defined by equations:
           V(xa) y = Gred1(xa) h1(x) (xa - x)
        OR V(xb) y = Gred1(xb) h1(x) (xb - x)
    If we multiply:
    * y^2 has coefficient V(xa)V(xb)
    * y has coefficient h1(x) * (V(xa) Gred1(xb) (x-xb) + V(xb) Gred1(xa) (x-xa))
      (x-degree 3)
    * 1 has coefficient Gred1(xa) Gred1(xb) h1(x)^2 (x-xa)(x-xb)
                      = Gred1(xa) Gred1(xb) h1(x)^2 U(x)
      (x-degree 4)
    """
    def __init__(self, G1, G2, H1, H2, hnew):
        """
        Initializes the Richelot correspondence.

        :param G1: First quadratic factor of the domain curve polynomial.
        :param G2: Second quadratic factor of the domain curve polynomial.
        :param H1: First quadratic factor of the codomain curve polynomial.
        :param H2: Second quadratic factor of the codomain curve polynomial.
        :param hnew: The polynomial of the codomain curve.
        """
        assert G1[2].is_one() and G2[2].is_one()
        self.G1 = G1
        self.G2 = G2
        self.H1 = H1
        self.H11 = H1*H1
        self.H12 = H1*H2
        self.H22 = H2*H2
        self.hnew = hnew
        self.jacobian = HyperellipticCurve(hnew).jacobian()
        self.x = hnew.parent().gen()

    def map(self, P):
        "Computes (non-monic) Mumford coordinates for the image of P"
        U, V = P
        if not U[2].is_one():
            U = U / U[2]
    
        if U == self.G1 or U == self.G2:
            return self.jacobian(0)
        
        V = V  % U
        # Sum and product of (xa, xb)
        s, p = -U[1], U[0]
        # Compute X coordinates (non reduced, degree 4)
        g1red = self.G1 - U
        g2red = self.G2 - U
        assert g1red[2].is_zero() and g2red[2].is_zero()
        g11, g10 = g1red[1], g1red[0]
        g21, g20 = g2red[1], g2red[0]
        # see above
        Px = (g11*g11*p + g11*g10*s + g10*g10) * self.H11 \
           + (2*g11*g21*p + (g11*g20+g21*g10)*s + 2*g10*g20) * self.H12 \
           + (g21*g21*p + g21*g20*s + g20*g20) * self.H22

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
        Py1 = (2*v1*g11*p + v1*g10*s + v0*g11*s + 2*v0*g10)*self.x \
          - (v1*g11*s*p + 2*v1*g10*p + v0*g11*(s*s-2*p) + v0*g10*s)
        Py1 *= self.H1
        # coefficient of 1 is Gred1(xa) Gred1(xb) h1(x)^2 U(x)
        Py0 = self.H11 * U * (g11*g11*p + g11*g10*s + g10*g10)

        # Now reduce the divisor, and compute Cantor reduction.
        # Py2 * y^2 + Py1 * y + Py0 = 0
        # y = - (Py2 * hnew + Py0) / Py1
        _, Py1inv, _ = Py1.xgcd(Px)
        Py = (- Py1inv * (Py2 * self.hnew + Py0)) % Px
        assert Px.degree() == 4
        assert Py.degree() <= 3

        Dx = ((self.hnew - Py ** 2) // Px)
        Dy = (-Py) % Dx
        return self.jacobian([Dx, Dy])

def FromJacToJac(G1, G2, G3):
    """
    Computes a Richelot isogeny between two Jacobians of hyperelliptic curves.

    :param G1: First quadratic factor of the curve polynomial.
    :param G2: Second quadratic factor of the curve polynomial.
    :param G3: Third quadratic factor of the curve polynomial. Cannot be a linear combination of G1 and G2.
    :return: A tuple (h_new, phi) where phi is a function that computes the image of a RichelotJacobianPoint
      under the isogeny with kernel <P, Q> where P, Q correspond to the splittings G1 G2.
    """
    h = G1 * G2 * G3
    Rx = h.parent()
    x = Rx.gen()
    # G3, r3 = h.quo_rem(G1 * G2)
    # assert r3 == 0, f"h: {h} \n G1: {G1} \n G2: {G2} \n r3: {r3}"
    delta = Matrix(G.padded_list(3) for G in (G1,G2,G3))
    # H1 = 1/det (G2[1]*G3[0] - G2[0]*G3[1])
    #        +2x (G2[2]*G3[0] - G3[2]*G2[0])
    #        +x^2(G2[1]*G3[2] - G3[1]*G2[2])
    # The coefficients correspond to the inverse matrix of delta.
    delta = delta.inverse()
    H1 = -delta[0][0]*x**2 + 2*delta[1][0]*x - delta[2][0]
    H2 = -delta[0][1]*x**2 + 2*delta[1][1]*x - delta[2][1]
    H3 = -delta[0][2]*x**2 + 2*delta[1][2]*x - delta[2][2]

    hnew = H1*H2*H3
    R = RichelotCorr(G1, G2, H1, H2, hnew)

    def isogeny(P):
        output = R.map(P)
        # output must be a divisor on the new Jacobian
        assert  (hnew - Rx(output[1])**2) % Rx(output[0]) == 0
        return output

    return hnew, isogeny

def FromJacToProd(G1, G2, G3):
    """
    Construct the "split" isogeny from Jac(y^2 = G1*G2*G3)
    to a product of elliptic curves.
    """
    h = G1*G2*G3
    R = h.parent()
    Fp2 = R.base()
    x = R.gen()
    print(f"---------------- FromJacToProd ----------------")
    print(f" G's: G1: {G1}, G2: {G2}, G3: {G3}")
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

    def morphE1(x, y):
        # from y^2=p1 to y^2=p1norm
        return (H11*H21*H31*x, H11*H21*H31*y)
    def morphE2(x, y):
        # from y^2=p1 to y^2=p2norm
        return (H10*H20*H30*x, H10*H20*H30*y)
    # The morphisms are:
    # inverse homography:
    # H->H2: x, y => ((b-dx) / (x-a), y/(x-a)^3)
    # then H2->E1:(x,y) => (x^2,y)
    #   or H2->E2:(x,y) => (1/x^2,y/x^3)

    def isogeny(D):
        # HyperellipticCurve(h).jacobian()([D[0], D[1]])
        # To map a divisor, perform the change of coordinates
        # on Mumford coordinates
        U, V = D # Lets call the roots x_a, x_b
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
        # print(f"v1: {v1}, v0: {v0}, s: {s}, p: {p}, homography_needed: {homography_needed}")
        # print(f"CHEATING: DO NOT ALLOW IN PRODUCTION CODE")
        # roots = U_.roots(ring=Fp2.extension(2))
        # xa = roots[0][0]
        # xb = roots[1][0] if len(roots) > 1 else roots[0][0]
        # print(f"U Roots: {xa}, {xb}")
        # za_1 = xa**2
        # zb_1 = xb**2
        # print(f"za_1: {za_1}, zb_1: {zb_1}")
        # za_2 = 1/xa**2
        # zb_2 = 1/xb**2
        # print(f"za_2: {za_2}, zb_2: {zb_2}")
        # wa_1 = v1 * xa + v0
        # wb_1 = v1 * xb + v0
        # print(f"wa_1: {wa_1}, wb_1: {wb_1}")
        # wa_2 = (v1 * xa + v0) / xa**3
        # wb_2 = (v1 * xb + v0) / xb**3
        # print(f"wa_2: {wa_2}, wb_2: {wb_2}")


        assert p != 0 # nether x_a not x_b can be zero. Perhaps a change of coordinates is needed.

        # For E1: x -> x^2 := z, and y -> y := w
        # We need to compute the divisor on E1, corresponding to the sum of the image of (x_a, y_a) and (x_b, y_b), then sum it to get the image on E1.
        if s.is_zero():
            # This is the case where the divisor on E1 is the sum of two points with the same x-coordinate: x_a^2 = (-x_b)^2
            if v0.is_zero():
                print("E1: Special case: s=0, v0=0")
                # Then v0 = 0 <=> y_a = - y_b
                E1_image = E1(0)
            else:
                print("E1: Special case: s=0, v1=0")
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
                print("E1: Special case when v0 = 0 and s != 0")
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
                print("E2 Special case: s=0, v1=0")
                E2_image = E2(0)
            # w_a = w_b <=> v0 = 0
            elif v0.is_zero():
                print("E2 Special case: s=0, v0=0")
                # In this case w_a = w_b = v1 x_a / x_a^3 = v1 / x_a^2
                E2_image = 2 * E2(morphE2(1/-p, - v1/p)) # This may be incorrect?
            else:
                raise ValueError("Cannot have s=0 with v1, v0 both non-zero. This is mathematically impossible unless x_a = 0...") # this is because y_a^2 = y_b^2 implies (v1 x_a + v0)^2 = (v1 (-x_a) + v0)^2 => 4 v1 v0 x_a = 0
        else:
            C = v1*(s*s-2*p) + v0*s
            U2 = x**2 - (s*s-2*p)/p**2*x + 1/p**2 # The roots are 1/x_a^2 := z_a, 1/x_b^2 := z_b
            
            V21 = x**2 * C # this is z^2 * (x_ay_a + x_by_b)
            V20 = p2 + x**4 * (p*(v1**2*p + v1*v0*s + v0**2)) # w^2 + z^4(x_a x_b y_ay_b)
            print(f"V21: {V21}, V20: {V20}, U2: {U2}")

            # V21 * w = V20 modulo U2 (for (z_a, w_a) and (z_b, w_b))
            _d, V21inv, _ = V21.xgcd(U2)
            print(f"GCD of V21 and U2: {_d}")
            print(f"C: {C}")
            assert _d.is_one(), f"GCD not 1: {_d}\n s: {s}, p: {p}\n V21: {V21}, U2: {U2}"
            V2 = (V21inv * V20) % U2
            assert V2**2 % U2 == p2 % U2

            # Reduce coordinates
            U2red = (p2 - V2**2) // U2
            xP2 = -U2red[0] / U2red[1]
            yP2 = V2(xP2)
            assert yP2**2 == p2(xP2)
            E2_image = E2(morphE2(xP2, yP2))

        return CouplePoint(E1_image, E2_image)

    return (E1, E2), isogeny



# def FromJacToProd(G1, G2, G3):
#     print(f"---------------- FromJacToProd ----------------")
#     h = G1*G2*G3
#     R = h.parent()
#     Fp2 = R.base()
#     x = R.gen()
#     print(f"h: {h}")
#     print(f"G1: {G1}, G2: {G2}, G3: {G3}")
#     print(f"G1 roots: {G1.roots()},\n G2 roots: {G2.roots()},\n G3 roots: {G3.roots()}")

#     M = Matrix(G.padded_list(3) for G in (G1,G2,G3))
#     # Find homography
#     need_homography = True
#     u, v, w = M.right_kernel().gen()
#     if u.is_zero():
#         H11, H21, H31 = M.column(2)
#         H10, H20, H30 = M.column(0)
#         need_homography = False
#     else:
#         d = u/2
#         (ad, _), (b, _) = (x**2 - v*x + w*d/2).roots()
#         a = ad / d

#         # Apply transform G(x) -> G((a*x+b)/(x+d))*(x+d)^2
#         # The coefficients of x^2 are M * (1, a, a^2)
#         # The coefficients of 1 are M * (d^2, b*d, b^2)
#         H11, H21, H31 = M * vector([1, a, a*a])
#         H10, H20, H30 = M * vector([d*d, b*d, b*b])
#         assert G1((a*x+b)/(x+d))*(x+d)**2 == H11*x**2+H10

#     # h2 = (H11*x**2+H10)*(H21*x**2+H20)*(H31*x**2+H30)
#     # H2 = HyperellipticCurve(h2)

#     p1 = (H11*x+H10)*(H21*x+H20)*(H31*x+H30)
#     p2 = (H11+H10*x)*(H21+H20*x)*(H31+H30*x)
#     # We will need to map to actual elliptic curve
#     p1norm = (x + H10*H21*H31)*(x + H20*H11*H31)*(x + H30*H11*H21)
#     p2norm = (x + H11*H20*H30)*(x + H21*H10*H30)*(x + H31*H10*H20)
#     E1 = EllipticCurve([0, p1norm[2], 0, p1norm[1], p1norm[0]])
#     E2 = EllipticCurve([0, p2norm[2], 0, p2norm[1], p2norm[0]])

#     def morphE1(x, y):
#         # from y^2=p1 to y^2=p1norm
#         return (H11*H21*H31*x, H11*H21*H31*y)
#     def morphE2(x, y):
#         # from y^2=p2 to y^2=p2norm
#         return (H10*H20*H30*x, H10*H20*H30*y)
#     # The morphisms are:
#     # inverse homography:
#     # H->H2: x, y => ((b-dx) / (x-a), y/(x-a)^3)
#     # then H2->E1:(x,y) => (x^2,y)
#     #   or H2->E2:(x,y) => (1/x^2,y/x^3)

#     def isogeny(D):
#         # HyperellipticCurve(h).jacobian()(D)
#         # To map a divisor, perform the change of coordinates
#         # on Mumford coordinates
#         U, V = D
#         print(f"U: {U}, V: {V}")
#         print(f"U roots: {U.roots()}")
#         # apply homography (if applied before)
#         # y = v1 x + v0 =>
#         if need_homography:
#             U_ = U[0] * (x+d)**2 + U[1]*(a*x+b)*(x+d) + U[2]*(a*x+b)**2
#             V_ = V[0] * (x+d)**3 + V[1]*(a*x+b)*(x+d)**2
#             V_ = V_ % U_
#         else:
#             U_, V_ = U, V
#         v1, v0 = V_[1], V_[0]
#         # prepare symmetric functions
#         s = - U_[1] / U_[2] # s = x1 + x2
#         p = U_[0] / U_[2] # p = x1 * x2
#         print(f"s: {s}\np: {p}")
#         # compute Mumford coordinates on E1
#         # Points x1, x2 map to x1^2, x2^2
#         U1 = x**2 - (s*s - 2*p)*x + p**2
#         print(f"U1 roots: {U1.roots()}")
#         # y = v1 x + v0 becomes (y - v0)^2 = v1^2 x^2
#         # so 2v0 y-v0^2 = p1 - v1^2 xH^2 = p1 - v1^2 xE1
#         # V1 = (p1 - v1**2 * x + v0**2) / (2*v0)
#         ### Added by Nic Swanson for Debugging
#         E1_image = None
#         if v0.is_zero() and s.is_zero():
#             E1_image = E1(0)
#         else:
#             if v0.is_zero():
#                 V1 = (v1 / s) * (x + p)
#             else:
#                 V1 = (p1 - v1**2 * x + v0**2) / (2*v0)
#             V1 = V1 % U1
#             U1red = (p1 - V1**2) // U1
#             xP1 = -U1red[0] / U1red[1]
#             yP1 = V1(xP1)
#             assert yP1**2 == p1(xP1)
#             E1_image = E1(morphE1(xP1, yP1))
#         ###
#         # Reduce Mumford coordinates to get a E1 point
        
#         # Same for E2
#         # Points x1, x2 map to 1/x1^2, 1/x2^2
#         U2 = x**2 - (s*s-2*p)/p**2*x + 1/p**2
#         print(f"U2 roots: {U2.roots()}")
#         E2_image = None
#         if v1.is_zero() and s.is_zero():
#             E2_image = E2(0)
#         else: 
#             # yE = y1/x1^3, xE = 1/x1^2
#             # means yE = y1 x1 xE^2
#             # (yE - y1 x1 xE^2)(yE - y2 x2 xE^2) = 0
#             # p2 - yE (x1 y1 + x2 y2) xE^2 + (x1 y1 x2 y2 xE^4) = 0
#             V21 = x**2 * (v1 * (s*s-2*p) + v0*s)
#             V20 = p2 + x**4 * (p*(v1**2*p + v1*v0*s + v0**2))
#             # V21 * y = V20
#             _, V21inv, _ = V21.xgcd(U2)
#             V2 = (V21inv * V20) % U2
#             assert V2**2 % U2 == p2 % U2
#             # Reduce coordinates
#             U2red = (p2 - V2**2) // U2
#             xP2 = -U2red[0] / U2red[1]
#             yP2 = V2(xP2)
#             assert yP2**2 == p2(xP2)
#             E2_image = E2(morphE2(xP2, yP2))

#         return CouplePoint(E1_image, E2_image)

#     return (E1, E2), isogeny