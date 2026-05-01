"""newton_birkhoff.py -- Newton-polygon / Birkhoff formal-solution
extractor for degree-2 PCF generating functions.

For a continued fraction with a_n = alpha1 n + alpha0,
b_n = beta2 n^2 + beta1 n + beta0, the partial-denominator sequence
Q_n satisfies the three-term linear recurrence

    Q_n = b_n Q_{n-1} + a_n Q_{n-2},     Q_0 = 1, Q_1 = b_1.

The generating function f(z) = sum_{n >= 0} Q_n z^n satisfies the
order-2 linear ODE (with theta = z d/dz)

    L f = 1,  where  L = 1
                       - z [ beta2 (theta+1)^2 + beta1 (theta+1) + beta0 ]
                       - z^2 [ alpha1 (theta+2) + alpha0 ].

Module API:
    build_op(a_coeffs, b_coeffs)    -> dict with operator coeffs c_{ij}
    newton(op)                       -> {slopes, multiplicities,
                                          char_eqs, char_roots}
    formal_solve(op, Kmax, dps)      -> two formal solutions
                                          f_i(u) = exp(mu_i/u) u^{rho_i}
                                          (1 + sum_{k>=1} a_{i,k} u^k)
                                        with z = u^2.
    borel_singularities(formal, dps) -> nearest singularities of the
                                          Gevrey-1 (in u) Borel transform
                                          of (1 + sum a_{i,k} u^k).
    detect_painleve(op, formal, dps) -> Painleve-class fingerprint
                                          (heuristic: compare indicial
                                          exponents and leading formal
                                          coefficients to standard P-III,
                                          P-V tau invariants).

Conventions for Newton polygon: in the (i, j) plane (i = z-power,
j = theta-power), the operator's points are
    (0, 0): coeff +1
    (1, 0): coeff -b(1) = -(beta2 + beta1 + beta0)
    (1, 1): coeff -(2*beta2 + beta1)
    (1, 2): coeff -beta2
    (2, 0): coeff -(2*alpha1 + alpha0)
    (2, 1): coeff -alpha1
The single nonzero slope of the lower-left convex hull at z=0 is 1/2,
edge (0,0)-(1,2). Characteristic equation along the edge is
    1 - beta2 (c^2/4) = 0  =>  c = +/- 2 / sqrt(beta2).
This means f admits a formal pair of solutions
    f_+/-(z) = exp(+/- 2 z^{-1/2} / sqrt(beta2)) z^{rho_+/-}
              (1 + O(z^{1/2})).
Setting z = u^2 turns this into a Gevrey-1 formal expansion at u=0.
"""

from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Dict, List, Tuple
import mpmath as mp


# ---------------------------------------------------------------------------

@dataclass
class PCFOp:
    """Linear ODE operator L = sum c_{ij} z^i theta^j  with theta = z d/dz."""
    alpha1: mp.mpf | int
    alpha0: mp.mpf | int
    beta2: mp.mpf | int
    beta1: mp.mpf | int
    beta0: mp.mpf | int
    family: str = "?"

    def points(self) -> Dict[Tuple[int, int], mp.mpf]:
        b1 = self.beta2 + self.beta1 + self.beta0
        return {
            (0, 0): mp.mpf(1),
            (1, 0): -mp.mpf(b1),
            (1, 1): -mp.mpf(2 * self.beta2 + self.beta1),
            (1, 2): -mp.mpf(self.beta2),
            (2, 0): -mp.mpf(2 * self.alpha1 + self.alpha0),
            (2, 1): -mp.mpf(self.alpha1),
        }


def build_op(family: str,
             alpha1: int, alpha0: int,
             beta2: int, beta1: int, beta0: int) -> PCFOp:
    return PCFOp(alpha1=alpha1, alpha0=alpha0,
                 beta2=beta2, beta1=beta1, beta0=beta0, family=family)


# ---------------------------------------------------------------------------
# Newton polygon at z = 0
# ---------------------------------------------------------------------------

def newton(op: PCFOp) -> Dict:
    """Compute the Newton polygon at z=0 of L.

    For our class of operators the polygon's lower-left convex hull
    has vertices at (0,0) and (1,2), giving a single slope-(1/2) edge
    of multiplicity 2 (in the j-direction).  The characteristic equation
    along that edge is

        c_{0,0} + c_{1,2} (c/2)^2 = 0   (after substituting theta -> -c/(2 u)
                                          and z -> u^2, leading-order WKB).

    The "+/- (c/2)^2" comes from theta acting on exp(c/u): theta(exp(c/u))
    = z * d/dz exp(c/(z^{1/2})) = z * (-c/(2 z^{3/2})) exp(c/z^{1/2})
    = -(c/2) z^{-1/2} exp(c/z^{1/2}) = -(c/2) (1/u) exp(c/u). So
    theta^2 acts as (c/2)^2 (1/u^2) exp(c/u).  Substituting back into L
    and balancing on the edge (i - j/2 = 0): coefficient at (0,0)
    contributes c_{0,0} = 1, at (1,2) contributes c_{1,2} (c/2)^2
    z^{1-1} = (-beta2)(c^2/4).

    So the characteristic polynomial is

        chi(c) = 1 - (beta2/4) c^2

    with roots c = +/- 2 / sqrt(beta2).
    """
    pts = op.points()
    c00 = pts[(0, 0)]
    c12 = pts[(1, 2)]
    # chi(c) = c00 + c12 * (c/2)^2
    # = 1 - (beta2/4) c^2 = 0
    c_squared = -mp.mpf(4) * c00 / c12  # = 4 / beta2
    c_root = mp.sqrt(c_squared)
    # indicial exponents at z=0 from the (1,1) and (0,0) Newton sub-edge:
    # subleading term gives the rho exponent.  Using ansatz
    #   f = exp(c/u) u^rho (1 + a1 u + ...),  z = u^2
    # and substituting into L f = 0 (homogeneous), the next-order matching
    # determines rho; see formal_solve for the explicit recurrence.
    return {
        "slopes": [{"slope": "1/2",
                    "edge": "(0,0) -- (1,2)",
                    "multiplicity": 2,
                    "char_poly": f"1 - ({op.beta2}/4) c^2",
                    "char_roots_c": [c_root, -c_root],
                    "abs_action_xi0": c_root}],
        "gevrey_class_in_z": 2,
        "gevrey_class_in_u_with_z_eq_u_sq": 1,
    }


# ---------------------------------------------------------------------------
# Formal solution f_i(u) = exp(c_i/u) u^{rho_i} (1 + sum_{k>=1} a_k u^k)
# in the variable u = sqrt(z).  We substitute the ansatz into the
# homogeneous ODE  L f = 0, expand in powers of u, and solve for
# (rho_i, a_1, a_2, ...).
# ---------------------------------------------------------------------------

def _theta_pow_on_ansatz(j: int, c: mp.mpf, rho_var, a_series: List[mp.mpf],
                         u: mp.mpf, K: int) -> List[mp.mpf]:
    """Compute the u-series of theta^j applied to
       g(u) = exp(c/u) u^{rho_var} (1 + a_1 u + ... + a_K u^K)
    truncated to u^K times g0 := exp(c/u) u^rho_var.

    theta = z d/dz = (u/2) d/du   (since z = u^2 -> dz = 2u du,
                                    z d/dz = u^2 * (1/(2u)) d/du = (u/2) d/du).

    Apply (u/2 d/du) to g.  Let g(u) = exp(c/u) u^rho_var * S(u) where
    S(u) = 1 + a_1 u + a_2 u^2 + ... .  Then
        d/du g = exp(c/u) [(-c/u^2) u^rho S + rho u^{rho-1} S + u^rho S']
        u/2 * d/du g = exp(c/u) u^rho [ (-c/(2u)) S + (rho/2) S + (u/2) S' ]
                     = (1/u) * (exp(c/u) u^rho) * [ (-c/2) S + (rho/2) u S + (u^2/2) S' ]

    So if we factor out g0 := exp(c/u) u^rho, then
        theta(g)/g0 = (1/u) * [ (-c/2) S(u) + (rho/2) u S(u) + (u^2/2) S'(u) ]
    Apply repeatedly for theta^j.  We carry the prefactor 1/u^j explicitly
    by tracking series in u of 'P_j(u) = u^j * (theta^j g)/g0'; then
    P_j is a polynomial in u of degree <= K + j (we truncate at u^K relative
    to the leading 1/u^j-pole structure).

    Return: coefficients of P_j as list of length K+1 (powers u^0 ... u^K),
    where the contribution to (theta^j g)/g0 is P_j(u) / u^j.
    """
    # Represent series S(u) = 1 + a_1 u + ... up to u^K
    S = [mp.mpf(0)] * (K + 1)
    S[0] = mp.mpf(1)
    for k in range(1, K + 1):
        S[k] = a_series[k - 1] if k - 1 < len(a_series) else mp.mpf(0)
    # P_0 = S; we want P_j of length K+j+1 but truncated.
    # Recurrence: P_{j+1}(u) = (-c/2) P_j(u) + (1/2) (theta_eff P_j)(u)
    # where the recurrence on P_j comes from
    #   (theta g)/g0 = (1/u) [ (-c/2) S + (rho/2) u S + (u^2/2) S' ]
    # generalised: applying theta to (1/u^j) Q(u) gives a similar result.
    # Actually let R_j(u) := (theta^j g)/g0.  Then
    #   theta R_j = (-c/(2u)) R_j + (1/2) [u d/du + ?] R_j
    # but R_j has its own u-structure.  Let's recompute carefully:
    # if R(u) := (theta^j g)/g0 = u^{-j} * Q(u)  (Q polynomial-ish), then
    #   theta(R g0) = theta_acting_on_R_times_g0 = [theta R + (theta_g0/g0)R] g0
    # theta_g0/g0 = theta(exp(c/u) u^rho)/[exp(c/u) u^rho]
    #             = (u/2) d/du log(exp(c/u) u^rho) = (u/2) ( -c/u^2 + rho/u )
    #             = -c/(2u) + rho/2
    # And theta R = (u/2) d/du R.  So
    #   theta(R g0)/g0 = theta R + (-c/(2u) + rho/2) R
    # That is the operator we iterate to get R_{j+1} = theta R_j +
    #     [-c/(2u) + rho/2] R_j.
    # Multiplying by u^{j+1} (writing R_j = u^{-j} Q_j):
    #   u^{j+1} R_{j+1} = u^{j+1} [(u/2) d/du)(u^{-j} Q_j)] +
    #                     u^{j+1} [-c/(2u) + rho/2] u^{-j} Q_j
    #                  = (u^{j+1}/2) [ -j u^{-j-1} Q_j + u^{-j} Q_j' ]
    #                       + u [-c/(2u) + rho/2] Q_j
    #                  = (1/2) [ -j Q_j + u Q_j' ] + (-c/2 + (rho/2) u) Q_j
    #                  = (-c/2 - j/2) Q_j + (rho/2) u Q_j + (u/2) Q_j'
    # So the recurrence on Q_j(u) = u^j * (theta^j g)/g0 is:
    #     Q_0(u) = S(u)
    #     Q_{j+1}(u) = (-c/2 - j/2) Q_j(u) + (rho_var/2) u Q_j(u)
    #                  + (u/2) Q_j'(u).
    #
    # All operations preserve a polynomial truncation to degree K
    # (because differentiation lowers degree, multiplication by u raises by 1
    #  but we truncate).
    Q = list(S)  # length K+1
    for jj in range(j):
        Qp = [mp.mpf(0)] * (K + 1)
        # term A = (-c/2 - jj/2) * Q
        A = mp.mpf(-c) / 2 - mp.mpf(jj) / 2
        for k in range(K + 1):
            Qp[k] += A * Q[k]
        # term B = (rho/2) u * Q  -> shifts up by 1
        for k in range(K):
            Qp[k + 1] += (mp.mpf(rho_var) / 2) * Q[k]
        # term C = (u/2) * Q'  : Q'(u) has coefficients k * Q[k] at u^{k-1}
        for k in range(1, K + 1):
            # u/2 * (k * Q[k] u^{k-1}) = (k/2) Q[k] u^k
            Qp[k] += (mp.mpf(k) / 2) * Q[k]
        Q = Qp
    return Q


def formal_solve(op: PCFOp, K: int, dps: int = 200,
                 sign: int = +1) -> Dict:
    """Compute the formal solution of L f = 0 of the form
       f(u) = exp(sign * c0 / u) * u^{rho} * (1 + a_1 u + ... + a_K u^K)
    where c0 = 2 / sqrt(beta2).

    Strategy: treat rho and a_1, ..., a_K as unknowns.  Substitute
    into the homogeneous part L_h f = - z^2 (alpha1 (theta+2) + alpha0) f
    - z (beta2 (theta+1)^2 + beta1 (theta+1) + beta0) f + 0*f
    (the +1 part of L is the inhomogeneity for the particular solution,
    not the homogeneous formal solutions).

    We need the homogeneous operator L_h such that L = 1 - L_h.  Then
    L_h f = f.  Equivalently, formal solutions of L_h - 1 = 0.

    Substitute into L f and balance powers of u (note z = u^2):
        z^i theta^j f / g0 = u^{2i} * u^{-j} Q_j(u)
        where Q_j is the polynomial computed above.

    L f / g0 = sum_{(i,j)} c_{ij} u^{2i - j} Q_j(u)
            = sum_k (coefficient of u^k) u^k.

    The leading-order balance is at the smallest 2i-j (most negative);
    setting that to zero is the characteristic equation, which is
    automatic for c = +/- 2/sqrt(beta2).  The next-order balance
    determines rho.  Successive orders determine a_1, a_2, ...

    Returns the formal coefficients and meta.
    """
    mp.mp.dps = dps
    beta2 = mp.mpf(op.beta2)
    c0 = mp.mpf(2) / mp.sqrt(beta2)
    c = sign * c0

    pts = op.points()  # {(i,j): coeff}

    # Determine rho by leading-subleading balance.
    # Use sympy-free symbolic: rho is determined linearly.
    # Build expansion of L f / g0 in powers of u, treating rho and a_1..a_K
    # symbolically (we use a numerical solve: rho is a single unknown
    # whose coefficient is linear; a_k are determined recursively).

    # Practical approach: solve order by order.  For a given rho candidate
    # (initially symbolic), the equations at order u^{-2} are automatic,
    # at u^{-1} give a linear equation in rho, at u^0 give a_1 in terms
    # of rho, etc.  We solve rho first, then a_1, a_2, ...

    # Step 1: solve for rho.
    # Need coefficients of u^{-2}, u^{-1} in L f / g0.
    # Pre-compute Q_j(u; rho, a_*) at arbitrary K-truncation, but with
    # a_1, ..., a_K = 0 first to get the leading rho equation (since
    # a_k contribute only to higher orders).
    #
    # The lowest order coefficient comes from (i, j) with 2i - j minimum.
    # For our points: 2i - j values:
    #   (0,0): 0
    #   (1,0): 2
    #   (1,1): 1
    #   (1,2): 0
    #   (2,0): 4
    #   (2,1): 3
    # min = 0, achieved at (0,0) and (1,2).  Good, that's the slope-1/2
    # edge.  The contribution to u^0 from these vanishes by the
    # characteristic equation.  But Q_j(u) has u-expansion, so each
    # (i,j) contributes to multiple u-powers.  The lowest u-power
    # contributed by (i,j) is 2i - j (from the constant term Q_j(0)).
    # So u^{-something} appears only if 2i - j < 0.  In our case all
    # 2i - j >= 0, so no negative u-powers.  Good.
    #
    # Therefore the formal series f / g0 has no negative u-powers in
    # L f, and the leading equation is at u^0.  At u^0:
    #
    # contributions from each (i, j) at u^{2i - j} with multiplier
    # Q_j(0) = constant term:
    #   (0,0): coeff * Q_0(0) = 1 * 1 = 1
    #   (1,2): -beta2 * Q_2(0)
    # Q_2(0) is computed via the recurrence with S(u) = 1 (a_k=0):
    #   Q_0 = [1, 0, 0, ...]
    #   Q_1 = (-c/2 - 0/2)*Q_0 + (rho/2) u Q_0 + (u/2) Q_0' (Q_0' = 0)
    #       = [-c/2, rho/2, 0, ...]
    #   Q_2 = (-c/2 - 1/2)*Q_1 + (rho/2) u Q_1 + (u/2) Q_1'
    #       Q_1' = [rho/2, 0, 0, ...] -> (u/2) Q_1' contributes (rho/4) u
    #       Q_2[0] = (-c/2 - 1/2) * (-c/2) = c^2/4 + c/4
    #       Q_2[1] = (-c/2 - 1/2)*(rho/2) + (rho/2)*(-c/2) + (rho/4)
    #             = -rho c/4 - rho/4 - rho c/4 + rho/4
    #             = -rho c / 2
    # So u^0-coefficient of L f at a_*=0:
    #   1 + (-beta2) * (c^2/4 + c/4)  + ... contributions from (1,1),(1,0)
    # but (1,1) gives 2i-j = 1 contribution at u^1, not u^0.
    # (1,0) gives 2i-j = 2 contribution at u^2.
    # (2,0) and (2,1) give >= 3.
    # So at u^0:
    #   1 + (-beta2)(c^2/4 + c/4) = 0
    #   beta2 c^2 / 4 = 1  -> c^2 = 4/beta2  (CHARACTERISTIC EQUATION).
    # plus a c-linear residual: -beta2 * c/4.  At c = +/- 2/sqrt(beta2),
    # this is -beta2 * (+/- 1/sqrt(beta2)) / 2 = -/+ sqrt(beta2)/2.
    # Hmm that's nonzero.  So the u^0-eq is NOT automatic!
    # That means I need to think again about what's at u^0.
    #
    # Wait: I computed (theta^j g)/g0 = (1/u^j) Q_j(u).  So z^i theta^j f
    # contributes u^{2i} / u^j * Q_j(u) = u^{2i - j} Q_j(u).  For (1,2):
    # u^{2-2} Q_2(u) = Q_2(u).  Q_2(0) = c^2/4 + c/4 (from above).
    # Hmm.  So Q_2(0) is not just c^2/4 -- there's a c/4 from the recursion.
    # Let me recheck: when I derived Q_{j+1} = (-c/2 - j/2)Q_j + ... I had
    # the (j/2) coming from differentiating u^{-j} d/du.  Let me re-derive.
    #
    # theta = (u/2) d/du.  Let R_j = (theta^j g)/g0 = u^{-j} Q_j(u).
    # theta R_j = (u/2) d/du (u^{-j} Q_j) = (u/2)[ -j u^{-j-1} Q_j + u^{-j} Q_j']
    #          = (1/2)(-j u^{-j} Q_j + u^{-j+1} Q_j')
    # Hmm wait, (u/2)(-j u^{-j-1} Q_j) = (-j/2) u^{-j} Q_j.
    # (u/2)(u^{-j} Q_j') = (1/2) u^{-j+1} Q_j'.
    # So theta R_j = (-j/2) u^{-j} Q_j + (1/2) u^{-j+1} Q_j'.
    #
    # And theta(R_j g0) / g0 = theta R_j + (theta g0 / g0) R_j.
    # Now theta g0 / g0 with g0 = exp(c/u) u^rho:
    # log g0 = c/u + rho log u.  d/du(log g0) = -c/u^2 + rho/u.
    # theta(log g0) = (u/2)(-c/u^2 + rho/u) = -c/(2u) + rho/2.
    # So theta(g0)/g0 = -c/(2u) + rho/2.
    # Therefore theta(R_j g0)/g0 = (-j/2 - c/(2u) + rho/2) u^{-j} Q_j + (1/2) u^{-j+1} Q_j'.
    # We want this = u^{-(j+1)} Q_{j+1}, i.e. multiply by u^{j+1}:
    # Q_{j+1} = (-j/2 - c/(2u) + rho/2) u Q_j + (1/2) u^2 Q_j'
    #       = (-j/2 + rho/2) u Q_j  - (c/2) Q_j + (u^2/2) Q_j'
    # So Q_{j+1}(u) = -(c/2) Q_j(u) + ((rho - j)/2) u Q_j(u) + (u^2/2) Q_j'(u).
    #
    # OK this is different from what I had before.  Let me redo Q_2 carefully:
    # Q_0 = [1, 0, 0, ...]
    # Q_1 = -(c/2) Q_0 + ((rho-0)/2) u Q_0 + (u^2/2) Q_0' = -(c/2) + (rho/2) u
    # Q_1 = [-c/2, rho/2, 0, ...]
    # Q_1' = [rho/2, 0, 0, ...]
    # Q_2 = -(c/2) Q_1 + ((rho-1)/2) u Q_1 + (u^2/2) Q_1'
    # Q_2 const: -(c/2)*(-c/2) = c^2/4
    # Q_2 u^1: -(c/2)*(rho/2) + ((rho-1)/2) * (-c/2) + 0
    #        = -rho c/4 + (-rho-1+1)... let me just compute:
    #        = -(c/2)(rho/2) + ((rho-1)/2)(-c/2)
    #        = -c rho/4 + (-c/4)(rho - 1)
    #        = -c rho/4 - c rho/4 + c/4
    #        = -c rho/2 + c/4
    # Q_2 u^2: -(c/2)*0 + ((rho-1)/2)*(rho/2) + (1/2) * (rho/2)*1*u**0 from u^2/2 * Q_1'
    #   wait Q_1'[k] for k=0 is rho/2 (coefficient of u^0 in dQ_1/du).
    #   (u^2/2) * Q_1' has u^2 * (rho/2) at lowest order, so u^2 coefficient is rho/4
    # Q_2[2] = ((rho-1)/2)(rho/2) + rho/4 = rho(rho-1)/4 + rho/4 = rho^2/4 - rho/4 + rho/4 = rho^2/4.
    # Hmm interesting -- Q_2[2] = rho^2/4 = (theta^2 acting on u^rho at leading).
    #
    # So Q_2(0) = c^2/4 (NOT c^2/4 + c/4 as I had before).  Good.
    # That means at u^0, the L f equation gives:
    #    sum over (i,j) with 2i-j=0 of c_{ij} Q_j(0)
    #    = c_{00}*Q_0(0) + c_{12}*Q_2(0)
    #    = 1*1 + (-beta2)*(c^2/4)
    #    = 1 - beta2 c^2/4
    # = 0 iff c^2 = 4/beta2. GOOD.  Characteristic equation correct.
    #
    # Now at u^1: contributions from
    #  (0,0): Q_0[1] = 0
    #  (1,1): u^{2-1}=u^1 * Q_1(u); coef of u^0 in Q_1 = -c/2; so u^1-contribution
    #         from (1,1) is c_{11} * (-c/2) = -(2 beta2 + beta1) * (-c/2)
    #                                       = (2 beta2 + beta1) * c / 2
    #  (1,2): Q_2[1] = -c rho/2 + c/4; * c_{12} = -beta2 * (-c rho/2 + c/4)
    #                                          = beta2 c rho/2 - beta2 c/4
    # Setting total to 0:
    #   (2 beta2 + beta1) c/2 + beta2 c rho/2 - beta2 c/4 = 0
    #   c [ (2 beta2 + beta1)/2 + beta2 rho/2 - beta2/4 ] = 0
    # since c != 0:
    #   beta2 rho / 2 = -(2 beta2 + beta1)/2 + beta2/4 = -beta2 - beta1/2 + beta2/4
    #                = -3 beta2/4 - beta1/2
    #   rho = -3/2 - beta1 / beta2.
    # So rho = -3/2 - beta1/beta2 (independent of sign of c).
    return _formal_solve_numeric(op, K, dps, sign)


def _formal_solve_numeric(op: PCFOp, K: int, dps: int,
                          sign: int) -> Dict:
    """Numerical formal solution: substitute ansatz with unknown rho, a_k
    and solve order by order from u^0 = char eq up to u^{K+2}."""
    mp.mp.dps = dps
    beta2 = mp.mpf(op.beta2)
    beta1 = mp.mpf(op.beta1)
    beta0 = mp.mpf(op.beta0)
    alpha1 = mp.mpf(op.alpha1)
    alpha0 = mp.mpf(op.alpha0)

    c = mp.mpf(sign) * mp.mpf(2) / mp.sqrt(beta2)
    rho = mp.mpf(-3) / 2 - beta1 / beta2
    pts = {
        (0, 0): mp.mpf(1),
        (1, 0): -(beta2 + beta1 + beta0),
        (1, 1): -(2 * beta2 + beta1),
        (1, 2): -beta2,
        (2, 0): -(2 * alpha1 + alpha0),
        (2, 1): -alpha1,
    }

    # We will compute coefficients a_1, a_2, ..., a_K of the formal series
    # S(u) = 1 + a_1 u + a_2 u^2 + ... .  At each order m of u, the
    # equation is linear in a_m (and depends on a_1, ..., a_{m-1}).
    #
    # Plan: maintain Q_j(u) as a list of length K+3 (in u-powers), where
    # each entry is a polynomial in a_1..a_K but really we'll do the
    # numeric substitution iteratively: solve a_1, then update Q_j and
    # solve a_2, etc.  Easier: treat S as the unknown polynomial, build
    # the operator's action on S to obtain a triangular system.
    #
    # Concretely: define for each j the operator T_j taking a polynomial
    # S (in u) to Q_j(u) via the recurrence.  Then L f / g0 = sum c_{ij}
    # u^{2i-j} Q_j(u).  We want this = 0 mod u^{K+1}.  Each a_k appears
    # linearly (because Q_j is linear in S, hence linear in {a_k}).
    #
    # So we compute the operator's matrix M (rows: u-power equation, cols:
    # unknown a_k), with column 0 being the constant residual (no a's).
    # Then a_k = - residual_k / diag_k for the diagonal, with back-sub
    # for off-diagonal.  In our case the structure is upper-triangular
    # if we look at how a_m enters the u^m equation.

    # Build Q_j operator action symbolically as a list-of-list-of-coeffs:
    # Q_j[k][m] = coefficient of a_m in Q_j at u^k, with a_0 := 1.
    Klen = K + 4  # working precision in u-powers
    a_dim = K + 1  # a_0 = 1, a_1, ..., a_K

    def zero_mat():
        return [[mp.mpf(0)] * a_dim for _ in range(Klen)]

    # Q_0 = S = sum_{m=0..K} a_m u^m
    Q = []
    Q0 = zero_mat()
    for m in range(a_dim):
        if m < Klen:
            Q0[m][m] = mp.mpf(1)
    Q.append(Q0)

    # Recurrence Q_{j+1}[k][m] from Q_j:
    # Q_{j+1}(u) = -(c/2) Q_j(u) + ((rho - j)/2) u Q_j(u) + (u^2/2) Q_j''_in_u? No,
    # recurrence above: Q_{j+1}(u) = -(c/2) Q_j(u) + ((rho - j)/2) u Q_j(u) + (u^2/2) Q_j'(u).
    for j in range(2):
        Qj = Q[-1]
        Qj1 = zero_mat()
        for k in range(Klen):
            for m in range(a_dim):
                # term -(c/2) Q_j[k][m]
                Qj1[k][m] += -(c / 2) * Qj[k][m]
            if k - 1 >= 0:
                for m in range(a_dim):
                    # ((rho - j)/2) * u * Q_j  shifts k-1 -> k
                    Qj1[k][m] += ((rho - mp.mpf(j)) / 2) * Qj[k - 1][m]
            if k - 1 >= 0:
                # (u^2/2) Q_j'(u): Q_j' has coefficient k * Q_j[k] at u^{k-1};
                # u^2/2 * (k_old * Q_j[k_old]) lives at u^{k_old + 1}.
                # So Qj1[k] gets (k_old/2) Q_j[k_old] where k_old = k - 1.
                k_old = k - 1
                if k_old >= 0:
                    for m in range(a_dim):
                        Qj1[k][m] += (mp.mpf(k_old) / 2) * Qj[k_old][m]
        Q.append(Qj1)

    # Now build E[k][m] = coefficient of a_m in (L f / g0)[u^k]
    E = zero_mat()
    for (i, jj), cij in pts.items():
        shift = 2 * i - jj
        Qj = Q[jj]
        for k_src in range(Klen):
            k_dst = k_src + shift
            if 0 <= k_dst < Klen:
                for m in range(a_dim):
                    E[k_dst][m] += cij * Qj[k_src][m]

    # Verify char eq at k=0: E[0][0] should be 0 (no a's contribute at u^0;
    # they appear only at u^>=1).
    char_residual = E[0][0]
    if abs(char_residual) > mp.mpf(10) ** (-dps + 10):
        raise RuntimeError(f"Char eq residual not zero: {char_residual}")

    # Solve for a_1, a_2, ..., a_K.  Structure: at row u^k the
    # "newest" coefficient appearing with NON-ZERO multiplier is a_{k-1}
    # (the column a_k vanishes at row k because the (0,0) and (1,2)
    # contributions cancel via the characteristic equation, the same
    # cancellation that makes c = +/- 2/sqrt(beta2) the slope-1/2
    # exponent).  Therefore equation at u^k determines a_{k-1}, with
    # diagonal E[k][k-1].
    a = [mp.mpf(0)] * a_dim
    a[0] = mp.mpf(1)
    for k in range(1, K + 1):
        # equation at row u^{k+0} -> determines a_k via diagonal E[k+1][k]
        # i.e. use row r = k + 1 to solve a_k.
        r = k + 1
        if r >= Klen:
            raise RuntimeError(f"Klen={Klen} too small for K={K}")
        rhs = E[r][0]  # constant part (a_0 = 1 already absorbed)
        for m in range(1, k):
            rhs += E[r][m] * a[m]
        diag = E[r][k]
        if abs(diag) < mp.mpf(10) ** (-dps + 20):
            raise RuntimeError(f"Shifted diagonal vanishes at r={r},col={k}: {diag}")
        a[k] = -rhs / diag

    return {
        "family": op.family,
        "sign": sign,
        "c": c,
        "rho": rho,
        "a": a,
        "K": K,
        "dps": dps,
        "pts": {f"{i},{j}": str(v) for (i, j), v in pts.items()},
    }


# ---------------------------------------------------------------------------
# Borel-plane radius extraction (Gevrey-1-in-u)
# ---------------------------------------------------------------------------

def borel_radius(formal: Dict, n_use: int = None) -> Dict:
    """Estimate the radius |w*| of the nearest singularity of the
    Gevrey-1 Borel transform B[S](w) = sum a_k w^k / k! using
    Domb-Sykes (ratio of consecutive coefficients of S/k! style) or
    direct |a_k / k!|^{1/k} -> 1/|w*|.

    We compute the limit of |b_k|^{1/k} where b_k = a_k / k! for large k.
    We expect 1/|w*| to converge as k grows.  Returns the inverse, |w*|.
    """
    a = formal["a"]
    K = formal["K"]
    if n_use is None:
        n_use = max(8, K // 2)
    # b_k = a_k / k!
    bs = []
    for k in range(1, K + 1):
        ak = a[k]
        kfact = mp.factorial(k)
        bs.append(ak / kfact)
    # |b_k|^{1/k} for k in [n_use, K]
    radii = []
    for k in range(n_use, K + 1):
        bk = bs[k - 1]
        if abs(bk) == 0:
            continue
        radii.append((k, mp.power(abs(bk), mp.mpf(1) / k)))
    # tail mean (last quartile) is more stable
    tail = radii[-max(1, len(radii) // 4):]
    inv_w = sum(r for _, r in tail) / mp.mpf(len(tail))
    w_star = mp.mpf(1) / inv_w if inv_w != 0 else mp.inf
    return {
        "n_used_range": [n_use, K] if radii else None,
        "radii_sample": [(k, str(r)) for (k, r) in radii[-8:]],
        "w_star_abs_estimate": w_star,
    }


# ---------------------------------------------------------------------------
# Painleve fingerprint (heuristic).
# We compare (xi_0, rho, leading 4 formal coefficients) across families.
# A genuine P-III(D_6) match would require monodromy data; here we only
# compare the *invariants captured by the formal solution* (xi_0, rho,
# normalised a_1..a_4).  Two families share these invariants iff they
# admit the same isomonodromic reduction.
# ---------------------------------------------------------------------------

def painleve_fingerprint(formal: Dict, ref_formal: Dict, threshold: int = 15) -> Dict:
    """Compare formal-solution invariants of two families.

    Returns digit-of-agreement on (xi_0, rho, a_1..a_4) and a verdict
    "match" if all >= threshold digits, "no-match" otherwise.
    """
    digits = {}
    for key in ("c", "rho"):
        v1, v2 = formal[key], ref_formal[key]
        if v2 == 0:
            d = mp.mpf("inf") if v1 == 0 else mp.mpf(0)
        else:
            err = abs(v1 - v2) / abs(v2) if abs(v2) > 0 else abs(v1 - v2)
            d = -mp.log10(err) if err > 0 else mp.mpf(2 * formal["dps"])
        digits[key] = float(d)
    for k in range(1, 5):
        v1 = formal["a"][k]
        v2 = ref_formal["a"][k]
        if abs(v2) > 0:
            err = abs(v1 - v2) / abs(v2)
        else:
            err = abs(v1 - v2)
        d = -mp.log10(err) if err > 0 else mp.mpf(2 * formal["dps"])
        digits[f"a_{k}"] = float(d)
    min_digits = min(digits.values())
    verdict = "match" if min_digits >= threshold else "no-match"
    return {"digits": digits, "min_digits": min_digits, "verdict": verdict,
            "threshold": threshold}
