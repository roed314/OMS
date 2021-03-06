
#        Copyright (C) 2012 Robert Pollack <rpollack@math.bu.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

import operator

from sage.structure.element import ModuleElement
from sage.matrix.matrix_integer_2x2 import MatrixSpace_ZZ_2x2
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.misc.cachefunc import cached_method
from sage.rings.padics.factory import Qp
from sage.rings.polynomial.all import PolynomialRing
from sage.rings.padics.padic_generic import pAdicGeneric
from sage.rings.arith import next_prime
from sage.rings.infinity import infinity
from sage.misc.misc import verbose
from sage.rings.padics.precision_error import PrecisionError

from sage.categories.action import Action
from fund_domain import Id
from manin_map import ManinMap, M2Z
from padic_lseries import pAdicLseries
from sigma0 import Sigma0

minusproj = [1,0,0,-1]


class PSModSymAction(Action):
    def __init__(self, actor, MSspace):
        Action.__init__(self, actor, MSspace, False, operator.mul)

    def _call_(self, sym, g):
        return sym.__class__(sym._map * g, sym.parent(), construct=True)

class PSModularSymbolElement(ModuleElement):
    def __init__(self, map_data, parent, construct=False):
        ModuleElement.__init__(self, parent)
        if construct:
            self._map = map_data
        else:
            self._map = ManinMap(parent._coefficients, parent._source, map_data)

    def _repr_(self):
        r"""
        Returns the print representation of the symbol.

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi._repr_()
            'Modular symbol of level 11 with values in Sym^0 Q^2'
        """
        return "Modular symbol of level %s with values in %s"%(self.parent().level(),self.parent().coefficient_module())
    
    def dict(self):
        r"""
        Returns dictionary on the modular symbol self, where keys are generators and values are the corresponding values of self on generators

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi.dict()
            {[1 0]
            [0 1]: -1/5, [ 0 -1]
            [ 1  3]: 3/2, [-1 -1]
            [ 3  2]: -1/2}
        """
        D = {}
        for g in self.parent().source().gens():
            D[g] = self._map[g]
        return D

    def weight(self):
        r"""
        Returns the weight of this Pollack-Stevens modular symbol.

        This is `k-2`, where `k` is the usual notion of weight for modular
        forms!

        EXAMPLES::
            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi.weight()
            0

        """
        return self.parent().weight()

    def values(self):
        r"""
        Returns the values of the symbol self on our chosen generators (generators are listed in self.dict().keys())

        EXAMPLES::

             sage: E = EllipticCurve('11a')
             sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
             sage: phi = ps_modsym_from_elliptic_curve(E)
             sage: phi.values()
             [-1/5, 3/2, -1/2]
             sage: phi.dict().keys()
             [
             [1 0]  [ 0 -1]  [-1 -1]
             [0 1], [ 1  3], [ 3  2]
             ]
             sage: phi.values() == phi.dict().values()
             True
        """
        return [self._map[g] for g in self.parent().source().gens()]

    def _normalize(self):
        """
        Normalizes all of the values of the symbol self

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi._normalize()
            Modular symbol of level 11 with values in Sym^0 Q^2
            sage: phi._normalize().values()
            [-1/5, 3/2, -1/2]
        """
        for val in self._map:
            val.normalize()
        return self

    def __cmp__(self, other):
        """
        Checks if self == other. Here self and other have the same parent.

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi == phi
            True
            sage: phi == 2*phi
            False
            sage: psi = ps_modsym_from_elliptic_curve(EllipticCurve('37a'))
            sage: psi == phi
            False
        """
        gens = self.parent().source().gens()
        for g in gens:
            c = cmp(self._map[g], other._map[g])
            if c: return c
        return 0

    def _add_(self, right):
        """
        Returns self + right

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
            [-1/5, 3/2, -1/2]
            sage: phi + phi
            Modular symbol of level 11 with values in Sym^0 Q^2
            sage: (phi + phi).values()
            [-2/5, 3, -1]
        """
        return self.__class__(self._map + right._map, self.parent(), construct=True)

    def _lmul_(self, right):
        """
        Returns self * right

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
            [-1/5, 3/2, -1/2]
            sage: 2*phi
            Modular symbol of level 11 with values in Sym^0 Q^2
            sage: (2*phi).values()
            [-2/5, 3, -1]
        """
        return self.__class__(self._map * right, self.parent(), construct=True)

    def _rmul_(self, right):
        """
        Returns self * right

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
            [-1/5, 3/2, -1/2]
            sage: phi*2
            Modular symbol of level 11 with values in Sym^0 Q^2
            sage: (phi*2).values()
            [-2/5, 3, -1]
        """
        return self.__class__(self._map * right, self.parent(), construct=True)

    def _sub_(self, right):
        """
        Returns self - right

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
            [-1/5, 3/2, -1/2]
            sage: phi - phi
            Modular symbol of level 11 with values in Sym^0 Q^2
            sage: (phi - phi).values()
            [0, 0, 0]
        """
        return self.__class__(self._map - right._map, self.parent(), construct=True)

    def _get_prime(self, p=None, alpha = None, allow_none=False):
        """
        Combines a prime specified by the user with the prime from the parent.

        INPUT:

        - ``p`` -- an integer or None (default None); if specified
          needs to match the prime of the parent.

        - ``alpha`` -- an element or None (default None); if p-adic
          can contribute a prime.

        - ``allow_none`` -- boolean (default False); whether to allow
          no prime to be specified.

        OUTPUT:

        - a prime or None.  If ``allow_none`` is False then a
          ValueError will be raised rather than returning None if no
          prime can be determined.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
            sage: D = Distributions(0, 5, 10);  M = PSModularSymbols(Gamma0(5), coefficients=D)
            sage: f = M(1); f._get_prime()
            5
            sage: f._get_prime(5)
            5
            sage: f._get_prime(7)
            Traceback (most recent call last):
            ...
            ValueError: inconsistent prime
            sage: f._get_prime(alpha=Qp(5)(1))
            5
            sage: D = Symk(0);  M = PSModularSymbols(Gamma0(2), coefficients=D)
            sage: f = M(1); f._get_prime(allow_none=True) is None
            True
            sage: f._get_prime(alpha=Qp(7)(1))
            7
            sage: f._get_prime(7,alpha=Qp(7)(1))
            7
            sage: f._get_prime()
            Traceback (most recent call last):
            ...
            ValueError: you must specify a prime
        """
        pp = self.parent().prime()
        ppp = ((alpha is not None) and hasattr(alpha.parent(),'prime') and alpha.parent().prime()) or None
        p = ZZ(p) or pp or ppp
        if not p:
            if not allow_none:
                raise ValueError("you must specify a prime")
        elif (pp and p != pp) or (ppp and p != ppp):
            raise ValueError("inconsistent prime")
        return p

    def plus_part(self):
        r"""
        Returns the plus part of self -- i.e. self + self | [1,0,0,-1].

        Note that we haven't divided by 2.  Is this a problem?

        OUTPUT:

        - self + self | [1,0,0,-1]

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
            [-1/5, 3/2, -1/2]
            sage: (phi.plus_part()+phi.minus_part()) == 2 * phi
            True
        """
        S0N = Sigma0(self.parent().level())
        return self + self * S0N(minusproj)

    def minus_part(self):
        r"""
        Returns the minus part of self -- i.e. self - self | [1,0,0,-1]

        Note that we haven't divided by 2.  Is this a problem?

        OUTPUT:

        - self - self | [1,0,0,-1]

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
            [-1/5, 3/2, -1/2]
            sage: (phi.plus_part()+phi.minus_part()) == phi * 2
            True
        """
        S0N = Sigma0(self.parent().level())
        return self - self * S0N(minusproj)

    def hecke(self, ell, algorithm="prep"):
        r"""
        Returns self | `T_{\ell}` by making use of the precomputations in
        self.prep_hecke()

        INPUT:

        - ``ell`` -- a prime

        - ``algorithm`` -- a string, either 'prep' (default) or
          'naive'

        OUTPUT:

        - The image of this element under the hecke operator
          `T_{\ell}`

        ALGORITHMS:

        - If ``algorithm == 'prep'``, precomputes a list of matrices
          that only depend on the level, then uses them to speed up
          the action.

        - If ``algorithm == 'naive'``, just acts by the matrices
          defining the Hecke operator.  That is, it computes
          sum_a self | [1,a,0,ell] + self | [ell,0,0,1],
          the last term occurring only if the level is prime to ell.

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
            [-1/5, 3/2, -1/2]
            sage: phi.hecke(2) == phi * E.ap(2)
            True
            sage: phi.hecke(3) == phi * E.ap(3)
            True
            sage: phi.hecke(5) == phi * E.ap(5)
            True
            sage: phi.hecke(101) == phi * E.ap(101)
            True

            sage: all([phi.hecke(p, algorithm='naive') == phi * E.ap(p) for p in [2,3,5,101]])
            True
        """
        return self.__class__(self._map.hecke(ell, algorithm), self.parent(), construct=True)

    def valuation(self, p=None):
        r"""
        Returns the valuation of self at `p`.

        Here the valuation is the minimum of the valuations of the values of self.

        INPUT:

        - ``p`` - prime

        OUTPUT:

        - The valuation of self at `p`

        EXAMPLES::

           sage: E = EllipticCurve('11a')
           sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
           sage: phi = ps_modsym_from_elliptic_curve(E)
           sage: phi.values()
           [-1/5, 3/2, -1/2]
           sage: phi.valuation(2)
           -1
           sage: phi.valuation(3)
           0
           sage: phi.valuation(5)
           -1
           sage: phi.valuation(7)
           0
           sage: phi.valuation()
           Traceback (most recent call last):
           ...
           ValueError: you must specify a prime

           sage: phi2 = phi.lift(11, M=2)
           sage: phi2.valuation()
           0
           sage: phi2.valuation(3)
           Traceback (most recent call last):
           ...
           ValueError: inconsistent prime
           sage: phi2.valuation(11)
           0
        """
        q = self._get_prime(p)
        return min([val.valuation(q) for val in self._map])

    def diagonal_valuation(self, p):
        """
        Retuns the minimum of the diagonal valuation on the values of self

        INPUT:

        - ``p`` -- a positive integral prime

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi.values()
            [-1/5, 3/2, -1/2]
            sage: phi.diagonal_valuation(2)
            -1
            sage: phi.diagonal_valuation(3)
            0
            sage: phi.diagonal_valuation(5)
            -1
            sage: phi.diagonal_valuation(7)
            0
        """
        return min([val.diagonal_valuation(p) for val in self._map])

    @cached_method
    def is_Tq_eigensymbol(self,q,p=None,M=None):
        r"""
        Determines if self is an eigenvector for `T_q` modulo `p^M`

        INPUT:

        - ``q`` -- prime of the Hecke operator

        - ``p`` -- prime we are working modulo

        - ``M`` -- degree of accuracy of approximation

        OUTPUT:

        - True/False

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi.values()
            [-1/5, 3/2, -1/2]
            sage: phi_ord = phi.p_stabilize(p = 3, ap = E.ap(3), M = 10, ordinary = True)
            sage: phi_ord.is_Tq_eigensymbol(2,3,10)
            True
            sage: phi_ord.is_Tq_eigensymbol(2,3,100)
            False
            sage: phi_ord.is_Tq_eigensymbol(2,3,1000)
            False
            sage: phi_ord.is_Tq_eigensymbol(3,3,10)
            True
            sage: phi_ord.is_Tq_eigensymbol(3,3,100)
            False
        """
        try:
            aq = self.Tq_eigenvalue(q, p, M)
            return True
        except ValueError:
            return False

    # what happens if a cached method raises an error?  Is it recomputed each time?
    @cached_method
    def Tq_eigenvalue(self, q, p=None, M=None, check=True):
        r"""
        Eigenvalue of `T_q` modulo `p^M`

        INPUT:

        - ``q`` -- prime of the Hecke operator

        - ``p`` -- prime we are working modulo (default: None)

        - ``M`` -- degree of accuracy of approximation (default: None)

        - ``check`` -- check that `self` is an eigensymbol

        OUTPUT:

        - Constant `c` such that `self|T_q - c * self` has valuation greater than
          or equal to `M` (if it exists), otherwise raises ValueError

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi.values()
            [-1/5, 3/2, -1/2]
            sage: phi_ord = phi.p_stabilize(p = 3, ap = E.ap(3), M = 10, ordinary = True)
            sage: phi_ord.Tq_eigenvalue(2,3,10) + 2
            O(3^10)

            sage: phi_ord.Tq_eigenvalue(3,3,10)
            2 + 3^2 + 2*3^3 + 2*3^4 + 2*3^6 + 3^8 + 2*3^9 + O(3^10)
            sage: phi_ord.Tq_eigenvalue(3,3,100)
            Traceback (most recent call last):
            ...
            ValueError: not a scalar multiple
        """
        qhecke = self.hecke(q)
        gens = self.parent().source().gens()
        if p is None:
            p = self.parent().prime()
        i = 0
        g = gens[i]
        verbose("Computing eigenvalue")
        while self._map[g].is_zero(p, M):
            if not qhecke._map[g].is_zero(p, M):
                raise ValueError("not a scalar multiple")
            i += 1
            try:
                g = gens[i]
            except IndexError:
                raise ValueError("self is zero")
        aq = self._map[g].find_scalar(qhecke._map[g], p, M, check)
        verbose("Found eigenvalues of %s"%(aq))
        if check:
            verbose("Checking that this is actually an eigensymbol")
            if p is None or M is None:
                for g in gens[1:]:
                    if qhecke._map[g] != aq * self._map[g]:
                        raise ValueError("not a scalar multiple")
            elif (qhecke - aq * self).valuation(p) < M:
                raise ValueError("not a scalar multiple")
        return aq

    def is_ordinary(self,p=None,P=None):
        r"""
        Returns true if the p-th eigenvalue is a p-adic unit.

        INPUT:
        
        - ``p`` - a positive integral prime, or None (default None)
        - ``P`` - a prime of the base ring above `p`, or None. This is ignored
          unless the base ring is a number field.

        OUTPUT:

        - True/False

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: E = EllipticCurve('11a1')
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi.is_ordinary(2)
            False
            sage: E.ap(2)
            -2
            sage: phi.is_ordinary(3)
            True
            sage: E.ap(3)
            -1
            sage: phip = phi.p_stabilize(3,20)
            sage: phip.is_ordinary()
            True

        A number field example. Here there are multiple primes above `p`, and
        `\phi` is ordinary at one but not the other.::

            sage: f = Newforms(32, 8, names='a')[1]
            sage: K = f.hecke_eigenvalue_field()
            sage: a = f[3]
            sage: phi = f.PS_modular_symbol()
            sage: phi.is_ordinary(K.ideal(3, 1/16*a + 3/2))
            False
            sage: phi.is_ordinary(K.ideal(3, 1/16*a + 5/2))
            True
            sage: phi.is_ordinary(3)
            Traceback (most recent call last):
            ...
            TypeError: P must be an ideal

        """
        # q is the prime below p, if base is a number field; q = p otherwise
        if p == None:
            if self.parent().prime() == 0:
                raise ValueError("need to specify a prime")
            q = p = self.parent().prime()
        elif p in ZZ:
            q = p
        else:
            q = p.smallest_integer()
        if not q.is_prime():
            raise ValueError("p is not prime")
        if (self.parent().prime() != q) and (self.parent().prime() != 0):
            raise ValueError("prime does not match coefficient module's prime")
        aq = self.Tq_eigenvalue(q)
        return aq.valuation(p) == 0

    def _consistency_check(self):
        """
        Check that the map really does satisfy the Manin relations loop (for debugging).
        The two and three torsion relations are checked and it is checked that the symbol
        adds up correctly around the fundamental domain

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: E = EllipticCurve('37a1')
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi._consistency_check()
            This modular symbol satisfies the manin relations

        """

        f = self._map
        MR = self._map._manin
        ## Test two torsion relations
        for g in MR.reps_with_two_torsion():
            gamg = MR.two_torsion_matrix(g)
            if not (f[g]*gamg + f[g]).is_zero():
                raise ValueError("Two torsion relation failed with",g)

        ## Test three torsion relations
        for g in MR.reps_with_three_torsion():
            gamg = MR.three_torsion_matrix(g)
            if not (f[g]*(gamg**2) + f[g]*gamg + f[g]).is_zero():
                raise ValueError("Three torsion relation failed with",g)

        ## Test that the symbol adds to 0 around the boundary of the fundamental domain
        t = self.parent().coefficient_module().zero_element()
        for g in MR.gens()[1:]:
            if (not g in MR.reps_with_two_torsion()) and (not g in MR.reps_with_three_torsion()):
                t += f[g] * MR.gammas[g] - f[g]
            else:
                if g in MR.reps_with_two_torsion():
                    t -= f[g] 
                else:
                    t -= f[g]
                    
        id = MR.gens()[0]
        if f[id]*MR.gammas[id] - f[id] != -t:
            print t
            print f[id]*MR.gammas[id] - f[id]
            raise ValueError("Does not add up correctly around loop")

        print "This modular symbol satisfies the manin relations"

class PSModularSymbolElement_symk(PSModularSymbolElement):
    def _find_M(self, M):
        """
        Determines `M` from user input. ?????

        INPUT:

        - ``M`` -- an integer at least 2 or None.  If None, sets `M` to
          be one more than the precision cap of the parent (the
          minimum amount of lifting).

        OUTPUT:

        - An updated ``M``.

        EXAMPLES::

            sage: pass
        """
        if M is None:
            M = ZZ(20)
        elif M <= 1:
            raise ValueError("M must be at least 2")
        else:
            M = ZZ(M)
        return M

    def _find_alpha(self, p, k, M=None, ap=None, new_base_ring=None, ordinary=True, check=True, find_extraprec=True):
        r"""
        Finds `alpha`, a `U_p` eigenvalue, which is found as a root of
        the polynomial `x^2 - ap * x + p^(k+1)*chi(p)`.

        INPUT:

        - ``p`` -- prime

        - ``k`` -- Pollack-Stevens weight

        - ``M`` -- precision (default: None) of `Q_p`

        - ``ap`` -- Hecke eigenvalue at p (default: None)

        - ``new_base_ring`` -- field of definition of `alpha` (default: None)

        - ``ordinary`` -- True if the prime is ordinary (default: True)

        - ``check`` -- check to see if the prime is ordinary (default: True)

        - ``find_extraprec`` -- setting this to True finds extra precision (default: True)

        OUTPUT:

        The output is a tuple (`alpha`, `new_base_ring`, `newM`, `eisenloss`,`q`,`aq`), with

        - ``alpha`` --  `U_p` eigenvalue

        - ``new_base_ring`` -- field of definition of `alpha` with precision at least `newM`

        - ``newM`` -- new precision

        - ``eisenloss`` -- loss of precision

        - ``q`` -- a prime not equal to p which was used to find extra precision

        - ``aq`` -- the Hecke eigenvalue `aq` corresponding to `q`

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: E = EllipticCurve('11a')
            sage: p = 5
            sage: M = 10
            sage: k = 0
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi._find_alpha(p,k,M)
            (1 + 4*5 + 3*5^2 + 2*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 3*5^7 + 2*5^8 + 3*5^9 + 3*5^10 + 3*5^12 + O(5^13), 5-adic Field with capped relative precision 13, 12, 1, 2, -2)
        """
        if ap is None:
            ap = self.Tq_eigenvalue(p, check=check)
        if check and ap.valuation(p) > 0:
            raise ValueError("p is not ordinary")

        chi = self._map._codomain._character
        if chi is not None:
            eps = chi(p)
        else:
            eps = 1
        poly = PolynomialRing(ap.parent(), 'x')([p**(k+1) * eps, -ap, 1])
        if new_base_ring is None:
            # These should actually be completions of disc.parent()
            if p == 2:
                # is this the right precision adjustment for p=2?
                new_base_ring = Qp(2, M+1)
            else:
                new_base_ring = Qp(p, M)
            set_padicbase = True
        else:
            set_padicbase = False
        try:
            verbose("finding alpha: rooting %s in %s"%(poly, new_base_ring))
            (v0,e0),(v1,e1) = poly.roots(new_base_ring)
        except (TypeError, ValueError):
            raise ValueError("new base ring must contain a root of x^2 - ap * x + p^(k+1)")
        if v0.valuation(p) > 0:
            v0, v1 = v1, v0
        if ordinary:
            alpha = v0
        else:
            alpha = v1
        if find_extraprec:
            newM, eisenloss, q, aq = self._find_extraprec(p, M, alpha, check)
        else:
            newM, eisenloss, q, aq = M, None, None, None
        if set_padicbase:
            # We want to ensure that the relative precision of alpha and (alpha-1) are both at least *newM*,
            # where newM is obtained from self._find_extraprec
            prec_cap = None
            verbose("testing prec_rel: newM = %s, alpha = %s"%(newM, alpha), level=2)
            if alpha.precision_relative() < newM:
                prec_cap = newM + alpha.valuation(p) + (1 if p == 2 else 0)
            if ordinary:
                a1val = (alpha - 1).valuation(p)
                verbose("a1val = %s"%a1val, level=2)
                if a1val > 0 and ap != 1 + p**(k+1): # if ap = 1 + p**(k+1) then alpha = 1 and we need to give up.
                    if prec_cap is None:
                        prec_cap = newM + a1val + (1 if p == 2 else 0)
                    else:
                        prec_cap = max(prec_cap, newM + a1val + (1 if p == 2 else 0))
            verbose("prec_cap = %s"%(prec_cap), level=2)
            if prec_cap is not None:
                new_base_ring = Qp(p, prec_cap)
                return self._find_alpha(p=p, k=k, M=M, ap=ap, new_base_ring=new_base_ring, ordinary=ordinary, check=False, find_extraprec=find_extraprec)
        return alpha, new_base_ring, newM, eisenloss, q, aq

    def p_stabilize(self, p=None, M=None, alpha=None, ap=None, new_base_ring=None, ordinary=True, check=True):
        r"""

        Returns the `p`-stablization of self to level `N*p` on which `U_p` acts by `alpha`.

        Note that since `alpha` is `p`-adic, the resulting symbol
        is just an approximation to the true `p`-stabilization
        (depending on how well `alpha` is approximated).

        INPUT:

        - ``p`` -- prime not dividing the level of self

        - ``M`` -- precision of `Q_p`

        - ``alpha`` -- `U_p` eigenvalue

        - ``ap`` -- Hecke eigenvalue

        - ``new_base_ring`` -- change of base ring

        OUTPUT:

        A modular symbol with the same Hecke eigenvalues as
        self away from `p` and eigenvalue `alpha` at `p`.
        The eigenvalue `alpha` depends on the parameter `ordinary`.

        If ordinary == True: the unique modular symbol of level
        `N*p` with the same Hecke eigenvalues as self away from
        `p` and unit eigenvalue at `p`; else  the unique modular
        symbol of level `N*p` with the same Hecke eigenvalues as
        self away from `p` and non-unit eigenvalue at `p`.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: E = EllipticCurve('11a')
            sage: p = 5
            sage: prec = 4
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phis = phi.p_stabilize(p,M = prec)
            sage: phis
            Modular symbol of level 55 with values in Sym^0 Q_5^2
            sage: phis.hecke(7) == phis*E.ap(7)
            True
            sage: phis.hecke(5) == phis*E.ap(5)
            False
            sage: phis.hecke(3) == phis*E.ap(3)
            True
            sage: phis.Tq_eigenvalue(5)
            1 + 4*5 + 3*5^2 + 2*5^3 + O(5^4)
            sage: phis = phi.p_stabilize(p,M = prec,ordinary=False)
            sage: phis.Tq_eigenvalue(5)
            5 + 5^2 + 2*5^3 + O(5^4)

        A complicated example (with nontrivial character)::

            sage: chi = DirichletGroup(24)([-1, -1, -1])
            sage: f = Newforms(chi,names='a')[0]
            sage: phi = f.PS_modular_symbol()
            sage: phi11, h11 = phi.completions(11,5)[0]
            sage: phi11s = phi11.p_stabilize()
            sage: phi11s.is_Tq_eigensymbol(11)         
            True
        """
        if check:
            p = self._get_prime(p, alpha)
        k = self.parent().weight()
        M = self._find_M(M)
        verbose("p stabilizing: M = %s"%M, level=2)
        if alpha is None:
            alpha, new_base_ring, newM, eisenloss, q, aq = self._find_alpha(p, k, M, ap, new_base_ring, ordinary, check, False)
        else:
            if new_base_ring is None:
                new_base_ring = alpha.parent()
            if check:
                if ap is None:
                    ap = self.base_ring()(alpha + p**(k+1)/alpha)
                elif alpha**2 - ap * alpha + p**(k+1) != 0:
                    raise ValueError("alpha must be a root of x^2 - a_p*x + p^(k+1)")
                if self.hecke(p) != ap * self:
                    raise ValueError("alpha must be a root of x^2 - a_p*x + p^(k+1)")
        verbose("found alpha = %s"%(alpha))
        V = self.parent()._p_stabilize_parent_space(p, new_base_ring)
        return self.__class__(self._map.p_stabilize(p, alpha, V), V, construct=True)

    def completions(self, p, M):
        r"""
        If `K` is the base_ring of self, this function takes all maps
        `K-->Q_p` and applies them to self return a list of
        (modular symbol,map: `K-->Q_p`) as map varies over all such maps.

        .. NOTE::

            This only returns all completions when `p` splits completely in `K`

        INPUT:

        - ``p`` -- prime

        - ``M`` -- precision

        OUTPUT:

        - A list of tuples (modular symbol,map: `K-->Q_p`) as map varies over all such maps

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_simple_modsym_space
            sage: D = ModularSymbols(67,2,1).cuspidal_submodule().new_subspace().decomposition()[1]
            sage: f = ps_modsym_from_simple_modsym_space(D)
            sage: S = f.completions(41,10); S
            [(Modular symbol of level 67 with values in Sym^0 Q_41^2, Ring morphism:
              From: Number Field in alpha with defining polynomial x^2 + 3*x + 1
              To:   41-adic Field with capped relative precision 10
              Defn: alpha |--> 5 + 22*41 + 19*41^2 + 10*41^3 + 28*41^4 + 22*41^5 + 9*41^6 + 25*41^7 + 40*41^8 + 8*41^9 + O(41^10)), (Modular symbol of level 67 with values in Sym^0 Q_41^2, Ring morphism:
              From: Number Field in alpha with defining polynomial x^2 + 3*x + 1
              To:   41-adic Field with capped relative precision 10
              Defn: alpha |--> 33 + 18*41 + 21*41^2 + 30*41^3 + 12*41^4 + 18*41^5 + 31*41^6 + 15*41^7 + 32*41^9 + O(41^10))]
            sage: TestSuite(S[0][0]).run()
        """
        K = self.base_ring()
        R = Qp(p,M+10)['x']
        x = R.gen()
        if K == QQ:
            f = x-1
        else:
            f = K.defining_polynomial()
        v = R(f).roots()
        if len(v) == 0:
            L = Qp(p,M).extension(f,names='a')
            a = L.gen()
            V = self.parent().change_ring(L)
            Dist = V.coefficient_module()
            psi = K.hom([K.gen()],L)
            embedded_sym = self.parent().element_class(self._map.apply(psi,codomain=Dist, to_moments=True),V, construct=True)
            ans = [embedded_sym,psi]
            return ans
        else:
            roots = [r[0] for r in v]
            ans = []
            V = self.parent().change_ring(Qp(p, M))
            Dist = V.coefficient_module()
            for r in roots:
                psi = K.hom([r],Qp(p,M))
                embedded_sym = self.parent().element_class(self._map.apply(psi, codomain=Dist, to_moments=True), V, construct=True)
                ans.append((embedded_sym,psi))
            return ans

    def lift(self, p=None, M=None, alpha=None, new_base_ring=None, algorithm='stevens', eigensymbol=False, check=True):
        r"""
        Returns a (`p`-adic) overconvergent modular symbol with
        `M` moments which lifts self up to an Eisenstein error

        Here the Eisenstein error is a symbol whose system of Hecke
        eigenvalues equals `ell+1` for `T_ell` when `ell`
        does not divide `Np` and 1 for `U_q` when `q` divides `Np`.

        INPUT:

        - ``p`` -- prime

        - ``M`` -- integer equal to the number of moments

        - ``alpha`` -- `U_p` eigenvalue

        - ``new_base_ring`` -- change of base ring

        - ``algorithm`` -- 'stevens' or 'greenberg' (default 'stevens')

        - ``eigensymbol`` -- if True, lifts to Hecke eigensymbol (self must be a `p`-ordinary eigensymbol)

        (Note: ``eigensymbol = True`` does *not* just indicate to the code that
        self is an eigensymbol; it solves a wholly different problem, lifting
        an eigensymbol to an eigensymbol.)

        OUTPUT:

        An overconvergent modular symbol whose specialization equals self, up
        to some Eisenstein error if ``eigensymbol`` is False. If ``eigensymbol
        = True`` then the output will be an overconvergent Hecke eigensymbol
        (and it will lift the input exactly, the Eisenstein error disappears).

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: E = EllipticCurve('11a')
            sage: f = ps_modsym_from_elliptic_curve(E)
            sage: g = f.lift(11,4,algorithm='stevens',eigensymbol=True)
            sage: g.is_Tq_eigensymbol(2)
            True
            sage: g.Tq_eigenvalue(3)
            10 + 10*11 + 10*11^2 + 10*11^3 + O(11^4)
            sage: g.Tq_eigenvalue(11)
            1 + O(11^4)

        We check that lifting and then specializing gives back the original symbol::

            sage: g.specialize() == f
            True
        """
        if p is None:
            p = self.parent().prime()
            if p == 0:
                raise ValueError("must specify a prime")
        elif (self.parent().prime() != 0) and p != self.parent().prime():
            raise ValueError("inconsistent prime")
        if M is None:
            M = self.parent().precision_cap() + 1
###  I don't understand this.  This might only make sense in weight 2.  Probably need a bound
###  on M related to the weight.
        elif M <= 1:
            raise ValueError("M must be at least 2")
        else:
            M = ZZ(M)
        if new_base_ring is None:
            if isinstance(self.parent().base_ring(), pAdicGeneric):
                new_base_ring = self.parent().base_ring()
            else:
                # We may need extra precision in solving the difference equation
                extraprec = (M-1).exact_log(p)
                # should eventually be a completion
                new_base_ring = Qp(p, M+extraprec)
        if algorithm is None:
            raise NotImplementedError
        if algorithm == 'stevens':
            if eigensymbol:
                # We need some extra precision due to the fact that solving
                # the difference equation can give denominators.
                if alpha is None:
                    alpha = self.Tq_eigenvalue(p, check=check)
                newM, eisenloss, q, aq = self._find_extraprec(p, M, alpha, check)
                return self._lift_to_OMS_eigen(p, M, new_base_ring, alpha, newM, eisenloss, q, aq, check)
            else:
                return self._lift_to_OMS(p, M, new_base_ring, check)
        elif algorithm == 'greenberg':
            return self._lift_greenberg(p, M, new_base_ring, check)
        else:
            raise ValueError("algorithm %s not recognized" % algorithm)

    
    def _lift_greenberg(self, p, M, new_base_ring=None, check=False):
        """
        This is the Greenberg algorithm for lifting a modular eigensymbol to
        an overconvergent modular symbol. One first lifts to any set of numbers
        (not necessarily satifying the Manin relations). Then one applies the U_p,
        and normalizes this result to get a lift satisfying the manin relations.
       
        
        INPUT:
        
        - ``p`` -- prime

        - ``M`` -- integer equal to the number of moments

        - ``new_base_ring`` -- new base ring

        - ``check`` -- THIS IS CURRENTLY NOT USED IN THE CODE!
        
        OUTPUT: 
        
        - an overconvergent modular symbol lifting the symbol that was input
        
        EXAMPLES:: 

            sage: E = EllipticCurve('11a')
            sage: phi = E.PS_modular_symbol()
            sage: Phi = phi.lift(11,5,algorithm='greenberg')
            sage: Phi2 = phi.lift(11,5,algorithm='stevens',eigensymbol=True)
            sage: Phi == Phi2
            True

        An example in higher weight::

            sage: f = Newforms(7, 4)[0].PS_modular_symbol()
            sage: fs = f.p_stabilize(5)
            sage: FsG = fs.lift(M=6, eigensymbol=True,algorithm='greenberg') 
            sage: FsG.values()[0]
            (2 + 5 + 3*5^2 + 4*5^3 + O(5^6), O(5^5), 2*5 + 3*5^2 + O(5^4), O(5^3), 5 + O(5^2), O(5))
            sage: FsS = fs.lift(M=6, eigensymbol=True,algorithm='stevens') 
            sage: FsS == FsG
            True
        """
        p = self._get_prime(p)

        #Right now this is actually slower than the Stevens algorithm. Probably due to bad coding.
        
        #get a lift that is not a modular symbol
        MS = self.parent()
        gens = MS.source().gens()
        if new_base_ring == None:
            new_base_ring = MS.base_ring()
        MSnew = MS._lift_parent_space(p, M, new_base_ring)
        CMnew = MSnew.coefficient_module()
        D = {}
        gens = MS.source().gens()
        for j in range(len(gens)):
            D[gens[j]] = CMnew( self.values()[j]._moments.list() + [0] )
        Phi1bad = MSnew(D)
        
        #fix the lift by applying a hecke operator
        Phi1 = Phi1bad.hecke(p)
        Phi1 = Phi1/Phi1.Tq_eigenvalue(p)
        
        #if you don't want to compute with good accuracy, stop
        #otherwise, keep lifting
        if M<=2:
            return Phi1
        else:
            padic_prec=M + 1
            R = Qp(p,padic_prec)

            def green_lift_once(Phi1, self, r):
                newvalues = []
                for adist in Phi1.values():
                    newdist = [R(moment).lift_to_precision(moment.precision_absolute()+1) for moment in adist._moments] + [0]
                    for s in xrange(self.weight()+1):
                        newdist[s] = R(self.values()[Phi1.values().index(adist)].moment(s), r+1)
                    newvalues.append(newdist)
                D2 = {}
                for j in range(len(gens)):
                    D2[ gens[j]] = CMnew( newvalues[j] )
                Phi2 = MSnew(D2)
                Phi2 = Phi2.hecke(p)
                return Phi2 / self.Tq_eigenvalue(p)
 
            for r in range(self.weight() + 2, M):
                Phi1 = green_lift_once(Phi1,self,r)
        
            return Phi1
    
        D0 = {}
        for j in range(num_gens):
            D0[gens[j]] = CM1( [zero_moms[j]] + (M-1)*[0])

        #hecke and divide by eigenvalue
        Phi=MS1(D0)
        Phi=Phi.hecke(p)/ap
        #fix first moments, hecke and divide by eigenvalues
        for k in range(M-1):
            D1 = {}
            for j in range(num_gens):
                vals = Phi.values()[j]
                newvals=[vals.moment(n) for n in range(M)]
                newvals[0] = K(zero_moms[j])
                D1[gens[j]] = CM1(vals)
            Phi = MS1(D1)
            Phi=Phi.hecke(p)/ap
           
        return Phi
    
    def _lift_greenberg2(self, p, M, new_base_ring=None, check=False):
    #this is a slower version of the _lift_greenberg that tries not to
    #instantiate a bunch of parents. It turns out to be actually slower.
    #This code actually only works for weight 2 too. 
        MS = self.parent()
        gens=MS.source().gens()
        num_gens=len(gens)
        K=Qp(p,M)
        zero_moms=self.values()
        ap = self.Tq_eigenvalue(p)
    
        if new_base_ring == None:
            new_base_ring = MS.base_ring()
        MS1 = MS._lift_parent_space(p,M,new_base_ring)
        CM1=MS1.coefficient_module()
        D0 = {}
        for j in range(num_gens):
            D0[gens[j]] = CM1( [zero_moms[j]] + (M-1)*[0])

        #hecke and divide by eigenvalue
        Phi=MS1(D0)
        Phi=Phi.hecke(p)/ap
        
        #keep fixing first moments, hecke and divide by eigenvalues
        for k in range(M-1):
            D1 = {}
            for j in range(num_gens):
                vals = Phi.values()[j]
                newvals=[vals.moment(n) for n in range(M)]
                newvals[0] = K(zero_moms[j])
                D1[gens[j]] = CM1(vals)
            Phi = MS1(D1)
            Phi = Phi.hecke(p)/ap
        
        return Phi

    

    def _lift_to_OMS(self, p, M, new_base_ring, check):
        r"""
        Returns a (`p`-adic) overconvergent modular symbol with
        `M` moments which lifts self up to an Eisenstein error

        Here the Eisenstein error is a symbol whose system of Hecke
        eigenvalues equals `ell+1` for `T_ell` when `ell`
        does not divide `Np` and 1 for `U_q` when `q` divides `Np`.

        INPUT:

        - ``p`` -- prime

        - ``M`` -- integer equal to the number of moments

        - ``new_base_ring`` -- new base ring

        - ``check`` -- THIS IS CURRENTLY NOT USED IN THE CODE!

        OUTPUT:

        - An overconvergent modular symbol whose specialization
        equals self up to some Eisenstein error.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: E = EllipticCurve('11a')
            sage: f = ps_modsym_from_elliptic_curve(E)
            sage: f._lift_to_OMS(11,4,Qp(11,4),True)
            Modular symbol of level 11 with values in Space of 11-adic distributions with k=0 action and precision cap 4

        """
        D = {}
        manin = self.parent().source()
        MSS = self.parent()._lift_parent_space(p, M, new_base_ring)
        verbose("Naive lifting: newM=%s, new_base_ring=%s"%(M, MSS.base_ring()))
        half = ZZ(1) / ZZ(2)
        for g in manin.gens()[1:]:
            twotor = g in manin.reps_with_two_torsion()
            threetor = g in manin.reps_with_three_torsion()
            if twotor:
                # See [PS] section 4.1
                gam = manin.two_torsion_matrix(g)
                mu = self._map[g].lift(p, M, new_base_ring)
                D[g] = (mu - mu * gam) * half
            elif threetor:
                # See [PS] section 4.1
                gam = manin.three_torsion_matrix(g)
                mu = self._map[g].lift(p, M, new_base_ring)
                D[g] = (2 * mu - mu * gam - mu * (gam**2)) * half
            else:
                # no two or three torsion
                D[g] = self._map[g].lift(p, M, new_base_ring)

        t = self.parent().coefficient_module().lift(p, M, new_base_ring).zero_element()
        ## This loops adds up around the boundary of fundamental domain except the two vertical lines
        for g in manin.gens()[1:]:
            twotor = g in manin.reps_with_two_torsion()
            threetor = g in manin.reps_with_three_torsion()
            if twotor or threetor:
               t = t - D[g]
            else:
               t = t + D[g] * manin.gammas[g] - D[g]
        ## t now should be sum Phi(D_i) | (gamma_i - 1) - sum Phi(D'_i) - sum Phi(D''_i)
        ## (Here I'm using the opposite sign convention of [PS1] regarding D'_i and D''_i)

        D[manin.gen(0)] = -t.solve_diff_eqn()  ###### Check this!

        return MSS(D)

    def _find_aq(self, p, M, check):
        r"""
        Helper function for finding Hecke eigenvalue `aq` for a prime `q`
        not equal to `p`. This is called in the case when `alpha = 1 (mod p^M)`
        (with `alpha` a `U_p`-eigenvalue), which creates the need to use
        other Hecke eigenvalues (and `alpha`s), because of division by `(alpha - 1)`.

        INPUT:

        - ``p`` -- working prime

        - ``M`` -- precision

        - ``check`` -- checks that `self` is a `Tq` eigensymbol

        OUTPUT:

        Tuple `(q, aq, eisenloss)`, with

        - ``q`` -- a prime not equal to `p`

        - ``aq`` -- Hecke eigenvalue at `q`

        - ``eisenloss`` -- the `p`-adic valuation of `aq - q^(k+1) - 1`

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: E = EllipticCurve('11a')
            sage: f = ps_modsym_from_elliptic_curve(E)
            sage: f._find_aq(5,10,True)
            (2, -2, 1)
        """
        N = self.parent().level()
        q = ZZ(2)
        k = self.parent().weight()
        aq = self.Tq_eigenvalue(q, check=check)
        eisenloss = (aq - q**(k+1) - 1).valuation(p)
        while ((q == p) or (N % q == 0) or (eisenloss >= M)) and (q<50):
            q = next_prime(q)
            aq = self.Tq_eigenvalue(q, check=check)
            if q != p:
                eisenloss = (aq - q**(k+1) - 1).valuation(p)
            else:
                eisenloss = (aq - 1).valuation(p)
        if q >= 50:
            raise ValueError("The symbol appears to be eisenstein -- not implemented yet")
        return q, aq, eisenloss

    def _find_extraprec(self, p, M, alpha, check):
        q, aq, eisenloss = self._find_aq(p, M, check)
        newM = M + eisenloss
        # We also need to add precision to account for denominators appearing while solving the difference equation.
        eplog = (newM -1).exact_log(p)
        while eplog < (newM + eplog).exact_log(p):
            eplog = (newM + eplog).exact_log(p)
            verbose("M = %s, newM = %s, eplog=%s"%(M, newM, eplog), level=2)
        newM += eplog

        # We also need to add precision to account for denominators that might be present in self
        s = self.valuation(p)
        if s < 0:
            newM += -s
        return newM, eisenloss, q, aq

    def _lift_to_OMS_eigen(self, p, M, new_base_ring, ap, newM, eisenloss, q, aq, check):
        r"""
        Returns Hecke-eigensymbol OMS lifting self -- self must be a
        `p`-ordinary eigensymbol

        INPUT:

        - ``p`` -- prime

        - ``M`` -- integer equal to the number of moments

        - ``new_base_ring`` -- new base ring

        - ``ap`` -- Hecke eigenvalue at `p`

        - ``newM`` --

        - ``eisenloss`` --

        - ``q`` -- prime

        - ``aq`` -- Hecke eigenvalue at `q`

        - ``check`` --

        OUTPUT:

        - Hecke-eigenvalue OMS lifting self.

        EXAMPLES::



        """
        if new_base_ring(ap).valuation() > 0:
            raise ValueError("Lifting non-ordinary eigensymbols not implemented (issue #20)")

        verbose("computing naive lift: M=%s, newM=%s, new_base_ring=%s"%(M, newM, new_base_ring))
        Phi = self._lift_to_OMS(p, newM, new_base_ring, check)

        ## Scale by a large enough power of p to clear denominators from solving difference equation
#        s = newM.exact_log(p)+1
#        Phi = Phi * p**s

        ## Act by Hecke to ensure values are in D and not D^dag after sovling difference equation        
        verbose("Applying Hecke")
        apinv = ~ap
        Phi = apinv * Phi.hecke(p)

        verbose(Phi._show_malformed_dist("naive lift"), level=2)
        verbose(Phi._show_malformed_dist("after reduction"), level=2)

        ## Killing eisenstein part
        verbose("Killing eisenstein part with q = %s"%(q))
        k = self.parent().weight()
        Phi = ((q**(k+1) + 1) * Phi - Phi.hecke(q))
        verbose(Phi._show_malformed_dist("Eisenstein killed"), level=2)

        ## Iterating U_p
        verbose("Iterating U_p")
        Psi = apinv * Phi.hecke(p)
        attempts = 0
        while (Phi != Psi) and (attempts < 2*newM):
            verbose("%s attempt"%(attempts+1))
            Phi = Psi
            Psi = Phi.hecke(p) * apinv 
            attempts += 1
        if attempts >= 2*newM:
            raise RuntimeError("Precision problem in lifting -- applied U_p many times without success")
        Phi =  ~(q**(k+1) + 1 - aq) * Phi

        return Phi.reduce_precision(M)
        
    def p_stabilize_and_lift(self, p=None, M=None, alpha=None, ap=None, new_base_ring=None, \
                               ordinary=True, algorithm=None, eigensymbol=False, check=True):
        """
        `p`-stabilizes and lifts self

        INPUT:

        - ``p`` -- (default: None)

        - ``M`` -- (default: None)

        - ``alpha`` -- (default: None)

        - ``ap`` -- (default: None)

        - ``new_base_ring`` -- (default: None)

        - ``ordinary`` -- (default: True)

        - ``algorithm`` -- (default: None)

        - ``eigensymbol`` -- (default: False)

        - ``check`` -- (default: True)

        OUTPUT:

        `p`-stabilized and lifted version of self.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: E = EllipticCurve('11a')
            sage: f = ps_modsym_from_elliptic_curve(E)
            sage: g = f.p_stabilize_and_lift(3,10)
            sage: g.Tq_eigenvalue(5)
            1 + O(3^10)
            sage: g.Tq_eigenvalue(7)
            1 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + O(3^10)
            sage: g.Tq_eigenvalue(3)
            2 + 3^2 + 2*3^3 + 2*3^4 + 2*3^6 + 3^8 + 2*3^9 + O(3^10)
        """
        if check:
            p = self._get_prime(p, alpha)
        k = self.parent().weight()
        M = self._find_M(M)
        # alpha will be the eigenvalue of Up
        if alpha is None:
            alpha, new_base_ring, newM, eisenloss, q, aq = self._find_alpha(p, k, M, ap, new_base_ring, ordinary, check)
        else:
            if new_base_ring is None:
                new_base_ring = alpha.parent()
            newM, eisenloss, q, aq = self._find_extraprec(p, M, alpha, check)
            if hasattr(new_base_ring, 'precision_cap') and newM > new_base_ring.precision_cap():
                raise ValueError("Not enough precision in new base ring")

        # Now we can stabilize
        self = self.p_stabilize(p=p, alpha=alpha,ap=ap, M=newM, new_base_ring = new_base_ring, check=check)
        # And use the standard lifting function for eigensymbols
        return self._lift_to_OMS_eigen(p=p, M=M, new_base_ring=new_base_ring, ap=alpha, newM=newM, eisenloss=eisenloss, q=q, aq=aq, check=check)

class PSModularSymbolElement_dist(PSModularSymbolElement):

    def _show_malformed_dist(self, location_str):
        malformed = []
        gens = self.parent().source().gens()
        for j, g in enumerate(gens):
            val = self._map[g]
            if val._is_malformed():
                malformed.append((j, val))
        return location_str + ": (%s/%s malformed)%s"%(len(malformed), len(gens), ", %s -- %s"%(malformed[0][0], str(malformed[0][1])) if len(malformed) > 0 else "")

    def reduce_precision(self, M):
        r"""
        Only holds on to `M` moments of each value of self
        """
        return self.__class__(self._map.reduce_precision(M), self.parent(), construct=True)

    def precision_absolute(self):
        r"""
        Returns the number of moments of each value of self
        """
        return min([a.precision_absolute() for a in self._map])

    def specialize(self, new_base_ring=None):
        r"""
        Returns the underlying classical symbol of weight `k` -- i.e.,
        applies the canonical map `D_k --> Sym^k` to all values of
        self.

        EXAMPLES::

            sage: D = Distributions(0, 5, 10);  M = PSModularSymbols(Gamma0(5), coefficients=D); M
            Space of overconvergent modular symbols for Congruence Subgroup Gamma0(5) with sign 0 and values in Space of 5-adic distributions with k=0 action and precision cap 10
            sage: f = M(1)
            sage: f.specialize()
            Modular symbol of level 5 with values in Sym^0 Z_5^2
            sage: f.specialize().values()
            [1 + O(5^10), 1 + O(5^10), 1 + O(5^10)]
            sage: f.values()
            [1, 1, 1]
            sage: f.specialize().parent()
            Space of modular symbols for Congruence Subgroup Gamma0(5) with sign 0 and values in Sym^0 Z_5^2
            sage: f.specialize().parent().coefficient_module()
            Sym^0 Z_5^2
            sage: f.specialize().parent().coefficient_module().is_symk()
            True

            sage: f.specialize(QQ)
            Modular symbol of level 5 with values in Sym^0 Q^2
            sage: f.specialize(QQ).values()
            [1, 1, 1]
            sage: f.specialize(QQ).parent().coefficient_module()
            Sym^0 Q^2
        """
        if new_base_ring is None:
            new_base_ring = self.base_ring()
        return self.__class__(self._map.specialize(new_base_ring),
                              self.parent()._specialize_parent_space(new_base_ring), construct=True)

    def padic_lseries(self,*args, **kwds):
        r"""
        Return the p-adic L-series of this modular symbol.

        EXAMPLE::

            sage: f = Newform("37a")
            sage: f.PS_modular_symbol().lift(37, M=6, algorithm="stevens").padic_lseries()
            37-adic L-series of Modular symbol of level 37 with values in Space of 37-adic distributions with k=0 action and precision cap 6
        """
        return pAdicLseries(self, *args, **kwds)
