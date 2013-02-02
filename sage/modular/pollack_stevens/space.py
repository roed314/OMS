r"""
Pollack-Stevens Modular Symbols Spaces

This module contains a class for spaces of modular symbols that use Glenn
Stevens' conventions.

There are two main differences between the modular symbols in this directory
and the ones in :mod:`sage.modular.modsym`:

- There is a shift in the weight: weight `k=0` here corresponds to weight `k=2`
  there.

- There is a duality: these modular symbols are functions from
  `Div^0(P^1(\QQ))` (cohomological objects), the others are formal linear
  combinations of `Div^0(P^1(\QQ))` (homological objects).
"""
#*****************************************************************************
#       Copyright (C) 2012 Robert Pollack <rpollack@math.bu.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.modules.module import Module
from sage.structure.factory import UniqueFactory
from distributions import Distributions, Symk
from sage.modular.dirichlet import DirichletCharacter
from sage.modular.arithgroup.all import Gamma0
from sage.modular.arithgroup.arithgroup_element import ArithmeticSubgroupElement
from sage.rings.arith import binomial
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from modsym import PSModularSymbolElement_symk, PSModularSymbolElement_dist, PSModSymAction
from fund_domain import ManinRelations, M2ZSpace
from sage.rings.padics.precision_error import PrecisionError
from sage.rings.infinity import infinity as oo

class PSModularSymbols_factory(UniqueFactory):
    r"""
    Create a space of Pollack-Stevens modular symbols.

    INPUT:

    - ``group`` -- integer or congruence subgroup
    - ``weight`` -- integer `\ge 2`, or None
    - ``sign`` -- integer; -1, 0, 1
    - ``base_ring`` -- ring or None
    - `p` -- prime or None
    - ``prec_cap`` -- positive integer or None
    - ``coefficients`` -- the coefficient module (a special type of module,
      typically distributions), or None.

    If an explicit coefficient module is given, then the arguments ``weight``,
    ``base_ring``, ``prec_cap`` and possibly also ``p`` are redundant and are
    ignored. They are only relevant if ``coefficients`` is None, in which case
    the coefficient module is inferred from the other data.

    EXAMPLES::

        sage: PSModularSymbols(Gamma0(7), weight=2, prec_cap = None)
        Space of modular symbols for Congruence Subgroup Gamma0(7) with sign 0 and values in Sym^0 Q^2

    An example with an explict coefficient module::

        sage: D = Distributions(3, 7, prec_cap=10)
        sage: PSModularSymbols(Gamma0(7), coefficients=D)
        Space of overconvergent modular symbols for Congruence Subgroup Gamma0(7) with sign 0 and values in Space of 7-adic distributions with k=3 action and precision cap 10 

    TESTS::

        sage: TestSuite(PSModularSymbols).run()
    """
    def create_key(self, group, weight=None, sign=0, base_ring=None, p=None, prec_cap=None, coefficients=None):
        r"""
        Sanitize input.

        EXAMPLES::

            sage: D = Distributions(3, 7, prec_cap=10)
            sage: M = PSModularSymbols(Gamma0(7), coefficients=D) # indirect doctest
            sage: TestSuite(PSModularSymbols).run()
        """
        if isinstance(group, (int, Integer)):
            group = Gamma0(group)
        if coefficients is None:
            if isinstance(group, DirichletCharacter):
                character = group.minimize_base_ring()
                group = Gamma0(character.modulus())
                character = (character, None)
            else:
                character = None
            if weight is None: raise ValueError("you must specify a weight or coefficient module")
            k = weight - 2
            if prec_cap is None:
                coefficients = Symk(k, base_ring, character)
            else:
                coefficients = Distributions(k, p, prec_cap, base_ring, character)
        else:
            # TODO: require other stuff to be None
            pass
        return (group, coefficients, sign)

    def create_object(self, version, key):
        r"""
        EXAMPLES::

            sage: D = Distributions(5, 7, 15); M = PSModularSymbols(Gamma0(7), coefficients=D) # indirect doctest
        """
        return PSModularSymbolSpace(*key)

PSModularSymbols = PSModularSymbols_factory('PSModularSymbols')

class PSModularSymbolSpace(Module):
    r"""
    A class for spaces of modular symbols that use Glenn Stevens'
    conventions. This class should not be instantiated directly by the user:
    this is handled by the factory object ``PSModularSymbols``.
    """

    def __init__(self, group, coefficients, sign=0):
        r"""
        INPUT:

        - ``group`` -- congruence subgroup
        - ``coefficients`` -- a coefficient module
        - ``sign`` -- (default: 0); 0, -1, 1

        EXAMPLES::

            sage: D = Distributions(2, 11);  M = PSModularSymbols(Gamma0(2), coefficients=D)
            sage: type(M)
            <class 'sage.modular.pollack_stevens.space.PSModularSymbolSpace_with_category'>
            sage: M.sign()
            0
            sage: M = PSModularSymbols(Gamma0(2), coefficients=D, sign=-1); M.sign()
            -1
            sage: M = PSModularSymbols(Gamma0(2), coefficients=D, sign=1); M.sign()
            1
            sage: TestSuite(M).run()
        """
        Module.__init__(self, coefficients.base_ring())
        if sign not in [0,-1,1]:
            # sign must be be 0, -1 or 1
            raise ValueError, "sign must be 0, -1, or 1"
        self._group = group
        self._coefficients = coefficients
        if coefficients.is_symk():
            self.Element = PSModularSymbolElement_symk
        else:
            self.Element = PSModularSymbolElement_dist
        self._sign = sign
        # should distingish between Gamma0 and Gamma1...
        self._source = ManinRelations(group.level())
        # We have to include the first action so that scaling by Z doesn't try to pass through matrices
        actions = [PSModSymAction(ZZ, self), PSModSymAction(M2ZSpace, self)]
        self._populate_coercion_lists_(action_list=actions)

    def _repr_(self):
        r"""
        Returns print representation.

        EXAMPLES::

            sage: D = Distributions(2, 11);  M = PSModularSymbols(Gamma0(2), coefficients=D); M._repr_()
            'Space of overconvergent modular symbols for Congruence Subgroup Gamma0(2) with sign 0 and values in Space of 11-adic distributions with k=2 action and precision cap 20'
        """
        if self.coefficient_module().is_symk():
            s = "Space of modular symbols for "
        else:
            s = "Space of overconvergent modular symbols for "
        s += "%s with sign %s and values in %s"%(self.group(), self.sign(), self.coefficient_module())
        return s

    def source(self):
        r"""
        Return object that represents the domain of the modular
        symbol elements.

        OUTPUT:

        - object of type fund_domain.PSModularSymbolsDomain

        EXAMPLES::

            sage: D = Distributions(2, 11);  M = PSModularSymbols(Gamma0(2), coefficients=D)
            sage: M.source()
            Manin Relations of level 2
        """
        return self._source

    def coefficient_module(self):
        r"""
        Returns the coefficient module of self.

        OUTPUT:

        - module

        EXAMPLES::

            sage: D = Distributions(2, 11);  M = PSModularSymbols(Gamma0(2), coefficients=D)
            sage: M.coefficient_module()
            Space of 11-adic distributions with k=2 action and precision cap 20
            sage: M.coefficient_module() == D
            True
            sage: M.coefficient_module() is D
            True
        """
        return self._coefficients

    def group(self):
        r"""
        Returns the congruence subgroup of self.

        EXAMPLES::

            sage: D = Distributions(2, 5)
            sage: G = Gamma0(23)
            sage: M = PSModularSymbols(G, coefficients=D)
            sage: M.group()
            Congruence Subgroup Gamma0(23)
            sage: D = Symk(4)
            sage: G = Gamma1(11)
            sage: M = PSModularSymbols(G, coefficients=D)
            sage: M.group()
            Congruence Subgroup Gamma1(11)
        """
        return self._group

    def sign(self):
        r"""
        Returns the sign of self.

        EXAMPLES::

            sage: D = Distributions(3, 17)
            sage: M = PSModularSymbols(Gamma(5), coefficients=D)
            sage: M.sign()
            0
            sage: D = Symk(4)
            sage: M = PSModularSymbols(Gamma1(8), coefficients=D, sign=-1)
            sage: M.sign()
            -1
        """
        return self._sign

    def ngens(self):
        r"""
        Returns the number of generators defining self.

        EXAMPLES::

            sage: D = Distributions(4, 29)
            sage: M = PSModularSymbols(Gamma1(12), coefficients=D)
            sage: M.ngens()
            5
            sage: D = Symk(2)
            sage: M = PSModularSymbols(Gamma0(2), coefficients=D)
            sage: M.ngens()
            2
        """
        return len(self._source.indices())

    def ncoset_reps(self):
        r"""
        Returns the number of coset representatives defining the
        full_data of self.

        OUTPUT:

        - The number of coset representatives stored in the manin
          relations. (Just the size of P^1(Z/NZ))

        EXAMPLES::

            sage: D = Symk(2)
            sage: M = PSModularSymbols(Gamma0(2), coefficients=D)
            sage: M.ncoset_reps()
            3
        """
        return len(self._source.reps())

    def level(self):
        r"""
        Returns the level `N`, where self is of level `Gamma_0(N)`.

        OUTPUT:

        - The level `N`, where self is of level `Gamma_0(N)`

        EXAMPLES::

            sage: D = Distributions(7, 11)
            sage: M = PSModularSymbols(Gamma1(14), coefficients=D)
            sage: M.level()
            14
        """
        return self._source.level()

    def _grab_relations(self):
        r"""
        This is used internally as part of a consistency check.

        EXAMPLES::

            sage: D = Distributions(4, 3)
            sage: M = PSModularSymbols(Gamma1(13), coefficients=D)
            sage: M._grab_relations()
            [[(1, [1 0]
            [0 1], 0)], [(-1, [-1 -1]
            [ 0 -1], 0)], [(1, [1 0]
            [0 1], 2)], [(1, [1 0]
                [0 1], 3)], [(1, [1 0]
                    [0 1], 4)], [(1, [1 0]
                        [0 1], 5)]]
        """
        v = []
        for r in range(len(self._source.gens())):
            for j in range(len(self._source.reps())):
                R = self._source.relations(j)
                if len(R) == 1 and R[0][2] == self._source.indices(r):
                    if R[0][0] != -1 or R[0][1] != M2ZSpace.one():
                        v = v + [R]
        return v



    def precision_cap(self):
        r"""
        Returns the number of moments of each value of self.

        OUTPUT:

        - Integer

        EXAMPLES::

            sage: D = Distributions(2, 5)
            sage: M = PSModularSymbols(Gamma1(13), coefficients=D)
            sage: M.precision_cap()
            20
            sage: D = Distributions(3, 7, prec_cap=10)
            sage: M = PSModularSymbols(Gamma0(7), coefficients=D)
            sage: M.precision_cap()
            10
        """
        return self.coefficient_module()._prec_cap

    def weight(self):
        r"""
        Returns the weight of self.

        We emphasize that in the Pollack-Stevens notation, this is the
        usual weight minus 2, so a classical weight 2 modular form corresponds
        to a modular symbol of "weight 0".

        OUTPUT:

        - Integer

        EXAMPLES::

            sage: D = Symk(5)
            sage: M = PSModularSymbols(Gamma1(7), coefficients=D)
            sage: M.weight()
            5
        """
        return self.coefficient_module()._k

    def prime(self):
        r"""
        Returns the prime of self.

        OUTPUT:

        - Integer

        EXAMPLES:
            sage: D = Distributions(2, 11)
            sage: M = PSModularSymbols(Gamma(2), coefficients=D)
            sage: M.prime()
            11
        """
        return self.coefficient_module()._p

    def _p_stabilize_parent_space(self, p, new_base_ring):
        r"""
        Returns the space of Pollack-Stevens modular symbols of level
        p * N, with changed base ring.  This is used internally when
        constructing the p-stabilization of a modular symbol.

        INPUT:

        - ``p`` -- prime number
        - ``new_base_ring`` -- the base ring of the result

        OUTPUT:

        - The space of modular symbols of level p * N, where N is the level of
          self, with precision M.

        EXAMPLES::

            sage: D = Distributions(2, 7); M = PSModularSymbols(Gamma(13), coefficients=D)
            sage: M._p_stabilize_parent_space(7, M.base_ring())
            Space of overconvergent modular symbols for Congruence Subgroup
            Gamma(91) with sign 0 and values in Space of 7-adic distributions
            with k=2 action and precision cap 20

            sage: D = Distributions(4, 17); M = PSModularSymbols(Gamma1(3), coefficients=D)
            sage: M._p_stabilize_parent_space(17, Qp(17))
            Space of overconvergent modular symbols for Congruence
            Subgroup Gamma1(51) with sign 0 and values in Space of
            17-adic distributions with k=4 action and precision cap 20
        """

        N = self.level()
        if N % p == 0:
            raise ValueError("the level isn't prime to p")
        from sage.modular.arithgroup.all import Gamma, is_Gamma, Gamma0, is_Gamma0, Gamma1, is_Gamma1
        G = self.group()
        if is_Gamma0(G):
            G = Gamma0(N*p)
        elif is_Gamma1(G):
            G = Gamma1(N*p)
        elif is_Gamma(G):
            G = Gamma(N*p)
        else:
            raise NotImplementedError
        return PSModularSymbols(G, coefficients=self.coefficient_module().change_ring(new_base_ring), sign=self.sign())

    def _specialize_parent_space(self, new_base_ring):
        r"""
        Internal function that is used by the specialize method on
        elements.  It returns a space with same parameters as this
        one, but over new_base_ring.

        EXAMPLES::

            sage: D = Distributions(7, 5);  M = PSModularSymbols(Gamma0(2), coefficients=D); M
            Space of overconvergent modular symbols for Congruence Subgroup Gamma0(2) with sign 0 and values in Space of 5-adic distributions with k=7 action and precision cap 20
            sage: M._specialize_parent_space(QQ)
            Space of modular symbols for Congruence Subgroup Gamma0(2) with sign 0 and values in Sym^7 Q^2
            sage: M.base_ring()
            5-adic Ring with capped absolute precision 20
            sage: M._specialize_parent_space(QQ).base_ring()
            Rational Field
        """
        return PSModularSymbols(self.group(), coefficients=self.coefficient_module().specialize(new_base_ring), sign=self.sign())

    def _lift_parent_space(self, p, M, new_base_ring):
        r"""
        Used internally when lifting modular symbols.

        INPUT:

        - `p` -- prime
        - `M` -- precision cap
        - ``new_base_ring`` -- ring

        EXAMPLES::

            sage: D = Distributions(4, 17, 2); M = PSModularSymbols(Gamma1(3), coefficients=D)
            sage: D.is_symk()
            False
            sage: M._lift_parent_space(17, 10, Qp(17))
            Traceback (most recent call last):
            ...
            TypeError: Coefficient module must be a Symk
            sage: PSModularSymbols(Gamma1(3), weight=3)._lift_parent_space(17,10,Qp(17))
            Space of overconvergent modular symbols for Congruence Subgroup Gamma1(3) with sign 0 and values in Space of 17-adic distributions with k=1 action and precision cap 10
        """
        if self.coefficient_module().is_symk():
            return PSModularSymbols(self.group(), coefficients=self.coefficient_module().lift(p, M, new_base_ring), sign=self.sign())
        else:
            raise TypeError("Coefficient module must be a Symk")

    def change_ring(self, new_base_ring):
        r"""
        Changes base ring of self to new_base_ring

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Symk
            sage: D = Symk(4)
            sage: M = PSModularSymbols(Gamma(6), coefficients=D); M
            Space of modular symbols for Congruence Subgroup Gamma(6) with sign 0 and values in Sym^4 Q^2
            sage: M.change_ring(Qp(5,8))
            Space of modular symbols for Congruence Subgroup Gamma(6) with sign 0 and values in Sym^4 Q_5^2
        """
        return PSModularSymbols(self.group(), coefficients=self.coefficient_module().change_ring(new_base_ring), sign=self.sign())

    def _an_element_(self):
        r"""
        Returns a "typical" element of self; in this case the constant map sending every element
        to an element of the coefficient module.

        EXAMPLES::

            sage: D = Symk(4)
            sage: M = PSModularSymbols(Gamma(6), coefficients=D)
            sage: x = M.an_element(); x       # indirect doctest
            Modular symbol with values in Sym^4 Q^2
            sage: x.values()
            [(0, 1, 2, 3, 4), (0, 1, 2, 3, 4), (0, 1, 2, 3, 4)]
            sage: D = Symk(2, Qp(11)); M = PSModularSymbols(Gamma0(2), coefficients=D)
            sage: x = M.an_element(); x.values()
            [(0, 1 + O(11^20), 2 + O(11^20)), (0, 1 + O(11^20), 2 + O(11^20))]
            sage: x in M
            True
        """
        return self(self.coefficient_module().an_element())


def cusps_from_mat(g):
    r"""
    Returns the cusps associated to an element of a congruence
    subgroup.

    INPUT:

    - `g` -- an element of a congruence subgroup or a matrix

    OUTPUT:

    - The cusps associated to `g`.

    EXAMPLES::

        sage: from sage.modular.pollack_stevens.space import cusps_from_mat
        sage: g = SL2Z.one()
        sage: cusps_from_mat(g)
        (+Infinity, 0)

    You can also just give the matrix of g::

        sage: type(g)
        <class 'sage.modular.arithgroup.arithgroup_element.ArithmeticSubgroupElement'>
        sage: cusps_from_mat(g.matrix())
        (+Infinity, 0)


    Another example::

        sage: from sage.modular.pollack_stevens.space import cusps_from_mat
        sage: g = GammaH(3, [2]).generators()[1].matrix(); g
        [-1  1]
        [-3  2]
        sage: cusps_from_mat(g)
        (1/3, 1/2)
    """
    if isinstance(g, ArithmeticSubgroupElement):
        g = g.matrix()
    a, b, c, d = g.list()
    if c: ac = a/c
    else: ac = oo
    if d: bd = b/d
    else: bd = oo
    return ac, bd

def ps_modsym_from_elliptic_curve(E):
    r"""
    Returns the PS modular symbol associated to an elliptic curve
    defined over the rationals.

    INPUT:

    - ``E`` -- an elliptic curve defined over the rationals

    OUTPUT:

    - the Pollack-Stevens modular symbol associated to ``E``

    EXAMPLES::

        sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
        sage: E = EllipticCurve('113a1')
        sage: symb = ps_modsym_from_elliptic_curve(E)
        sage: symb
        Modular symbol with values in Sym^0 Q^2
        sage: symb.values()
        [-1/2, 3/2, -2, 1/2, 0, 1, 2, -3/2, 0, -3/2, 0, -1/2, 0, 1, -2, 1/2, 0,
        0, 2, 0, 0]

        sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
        sage: E = EllipticCurve([0,1])
        sage: symb = ps_modsym_from_elliptic_curve(E)
        sage: symb.values()
        [-1/6, 7/12, 1, 1/6, -5/12, 1/3, -7/12, -1, -1/6, 5/12, 1/4, -1/6, -5/12]
    """
    if not (E.base_ring() is QQ):
        raise ValueError("The elliptic curve must be defined over the rationals.")
    N = E.conductor()
    V = PSModularSymbols(Gamma0(N), 2)
    D = V.coefficient_module()
    manin = V.source()
    plus_sym = E.modular_symbol(sign = 1)
    minus_sym = E.modular_symbol(sign = -1)
    val = {}
    for g in manin.gens():
        ac, bd = cusps_from_mat(g)
        val[g] = D([plus_sym(ac) + minus_sym(ac) - plus_sym(bd) - minus_sym(bd)])
    return V(val)

def ps_modsym_from_simple_modsym_space(A):
    r"""
    Returns some choice -- only well defined up a nonzero scalar (!)
    -- of a Pollack-Stevens modular symbol that corresponds to A.

    INPUT:

    - `A` -- nonzero simple Hecke equivariant new space of modular
      symbols, which need not be cuspidal.

    OUTPUT:

    - A choice of corresponding Pollack-Stevens modular symbols; when
      dim(A)>1, we make an arbitrary choice of defining polynomial for
      the codomain field.

    EXAMPLES::

        sage: from sage.modular.pollack_stevens.space import ps_modsym_from_simple_modsym_space

    The level 11 example::

        sage: A = ModularSymbols(11, sign=1, weight=2).decomposition()[0]
        sage: A.is_cuspidal()
        True
        sage: f = ps_modsym_from_simple_modsym_space(A); f
        Modular symbol with values in Sym^0 Q^2
        sage: f.values()
        [1, -5/2, -5/2]
        sage: f.weight()         # this is A.weight()-2  !!!!!!
        0

    And the -1 sign for the level 11 example::

        sage: A = ModularSymbols(11, sign=-1, weight=2).decomposition()[0]
        sage: f = ps_modsym_from_simple_modsym_space(A); f.values()
        [0, 1, -1]

    A does not have to be cuspidal; it can be Eisenstein::

        sage: A = ModularSymbols(11, sign=1, weight=2).decomposition()[1]
        sage: A.is_cuspidal()
        False
        sage: f = ps_modsym_from_simple_modsym_space(A); f
        Modular symbol with values in Sym^0 Q^2
        sage: f.values()
        [1, 0, 0]

    We create the simplest weight 2 example in which A has dimension
    bigger than 1::

        sage: A = ModularSymbols(23, sign=1, weight=2).decomposition()[0]
        sage: f = ps_modsym_from_simple_modsym_space(A); f.values()
        [1, 0, 0, 0, 0]
        sage: A = ModularSymbols(23, sign=-1, weight=2).decomposition()[0]
        sage: f = ps_modsym_from_simple_modsym_space(A); f.values()
        [0, 1, -alpha, alpha, -1]
        sage: f.base_ring()
        Number Field in alpha with defining polynomial x^2 + x - 1

    We create the +1 modular symbol attached to the weight 12 modular form Delta::

        sage: A = ModularSymbols(1, sign=+1, weight=12).decomposition()[0]
        sage: f = ps_modsym_from_simple_modsym_space(A); f
        Modular symbol with values in Sym^10 Q^2
        sage: f.values()
        [(-1620/691, 0, 1, 0, -9/14, 0, 9/14, 0, -1, 0, 1620/691), (1620/691, 1620/691, 929/691, -453/691, -29145/9674, -42965/9674, -2526/691, -453/691, 1620/691, 1620/691, 0), (0, -1620/691, -1620/691, 453/691, 2526/691, 42965/9674, 29145/9674, 453/691, -929/691, -1620/691, -1620/691)]

    And, the -1 modular symbol attached to Delta::

        sage: A = ModularSymbols(1, sign=-1, weight=12).decomposition()[0]
        sage: f = ps_modsym_from_simple_modsym_space(A); f
        Modular symbol with values in Sym^10 Q^2
        sage: f.values()
        [(0, 1, 0, -25/48, 0, 5/12, 0, -25/48, 0, 1, 0), (0, -1, -2, -119/48, -23/12, -5/24, 23/12, 3, 2, 0, 0), (0, 0, 2, 3, 23/12, -5/24, -23/12, -119/48, -2, -1, 0)]

    A consistency check with ps_modsym_from_simple_modsym_space::

        sage: from sage.modular.pollack_stevens.space import (ps_modsym_from_elliptic_curve, ps_modsym_from_simple_modsym_space)
        sage: E = EllipticCurve('11a')
        sage: f_E = ps_modsym_from_elliptic_curve(E); f_E.values()
        [-1/5, 3/2, -1/2]
        sage: A = ModularSymbols(11, sign=1, weight=2).decomposition()[0]
        sage: f_plus = ps_modsym_from_simple_modsym_space(A); f_plus.values()
        [1, -5/2, -5/2]
        sage: A = ModularSymbols(11, sign=-1, weight=2).decomposition()[0]
        sage: f_minus = ps_modsym_from_simple_modsym_space(A); f_minus.values()
        [0, 1, -1]

    We find that a linear combination of the plus and minus parts
    equals the Pollack-Stevens symbol attached to E. This illustrates
    how ps_modsym_from_simple_modsym_space is only well-defined up to
    a nonzero scalar.

        sage: (-1/5)*vector(QQ, f_plus.values()) + vector(QQ, f_minus.values())
        (-1/5, 3/2, -1/2)
        sage: vector(QQ, f_E.values())
        (-1/5, 3/2, -1/2)

    The next few examples all illustrate the ways in which exceptions are raised
    if A does not satisfy various constraints.

    First, A must be new::

        sage: A = ModularSymbols(33,sign=1).cuspidal_subspace().old_subspace()
        sage: ps_modsym_from_simple_modsym_space(A)
        Traceback (most recent call last):
        ...
        ValueError: A must be new

    A must be simple::

        sage: A = ModularSymbols(43,sign=1).cuspidal_subspace()
        sage: ps_modsym_from_simple_modsym_space(A)
        Traceback (most recent call last):
        ...
        ValueError: A must be simple

    A must have sign -1 or +1 in order to be simple::

        sage: A = ModularSymbols(11).cuspidal_subspace()
        sage: ps_modsym_from_simple_modsym_space(A)
        Traceback (most recent call last):
        ...
        ValueError: A must have sign +1 or -1 (otherwise it is not simple)

    The dimension must be positive::

        sage: A = ModularSymbols(10).cuspidal_subspace(); A
        Modular Symbols subspace of dimension 0 of Modular Symbols space of dimension 3 for Gamma_0(10) of weight 2 with sign 0 over Rational Field
        sage: ps_modsym_from_simple_modsym_space(A)
        Traceback (most recent call last):
        ...
        ValueError: A must positive dimension
    """
    if A.dimension() == 0:
        raise ValueError, "A must positive dimension"

    if A.sign() == 0:
        raise ValueError, "A must have sign +1 or -1 (otherwise it is not simple)"

    if not A.is_new():
        raise ValueError, "A must be new"

    if not A.is_simple():
        raise ValueError, "A must be simple"

    M = A.ambient_module()
    w = A.dual_eigenvector()
    K = w.base_ring()
    V = PSModularSymbols(A.group(), A.weight(), base_ring=K, sign=A.sign())
    D = V.coefficient_module()
    N = V.level()
    k = V.weight() # = A.weight() - 2
    manin = V.source()
    val = {}
    for g in manin.gens():
        ac, bd = cusps_from_mat(g)
        v = []
        for j in range(k+1):
            # TODO: The following might be backward: it should be the coefficient of X^j Y^(k-j)
            v.append(w.dot_product(M.modular_symbol([j, ac, bd]).element()) * (-1)**(k-j))
        val[g] = D(v)
    return V(val)
