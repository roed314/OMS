from sage.modular.btquotients.btquotient import *
from sage.modular.btquotients.ocmodule import *
#########################################################################
#       Copyright (C) 2011 Cameron Franc and Marc Masdeu
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################
from collections import namedtuple
from sage.structure.element import Element, ModuleElement
from sage.structure.parent import Parent
from sage.modules.module import Module
from sage.rings.all import Integer
from sage.structure.element import Element
from sage.matrix.constructor import Matrix, zero_matrix
from sage.rings.all import Qp
from sage.rings.all import RationalField
from sage.rings.number_field.all import NumberField
from copy import copy
from sage.quadratic_forms.quadratic_form import QuadraticForm
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.modular.hecke.all import (AmbientHeckeModule, HeckeSubmodule, HeckeModuleElement)
from sage.rings.infinity import Infinity
import sage.rings.arith as arith
import sage.modular.hecke.hecke_operator
from sage.modular.pollack_stevens.distributions import Distributions, Symk
from sage.misc.misc import verbose, cputime

# act_S0p = WeightKAction(self, character, tuplegen, act_on_left, padic = True)

# Need this to be pickleable
class _btquot_tuplegen(UniqueRepresentation):
    """
    Callable object that turns matrices into 4-tuples.

    EXAMPLES::

        sage: A = sage.modular.pollack_stevens.distributions._default_tuplegen(); A
        <sage.modular.pollack_stevens.distributions._default_tuplegen object at 0x...>
        sage: TestSuite(A).run()
    """
    def __call__(self, g):
        """
        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
            sage: T = sage.modular.pollack_stevens.distributions._default_tuplegen()
            sage: T(matrix(ZZ,2,[1..4]))
            (4, 2, 3, 1)
        """
        return g[1,1], g[0,1], g[1,0], g[0,0]

def act_left(g,v):
    group = MatrixSpace(v.base_ring(),2,2)
    return group(g) * v

class HarmonicCocycleElement(HeckeModuleElement):
    r"""
    Gamma-invariant harmonic cocycles on the Bruhat-Tits
    tree. Gamma-invariance is necessary so that the cocycle can be
    stored in terms of a finite amount of data.

    More precisely, given a BTQuotient T, harmonic cocycles are stored as
    a list of values in some coefficient module (e.g. for weight 2 forms
    can take Cp) indexed by edges of a fundamental domain for T in the
    Bruhat-Tits tree. Evaluate the cocycle at other edges using Gamma
    invariance (although the values may not be equal over an orbit of
    edges as the coefficient module action may be nontrivial).

    INPUT:

    - ``vec`` - (default: None)

    - ``from_values`` -  (default: False)

    EXAMPLES:

    Harmonic cocycles form a vector space, so they can be added::

        sage: X = BTQuotient(5,23)
        sage: H = HarmonicCocycles(X,2,prec=10)
        sage: v1 = H.basis()[0]; v2 = H.basis()[1] # indirect doctest
        sage: v3 = v1+v2
        sage: v1 == v3-v2
        True

    and rescaled::

        sage: v4 = 2*v1
        sage: v1 == v4 - v1
        True

    AUTHORS:

    - Cameron Franc (2012-02-20)
    - Marc Masdeu
    """
    def __init__(self,_parent,vec = None,from_values = False):
        """
        Create a harmonic cocycle element.

        EXAMPLES::

            sage: X = BTQuotient(31,7)
            sage: H = HarmonicCocycles(X,2,prec=10)
            sage: v = H.basis()[0] # indirect doctest
            sage: TestSuite(v).run()
        """
        HeckeModuleElement.__init__(self,_parent,None)
        self._parent = _parent

        if from_values:
            self._R = _parent._U.base_ring()
            self._wt = _parent._k
            self._nE = len(_parent._E)
            self._F = copy(vec)
            return
        if isinstance(vec, self.__class__):
            #The first argument is just a modular form that we copy as is
            self._R = vec._R
            self._nE = vec._nE
            self._wt = vec._wt
            self._F = [_parent._U(o) for o in vec._F]
            return
        # When vec contains coordinates for the basis
        self._R  =  _parent._R
        try:
            v = [self._R(x) for x in vec.list()]
        except AttributeError:
            v = [self._R(vec) for ii in range(_parent.dimension())]
        self._wt = _parent._k
        self._nE = len(_parent._E)
        vmat = Matrix(self._R,1,_parent.dimension(),v)
        tmp = (vmat*_parent.ambient_module().basis_matrix()).row(0)
        # self._F = [_parent._U(Matrix(self._R,self._wt-1,1,tmp[e*(self._wt-1):(e+1)*(self._wt-1)]),check = False) for e in range(self._nE)]
        self._F = [_parent._U(tmp[e*(self._wt-1):(e+1)*(self._wt-1)]) for e in range(self._nE)]
        return

    def _add_(self,g):
        r"""
        Add two cocycles componentwise.

        EXAMPLES::

            sage: X = BTQuotient(5,23)
            sage: H = HarmonicCocycles(X,2,prec=10)
            sage: v1 = H.basis()[0]; v2 = H.basis()[1]
            sage: v3 = v1+v2 # indirect doctest
            sage: v1 == v3-v2
            True
        """
        #Should ensure that self and g are modular forms of the same weight and on the same curve
        C = self.__class__
        return C(self.parent(),self.element()+g.element())

    def _sub_(self,g):
        r"""
        Computes the difference of two cocycles.

        EXAMPLES::

            sage: X = BTQuotient(5,23)
            sage: H = HarmonicCocycles(X,2,prec=10)
            sage: v1 = H.basis()[0]; v2 = H.basis()[1]
            sage: v3 = v1-v2 # indirect doctest
            sage: v1 == v3+v2
            True
        """
        #Should ensure that self and g are modular forms of the same weight and on the same curve
        C = self.__class__
        return C(self.parent(),self.element()-g.element())

    def _rmul_(self,a):
        r"""
        Multiplies a cocycle by a scalar.

        EXAMPLES::

            sage: X = BTQuotient(5,23)
            sage: H = HarmonicCocycles(X,2,prec=10)
            sage: v1 = H.basis()[0]
            sage: v2 = 2*v1 # indirect doctest
            sage: v1 == v2-v1
            True
        """
        #Should ensure that 'a' is a scalar
        C = self.__class__
        return C(self.parent(),a*self.element())


    def _repr_(self):
        r"""
        Retuns a string describing the cocycle.

        EXAMPLES::

            sage: X = BTQuotient(5,23)
            sage: H = HarmonicCocycles(X,2,prec=10)
            sage: v1 = H.basis()[0]
            'Element of Space of harmonic cocycles of weight 2 on Quotient of the Bruhat Tits tree of GL_2(QQ_5) with discriminant 23 and level 1'
        """
        tmp = 'Element of '+str(self.parent())
        return tmp

    def _eq_(self,other):
        r"""
        Test for equality with another cocycle. Two cocycles are
        equal if they take the same values on all the edges in a
        fundamental domain.

        EXAMPLES::

            sage: X = BTQuotient(5,23)
            sage: H = HarmonicCocycles(X,2,prec=10)
            sage: v1 = H.basis()[0]; v2 = H.basis()[1]
            sage: v1 == v2 # indirect doctest
            False
            sage: v1 == v1
            True
        """
        #This is needs work. For example, comparison with zero won't work.
        return all([self._F[e].__eq__(other._F[e]) for e in range(self._nE)])

    def __ne__(self,other):
        r"""
        Test for non-equality with another cocycle. Two cocycles
        are non-equal if they take the different values on at least
        one edge in a fundamental domain.

        EXAMPLES::

            sage: X = BTQuotient(3,11)
            sage: H = HarmonicCocycles(X,4,prec=10)
            sage: v1 = H.basis()[0]; v2 = H.basis()[1]
            sage: v1 != v2 # indirect doctest
            True
            sage: v1 != v1
            False
        """
        return any([self._F[e].__ne__(other._F[e]) for e in range(self._nE)])

    def __nonzero__(self):
        r"""
        Test for being non-zero.

        EXAMPLES::

            sage: X = BTQuotient(3,13)
            sage: H = HarmonicCocycles(X,2,prec=10)
            sage: v1 = H.basis()[0]
            sage: v1.__nonzero__()
            True
            sage: v1 == 0
            False
            sage: v1 == H(0)
            False
        """
        tmp = [0 for nnn in range(parent(self)._X._num_edges)]
        return any([self._F[e]==tmp for e in range(self._nE)])

    def valuation(self):
        r"""
        Returns the valuation of the cocycle, defined as the
        minimum of the values it takes on a set of representatives.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: X = BTQuotient(3,17)
            sage: H = HarmonicCocycles(X,2,prec=10)
            sage: b1 = H.basis()[0]
            sage: b2 = 3*b1
            sage: b1.valuation()
            0
            sage: b2.valuation()
            1
            sage: H(0).valuation()
            Infinity
        """
        if self.is_zero():
            return Infinity
        else:
            return min([self._F[e].valuation() for e in range(self._nE)])

    def _compute_element(self):
        r"""

        """
        R = self._R
        A = self.parent().basis_matrix().transpose()
        B = Matrix(R,self._nE*(self.parent()._k-1),1,[self._F[e].moment(ii)  for e in range(self._nE) for ii in range(self.parent()._k-1) ])
        res = (A.solve_right(B)).transpose()
        return self.parent().free_module()(res.row(0))

    #In HarmonicCocycle
    def evaluate(self,e1):
        r"""
        Evaluates a harmonic cocycle on an edge of the Bruhat-Tits tree.

        INPUT:

        - ``e1`` - a matrix corresponding to an edge of the Bruhat-Tits tree

        OUTPUT:

        - An element of the coefficient module of the cocycle which describes 
          the value of the cocycle on e1

        EXAMPLES:
        """
        X = self.parent()._X
        p = X._p
        u = DoubleCosetReduction(X,e1)
        if u.label < self._nE:
            val  =  self._F[u.label]
        else:
            val  =  -self._F[u.label-self._nE]

        # return val.l_act_by(u.igamma(self.parent().embed_quaternion) * (p**(-u.power)))
        return (u.igamma(self.parent().embed_quaternion) * (p**(-u.power))) * val
        #return (u.sign()*self._F[u.label])

    #In HarmonicCocycle
    def riemann_sum(self,f,center = 1,level = 0,E = None):
        r"""
        This function evaluates the integral of the funtion ``f`` with respect to the measure determined by ``self``.

        EXAMPLES::
        """
        R1 = LaurentSeriesRing(f.base_ring(),'r1')
        R1.set_default_prec(self.parent()._k-1)
        R2 = PolynomialRing(f.base_ring(),'r2')

        if E is None:
            E = self.parent()._X._BT.get_balls(center,level)
        else:
            E = self.parent()._X._BT.subdivide(E,level)
        value = 0
        ii = 0
        for e in E:
            ii += 1
            exp = R1([e[1,1],e[1,0]])**(self.parent()._k-2)*e.determinant()**(-(self.parent()._k-2)/2)*f(R1([e[0,1],e[0,0]])/R1([e[1,1],e[1,0]])).truncate(self.parent()._k-1)
            # new = self.evaluate(e).l_act_by(e.inverse()).evaluate(exp)
            new = (e.inverse() * self.evaluate(e) ).evaluate(exp)
            value += new
        return value

    def modular_form(self,z = None,level = 0):
        r"""

        EXAMPLES::
        """
        return self.derivative(z,level,order = 0)

    # In HarmonicCocycle
    def derivative(self,z = None,level = 0,order = 1):
        r"""
        The derivative of the modular form attached to ``self``.

        EXAMPLES::
        """
        def F(z):
            R = PolynomialRing(z.parent(),'x,y').fraction_field()
            Rx = PolynomialRing(z.parent(),'x1').fraction_field()
            x1 = Rx.gen()
            subst = R.hom([x1,z],codomain = Rx)
            x,y = R.gens()
            center = self.parent()._X._BT.find_containing_affinoid(z)
            zbar = z.trace()-z
            f = R(1)/(x-y)
            k = self.parent()._k
            V = [f]
            for ii in range(order):
                V = [v.derivative(y) for v in V]+[k/(y-zbar)*v for v in V]
                k += 2
            return sum([self.riemann_sum(subst(v),center,level) for v in V])
        if(z is None):
            return F
        else:
            return F(z)


class HarmonicCocycles(AmbientHeckeModule):
    Element = HarmonicCocycleElement
    r"""
    Represents a space of Gamma invariant harmonic
    cocycles valued in a cofficient module.

    INPUT:

    - ``X`` - A BTQuotient object

    - ``k`` - integer - The weight. It must be even.

    - ``prec`` - integer (Default: None). If specified, the precision
      for the coefficient module

    - ``basis_matrix`` - integer (Default: None)

    - ``base_field`` - (Default: None)

    EXAMPLES::

        sage: X = BTQuotient(3,23)
        sage: H = HarmonicCocycles(X,2,prec = 5)
        sage: H.dimension()
        3
        sage: X.genus()

    Higher even weights are implemented::

        sage: H = HarmonicCocycles(X,8, prec = 10)
        sage: H.dimension()
        26

    AUTHORS:

    - Cameron Franc (2012-02-20)
    - Marc Masdeu
    """
    def __init__(self,X,k,prec = None,basis_matrix = None,base_field = None):
        """
        Compute the space of harmonic cocycles.

        EXAMPLES:
            sage: X = BTQuotient(3,37)
            sage: H = HarmonicCocycles(X,4,prec=10)
            sage: TestSuite(H).run()
        """
        self._k = k
        self._X = X
        self._E = self._X.get_edge_list()
        self._V = self._X.get_vertex_list()

        if base_field is not None and not base_field.is_exact():
            prec = base_field.precision_cap()

        if prec is None:
            if base_field is None:
                try:
                    self._R =  X.get_splitting_field()
                except AttributeError:
                    raise ValueError, "It looks like you are not using Magma as backend...and still we don't know how to compute splittings in that case!"
            else:
                pol = X.get_splitting_field().defining_polynomial().factor()[0][0]
                self._R = base_field.extension(pol,pol.variable_name()).absolute_field(name = 'r')
            # self._U = OCVn(self._k-2,self._R)
            self._U = Symk(self._k-2,base = self._R,act_on_left = True,tuplegen = _btquot_tuplegen(),character = (None,-ZZ((self._k-2)/2))) #monoid = MatrixSpace(self._R,2,2))
        else:
            self._prec = prec
            if base_field is None:
                self._R = Qp(self._X._p,prec = prec)
            else:
                self._R = base_field
            # self._U = OCVn(self._k-2,self._R,self._k-1)
            self._U = Symk(self._k-2,base = self._R,act_on_left = True,tuplegen = _btquot_tuplegen(),character = (None,-ZZ((self._k-2)/2)))

        if basis_matrix is None:
            self.__rank = self._X.dimension_harmonic_cocycles(self._k)
        else:
            self.__rank = basis_matrix.nrows()

        if basis_matrix is not None:
            self.__matrix = basis_matrix
            self.__matrix.set_immutable()
            assert self.__rank == self.__matrix.nrows()

        AmbientHeckeModule.__init__(self, self._R, self.__rank, self._X.prime()*self._X.Nplus()*self._X.Nminus(), weight = self._k)
        self._populate_coercion_lists_()

    def base_extend(self,base_ring):
        r"""
        Extends the base ring of the coefficient module.

        INPUT:

        - ``base_ring`` - a ring that has a coerce map from the
          current base ring

        OUTPUT:

        A new space of HarmonicCocycles with the base extended.

        EXAMPLES::

            sage: X = BTQuotient(3,19)
            sage: H = HarmonicCocycles(X,2,10)
            sage: H.base_ring()
            3-adic Field with capped relative precision 10
            sage: H1 = H.base_extend(Qp(3,prec=15))
            sage: H1.base_ring()
            3-adic Field with capped relative precision 15
        """
        if not base_ring.has_coerce_map_from(self.base_ring()):
            raise ValueError, "No coercion defined"
        else:
            return self.change_ring(base_ring)

    def change_ring(self, new_base_ring):
        r"""
        Changes the base ring of the coefficient module.

        INPUT:

        - ``new_base_ring'' - a ring that has a coerce map from the
          current base ring

        OUTPUT:

        New space of HarmonicCocycles with different base ring

        EXAMPLES::

            sage: X = BTQuotient(5,17)
            sage: H = HarmonicCocycles(X,2,10)
            sage: H.base_ring()
            5-adic Field with capped relative precision 10
            sage: H1 = H.base_extend(Qp(5,prec=15)) # indirect doctest
            sage: H1.base_ring()
            5-adic Field with capped relative precision 15
        """
        if not new_base_ring.has_coerce_map_from(self.base_ring()):
            raise ValueError, "No coercion defined"

        else:
            return self.__class__(self._X,self._k,prec = None,basis_matrix = self.basis_matrix().change_ring(new_base_ring),base_field = new_base_ring)

    def rank(self):
        r"""
        Returns the rank (dimension) of ``self``.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: X = BTQuotient(7,11)
            sage: H = HarmonicCocycles(X,2,prec = 10)
            sage: X.genus() == H.rank()
            True
            sage: H1 = HarmonicCocycles(X,4,prec = 10)
            sage: H1.rank()
            16
        """
        return self.__rank

    def submodule(self,v,check = False):
        r"""
        Return the submodule of ``self`` spanned by ``v``.

        INPUT:

        - ``v`` - Submodule of self.free_module().

        - ``check`` - Boolean (Default = False).

        OUTPUT:

        Subspace of harmonic cocycles.

        EXAMPLES::

            sage: X = BTQuotient(3,17)
            sage: H = HarmonicCocycles(X,2,prec=10)
            sage: H.rank()
            3
            sage: v = H.an_element()
            sage: N = H.free_module().span([v.element()])
            sage: H1 = H.submodule(N)
            sage: H1.rank() == 1
            True
        """
        return HarmonicCocyclesSubmodule(self,v,dual = None,check = check)

    def is_simple(self):
        r"""
        Whether ``self`` is irreducible.

        OUTPUT:

        Boolean. True iff self is irreducible.

        EXAMPLES::

            sage: X = BTQuotient(3,29)
            sage: H = HarmonicCocycles(X,4,prec =10)
            sage: H.rank()
            14
            sage: H.is_simple()
            False
            sage: X = BTQuotient(7,2)
            sage: H = HarmonicCocycles(X,2,prec=10)
            sage: H.rank()
            1
            sage: H.is_simple()
            True
        """
        return self.rank() == 1

    def _repr_(self):
        r"""
        This returns the representation of self as a string.

        EXAMPLES::

            sage: X = BTQuotient(5,23)
            sage: H = HarmonicCocycles(X,2,prec=10)
            sage: H._repr_()
            Space of harmonic cocycles of weight 2 on Quotient of the Bruhat Tits tree of GL_2(QQ_5) with discriminant 23 and level 1
        """
        return 'Space of harmonic cocycles of weight %s on %s'%(self._k,self._X)

    def _latex_(self):
        r"""
        A LaTeX representation of ``self``.

        EXAMPLES::

            sage: X = BTQuotient(5,23)
            sage: H = HarmonicCocycles(X,2,prec=10)
            \text{Space of harmonic cocycles of weight } 2 \text{ on } X(5 \cdot 23,1)\otimes_{\mathbb{Z}} \mathbb{F}_{5}
        """
        s = '\\text{Space of harmonic cocycles of weight }'+latex(self._k)+'\\text{ on }'+latex(self._X)
        return s

    def _an_element_(self):
        r"""
        Returns an element of the ambient space

        OUTPUT:

        A harmonic cocycle in self.

        EXAMPLES:

            sage: X = BTQuotient(5,23)
            sage: H = HarmonicCocycles(X,2,prec=10)
            sage: H.an_element()
            Element of Space of harmonic cocycles of weight 2 on Quotient of the Bruhat Tits tree of GL_2(QQ_5) with discriminant 23 and level 1
        """
        return self.basis()[0]


    def _coerce_map_from_(self, S):
        r"""
        Can coerce from other HarmonicCocycles or from pAutomorphicForms

        OUTPUT:

        Boolean. True iff self is a space of HarmonicCocycles or
        pAutomorphicForms.

        EXAMPLES::

            sage: X = BTQuotient(3,17)
            sage: H = HarmonicCocycles(X,2,prec=10)
            sage: v = H.an_element()
            sage: N = H.free_module().span([v.element()])
            sage: H1 = H.submodule(N)
            sage: H1._coerce_map_from_(H)
            True
            sage: H._coerce_map_from_(H1)
            True
        """
        if isinstance(S,(HarmonicCocycles,pAutomorphicForms)):
            if S._k != self._k:
                return False
            if S._X != self._X:
                return False
            return True
        return False

    def __cmp__(self,other):
        r"""

        """
        try:
            res = (self.base_ring() == other.base_ring() and self._X == other._X and self._k == other._k)
            return res
        except:
            return False

    def _element_constructor_(self,x):
        r"""

        """
        #Code how to coherce x into the space
        #Admissible values of x?
        if isinstance(x,HarmonicCocycleElement):
            return HarmonicCocycleElement(self,x)
        elif isinstance(x,pAutomorphicFormElement):
            # tmp = [self._U(x._F[ii]).l_act_by(self._E[ii].rep) for ii in range(self._nE)]
            tmp = [self._E[ii].rep * self._U(x._F[ii]) for ii in range(self._nE)]
            return HarmonicCocycleElement(self,tmp,from_values = True)
        elif x == 0:
            tmp = [[0 for jj in range(self.weight()-1)] for ii in range(self._X._num_edges)]
            return HarmonicCocycleElement(self,tmp,from_values = True)

        #This will error if you try to coerce something illegal
        else:
            return HarmonicCocycleElement(self,x)


    def free_module(self):
        r"""
        This function returns the underlying free module

        EXAPLES::
        """
        try: return self.__free_module
        except AttributeError: pass
        V = self.base_ring()**self.dimension()
        self.__free_module = V
        return V

    def character(self):
        r"""
        Only implemented the trivial character so far.

        EXAMPLES::

        """
        return lambda x:x

    def embed_quaternion(self,g):
        r"""
        Embed the quaternion element ``g`` into the matrix algebra.

        EXAMPLES::
        """
        tmp = self._X.embed_quaternion(g,exact = self._R.is_exact(), prec = self._prec)
        tmp.set_immutable()
        return tmp

    def basis_matrix(self):
        r"""
        Returns a basis of ``self`` in matrix form.

        If the coefficient module `M` is of finite rank then the space
        of Gamma invariant `M` valued harmonic cocycles can be
        represented as a subspace of the finite rank space of all
        functions from the finitely many edges in the corresponding
        BTQuotient into `M`. This function computes this
        representation of the space of cocycles.

        OUTPUT:

        - A basis matrix describing the cocycles in the spaced of all
          `M` valued Gamma invariant functions on the tree.

        EXAMPLES::

            sage: X = BTQuotient(3,19)
            sage: C = HarmonicCocycles(X,4,prec = 5)
            sage: B = C.basis()
            Traceback (most recent call last):
            ...
            RuntimeError: The computed dimension does not agree with the expectation. Consider increasing precision!

        We try increasing the precision:

        ::

            sage: C = HarmonicCocycles(X,4,prec = 20)
            sage: B = C.basis()
            sage: len(B) == X.dimension_harmonic_cocycles(4)
            True

        AUTHORS:

        - Cameron Franc (2012-02-20)
        - Marc Masdeu (2012-02-20)
        """
        try: return self.__matrix
        except AttributeError: pass
        nV = len(self._V)
        nE = len(self._E)
        stab_conds = []
        S = self._X.get_edge_stabs()
        p = self._X._p
        d = self._k-1
        for e in self._E:
            try:
                g = filter(lambda g:g[2],S[e.label])[0]
                # C = self._U.l_matrix_representation(self.embed_quaternion(g[0]))
                # C -= self._U.l_matrix_representation(Matrix(QQ,2,2,p**g[1]))
                C = self._U.acting_matrix(self.embed_quaternion(g[0]),d).transpose()
                C -= self._U.acting_matrix(Matrix(QQ,2,2,p**g[1]),d).transpose()
                stab_conds.append([e.label,C])
            except IndexError: pass

        n_stab_conds = len(stab_conds)
        self._M = Matrix(self._R,(nV+n_stab_conds)*d,nE*d,0,sparse = True)
        for v in self._V:
            for e in filter(lambda e:e.parity == 0,v.leaving_edges):
                # C = sum([self._U.l_matrix_representation(self.embed_quaternion(x[0])) for x in e.links],Matrix(self._R,d,d,0))
                C = sum([self._U.acting_matrix(self.embed_quaternion(x[0]),d) for x in e.links],Matrix(self._R,d,d,0)).transpose()
                self._M.set_block(v.label*d,e.label*d,C)
            for e in filter(lambda e:e.parity == 0,v.entering_edges):
                # C = sum([self._U.l_matrix_representation(self.embed_quaternion(x[0])) for x in e.opposite.links],Matrix(self._R,d,d,0))
                C = sum([self._U.acting_matrix(self.embed_quaternion(x[0]),d) for x in e.opposite.links],Matrix(self._R,d,d,0)).transpose()
                self._M.set_block(v.label*d,e.opposite.label*d,C)

        for kk in range(n_stab_conds):
            v = stab_conds[kk]
            self._M.set_block((nV+kk)*d,v[0]*d,v[1])

        x1 = self._M.right_kernel().matrix()

        if x1.nrows() !=  self.rank():
            raise RuntimeError, 'The computed dimension does not agree with the expectation. Consider increasing precision!'

        K = [c for c in x1.rows()]

        if not self._R.is_exact():
            for ii in range(len(K)):
                s = min([t.valuation() for t in K[ii]])
                for jj in range(len(K[ii])):
                    K[ii][jj] = (p**(-s))*K[ii][jj]

        self.__matrix = Matrix(self._R,len(K),nE*d,K)
        self.__matrix.set_immutable()
        return self.__matrix

    def __apply_atkin_lehner(self,q,f):
        r"""
        This function applies an Atkin-Lehner involution to a
        harmonic cocycle

        INPUT:

        - ``q`` - an integer dividing the full level p*Nminus*Nplus

        - ``f`` - a harmonic cocycle

        OUTPUT:

        - The harmonic cocycle obtained by hitting f with the
          Atkin-Lehner at q

        EXAMPLES::
        """
        R = self._R
        Data = self._X._get_atkin_lehner_data(q)
        p = self._X._p
        tmp = [self._U(0) for jj in range(len(self._E))]
        d1 = Data[1]
        mga = self.embed_quaternion(Data[0])
        nE = len(self._E)
        for jj in range(nE):
            t = d1[jj]
            if t.label < nE:
                # tmp[jj] += (f._F[t.label]).l_act_by(p**(-t.power)*mga*t.igamma(self.embed_quaternion))
                tmp[jj] += (p**(-t.power)*mga*t.igamma(self.embed_quaternion)) * f._F[t.label]
            else:
                # tmp[jj] += (-f._F[t.label-nE]).l_act_by(p**(-t.power)*mga*t.igamma(self.embed_quaternion))
                tmp[jj] += (p**(-t.power)*mga*t.igamma(self.embed_quaternion)) * (-f._F[t.label-nE])


        return HarmonicCocycleElement(self,tmp,from_values = True)

    def __apply_hecke_operator(self,l,f):
        r"""
        This function applies a Hecke operator to a harmonic cocycle.

        INPUT:

        - ``l`` - an integer

        - ``f`` - a harmonic cocycle

        OUTPUT:

        - A harmonic cocycle which is the result of applying the lth
          Hecke operator to f

        EXAMPLES::

        """
        R = self._R
        HeckeData,alpha = self._X._get_hecke_data(l)
        if(self.level()%l == 0):
            factor = QQ(l**(Integer((self._k-2)/2))/(l+1))
        else:
            factor = QQ(l**(Integer((self._k-2)/2)))
        p = self._X._p
        alphamat = self.embed_quaternion(alpha)
        tmp = [self._U(0) for jj in range(len(self._E))]
        for ii in range(len(HeckeData)):
            d1 = HeckeData[ii][1]
            mga = self.embed_quaternion(HeckeData[ii][0])*alphamat
            nE = len(self._E)
            for jj in range(nE):
                t = d1[jj]
                if t.label < nE:
                    # tmp[jj] += f._F[t.label].l_act_by(p**(-t.power)*mga*t.igamma(self.embed_quaternion))
                    tmp[jj] += (p**(-t.power)*mga*t.igamma(self.embed_quaternion)) * f._F[t.label]
                else:
                    # tmp[jj] += (-f._F[t.label-nE]).l_act_by(p**(-t.power)*mga*t.igamma(self.embed_quaternion))
                    tmp[jj] += (p**(-t.power)*mga*t.igamma(self.embed_quaternion)) * (-f._F[t.label-nE])

        return HarmonicCocycleElement(self,[factor*x for x in tmp],from_values = True)

    def _compute_atkin_lehner_matrix(self,d):
        r"""
        When the underlying coefficient module is finite, this
        function computes the matrix of an Atkin-Lehner involution in
        the basis provided by the function basis_matrix

        INPUT:

        - ``d`` - an integer dividing p*Nminus*Nplus

        OUTPUT:

        - The matrix of the AL-involution at d in the basis given by
          self.basis_matrix

        EXAMPLES::

        """
        res = self.__compute_operator_matrix(lambda f:self.__apply_atkin_lehner(d,f))
        return res

    def _compute_hecke_matrix_prime(self,l):
        r"""
        When the underlying coefficient module is finite, this
        function computes the matrix of a (prime) Hecke operator in
        the basis provided by the function basis_matrix

        INPUT:

        - ``l`` - an integer prime

        OUTPUT:

        - The matrix of `T_l` acting on the cocycles in the basis given by
          self.basis_matrix

        EXAMPLES::

        """
        res = self.__compute_operator_matrix(lambda f:self.__apply_hecke_operator(l,f))
        return res

    def __compute_operator_matrix(self,T):
        r"""
        Compute the matrix of the operator `T`.

        EXAMPLES::

        """
        R = self._R
        A = self.basis_matrix().transpose()
        basis = self.basis()
        B = zero_matrix(R,len(self._E) * (self._k-1),self.dimension())
        for rr in range(len(basis)):
            g = T(basis[rr])
            B.set_block(0,rr,Matrix(R,len(self._E) * (self._k-1),1,[g._F[e].moment(ii)  for e in range(len(self._E)) for ii in range(self._k-1) ]))

        res = (A.solve_right(B)).transpose()
        res.set_immutable()
        return res

class HarmonicCocyclesSubmodule(sage.modular.hecke.submodule.HeckeSubmodule,HarmonicCocycles):
    r"""
    Submodule of a space of HarmonicCocycles.

    INPUT:

    - ``x`` - integer (default: 1) the description of the
      argument x goes here.  If it contains multiple lines, all
      the lines after the first need to be indented.

    - ``y`` - integer (default: 2) the ...

    EXAMPLES::

        sage: X = BTQuotient(3,17)
        sage: H = HarmonicCocycles(X,2,prec=10)
        sage: N = H.free_module().span([H.an_element().element()])
        sage: H1 = H.submodule(N) # indirect doctest
        sage: H1
        Subspace of Space of harmonic cocycles of weight 2 on Quotient of the Bruhat Tits tree of GL_2(QQ_3) with discriminant 17 and level 1 of dimension 1

    AUTHOR:

    - Marc Masdeu (2012-02-20)
    """
    def __init__(self, ambient_module, submodule, dual = None, check = False):
        """
        Submodule of harmonic cocycles.

        INPUT:

        - `ambient_module` - HarmonicCocycles

        - `submodule` - submodule of the ambient space.

        - `dual_module` - (default: None) ignored

        - `check` - (default: False) whether to check that the
          submodule is Hecke equivariant

        EXAMPLES::

            sage: X = BTQuotient(3,17)
            sage: H = HarmonicCocycles(X,2,prec=10)
            sage: N = H.free_module().span([H.an_element().element()])
            sage: H1=H.submodule(N)
            sage: TestSuite(H1).run()
        """
        A = ambient_module
        sage.modular.hecke.submodule.HeckeSubmodule.__init__(self, A, submodule, check = check)
        self.__rank = submodule.dimension()
        HarmonicCocycles.__init__(self,X = A._X,k = A._k,prec = A._prec,basis_matrix = submodule.basis_matrix()*A.basis_matrix())

    def rank(self):
        r"""
        Returns the rank (dimension) of the submodule.

        OUTPUT:

        Integer - The rank of ``self``.

        EXAMPLES::

            sage: X = BTQuotient(3,17)
            sage: H = HarmonicCocycles(X,2,prec=10)
            sage: N = H.free_module().span([H.an_element().element()])
            sage: H1=H.submodule(N)
            sage: H1.rank()
            1
        """
        return self.__rank

    def _repr_(self):
        r"""
        Returns the representation of self as a string.

        OUTPUT:

        String representation of self.

        EXAMPLES::

            sage: X = BTQuotient(3,17)
            sage: H = HarmonicCocycles(X,2,prec=10)
            sage: N = H.free_module().span([H.an_element().element()])
            sage: H1=H.submodule(N)
            sage: H1._repr_()
            'Subspace of Space of harmonic cocycles of weight 2 on Quotient of the Bruhat Tits tree of GL_2(QQ_3) with discriminant 17 and level 1 of dimension 1'
        """
        return "Subspace of %s of dimension %s"%(self.ambient(),self.dimension())


class pAutomorphicFormElement(ModuleElement):
    r"""
    Rudimentary implementation of a class for a p-adic
    automorphic form on a definite quaternion algebra over Q. These
    are required in order to compute moments of measures associated to
    harmonic cocycles on the BT-tree using the overconvergent modules
    of Darmon-Pollack and Matt Greenberg. See Greenberg's thesis for
    more details.

    INPUT:

    - ``vec`` - Quite flexible input
    - ``quick`` - boolean (default: False)

    EXAMPLES::

    REFERENCES:

    Matthew Greenberg's thesis (available on his webpage as of 02/12).

    AUTHORS:

    - Cameron Franc (2012-02-20)
    - Marc Masdeu

    """
    def __init__(self,parent,vec,quick = False):
        ModuleElement.__init__(self,parent)
        self._num_generators = len(parent._list)
        self._cached_values = dict()
        self._R = Qp(parent.prime(),prec = parent._prec)
        if quick:
            self._value = [ parent._U(v) for v in vec ]
        else:
            if isinstance(vec,pAutomorphicFormElement):
                self._value = parent._make_invariant([parent._U(vec._value[ii]) for ii in range(self._num_generators)])

            elif isinstance(vec,HarmonicCocycleElement):
                assert(parent._U.weight() == vec._wt-2)
                F = []
                assert(2*len(vec._F) == self._num_generators)
                # assert(isinstance(parent._U,OCVn))
                E = parent._list

                tmp = []
                for ii in range(len(vec._F)):
                    assert len(vec._F[ii].moments) > 0
                    # newtmp = vec.parent()(vec._F[ii]).l_act_by(E[ii].rep.inverse())
                    newtmp = (E[ii].rep.inverse()) * vec._F[ii]
                    tmp.append(newtmp)
                    F.append(newtmp)
                A = Matrix(QQ,2,2,[0,-1/parent.prime(),-1,0])
                for ii in range(len(vec._F)):
                    # F.append(parent._U(-1*tmp[ii].r_act_by(A))) #TBC
                    newF = - act_left(A.adjoint(),tmp[ii])
                    F.append(newF)
                self._value = parent._make_invariant(F)

            elif(isinstance(vec,list) and len(vec) == self._num_generators):
                try:
                    self._value = [parent._U(v) for v in vec]
                except:
                    try:
                        veczp = parent._U.base_ring()(vec)

                        self._value = [parent._U(veczp) for ii in range(self._num_generators)]
                    except:
                        print vec
                        assert 0
            else:
                try:
                    veczp = parent._U.base_ring()(vec)
                    self._value = [parent._U(veczp) for ii in range(self._num_generators)]
                except:
                    raise ValueError,"Cannot initialize a p-adic automorphic form with the given input = "+str(vec)


    def precision(self):
        r"""
        The precision of ``self``, which is the minimum among the
        precision of the values on a fundamental domain.

        EXAMPLES::

        """
        return min(x.precision() for x in self._value)

    def _add_(self,g):
        r"""
        This function adds two p-adic automorphic forms.

        INPUT:

        - ``g`` - a p-adic automorphic form

        OUTPUT:

        - the result of adding g to self

        EXAMPLES::
        """
        #Should ensure that self and g are of the same weight and on the same curve
        vec = [self._value[e]+g._value[e] for e in range(self._num_generators)]
        return pAutomorphicFormElement(self.parent(),vec,quick = True)

    def _sub_(self,g):
        r"""
        This function subtracts a p-adic automorphic form from another.

        INPUT:

        - ``g`` - a p-adic automorphic form

        OUTPUT:

        - the result of subtracting g from self

        EXAMPLES::

        """
        #Should ensure that self and g are of the same weight and on the same curve
        vec = [self._value[e]-g._value[e] for e in range(self._num_generators)]
        return pAutomorphicFormElement(self.parent(),vec,quick = True)

    def _getitem_(self,e1):
        r"""
        Evaluates a p-adic automorphic form on a matrix in GL2(Qp).

        INPUT:

        - ``e1`` - a matrix in GL2(Qp)

        OUTPUT:

        - the value of self evaluated on e1

        EXAMPLES::

        """
        return self.evaluate(e1)

    def evaluate(self,e1):
        r"""
        Evaluates a p-adic automorphic form on a matrix in GL2(Qp).

        INPUT:

        - ``e1`` - a matrix in GL2(Qp)

        OUTPUT:

        - the value of self evaluated on e1

        EXAMPLES::

        """
        X = self.parent()._source
        p = self.parent().prime()
        u = DoubleCosetReduction(X,e1)
        # return self._value[u.label].r_act_by((u.t())*p**(u.power)) #TBC
        return act_left(((u.t(self.precision_cap()+1))*p**(u.power)).adjoint(),self._value[u.label])

    def _rmul_(self,a):
        r"""
        Multiplies the automorphic form by a scalar.

        EXAMPLES::

        """
        #Should ensure that 'a' is a scalar
        return pAutomorphicFormElement(self.parent(),[a*self._value[e] for e in range(self._num_generators)],quick = True)

    def _repr_(self):
        r"""
        This returns the representation of self as a string.

        EXAMPLES::

        """
        tmp = 'p-adic automorphic form on '+str(self.parent())+':\n'
        tmp += '   e   |   c(e)'
        tmp += '\n'
        for e in range(Integer(self._num_generators/2)):
            tmp += '  '+str(e)+' | '+str(self._value[e])+'\n'
        return tmp

    def valuation(self):
        r"""
        The valuation of ``self``, defined as the minimum of the
        valuations of the values that it takes on a set of edge
        representatives.

        EXAMPLES::

        """
        if self == 0:
            return Infinity
        else:
            return(min([self._value[e].valuation() for e in range(self._num_generators)]))

    def __nonzero__(self):
        r"""
        """
        return(any([self._value[e].__nonzero__() for e in range(self._num_generators)]))

    def improve(self):
        r"""
        Repeatedly applies the `U_p` operator to a p-adic
        automorphic form. This is used to compute moments of a measure
        associated to a rigid modular form in the following way: lift
        a rigid modular form to an ``overconvergent'' `p`-adic
        automorphic form in any way, and then repeatedly apply `U_p`
        to project to the ordinary part.  The resulting form encodes
        the moments of the measure of the original rigid modular form
        (assuming it is ordinary).


        EXAMPLES::


        REFERENCES:

        For details see Matthew Greenberg's thesis (available on his
        webpage as of 02/12).  Alternatively check out Darmon-Pollack
        for the analogous algorithm in the case of modular symbols.

        AUTHORS:

        - Cameron Franc (2012-02-20)
        - Marc Masdeu

        """
        MMM = self.parent()
        h2 = MMM._apply_Up_operator(self,True)
        verbose("Applied Up 1 time")
        ii = 0
        current_val = 0
        old_val = -Infinity
        init_val = self.valuation()
        while(current_val>old_val):
            old_val = current_val
            ii += 1
            self._value = [self.parent()._U(c) for c in h2._value]
            h2 = MMM._apply_Up_operator(self,scale = True)
            current_val = (h2-self).valuation()-init_val
            # print 'val  = ',current_val
            if current_val is Infinity:
                break
            verbose('Applied Up %s times'%ii+1)
        self._value = [self.parent()._U(c) for c in h2._value]

    def integrate(self,f,center = 1,level = 0,method = 'moments'):
        r"""
        Calculate
        .. MATH::

            \int_{\PP^1(\QQ_p)} f(x)d\mu(x)

        were `\mu` is the measure associated to ``self``.

        INPUT:

        - ``f`` - An analytic function.

        - ``center`` - 2x2 matrix over Qp (default: 1)

        - ``level`` - integer (default: 0)

        - ``method`` - string (default: 'moments'). Which method of
          integration to use. Either 'moments' or 'riemann_sum'.


        EXAMPLES::

        AUTHORS:

        - Marc Masdeu (2012-02-20)
        - Cameron Franc (2012-02-20)

        """
        E = self.parent()._source._BT.get_balls(center,level)
        R1 = LaurentSeriesRing(f.base_ring(),'r1')
        R2 = PolynomialRing(f.base_ring(),'r2')
        value = 0
        ii = 0
        if(method == 'riemann_sum'):
            R1.set_default_prec(self.parent()._U.weight()+1)
            for e in E:
                ii += 1
                #print ii,"/",len(E)
                exp = ((R1([e[1,1],e[1,0]]))**(self.parent()._U.weight())*e.determinant()**(-(self.parent()._U.weight())/2))*f(R1([e[0,1],e[0,0]])/R1([e[1,1],e[1,0]]))
                #exp = R2([tmp[jj] for jj in range(self.parent()._k-1)])
                new = self.evaluate(e).evaluate(exp.truncate(self.parent()._U.weight()+1))
                value += new
        elif(method == 'moments'):
            R1.set_default_prec(self.parent()._U.precision_cap())
            for e in E:
                ii += 1
                #print ii,"/",len(E)
                tmp = ((R2([e[1,1],e[1,0]]))**(self.parent()._U.weight())*e.determinant()**(-(self.parent()._U.weight())/2))*f(R2([e[0,1],e[0,0]])/R2([e[1,1],e[1,0]]))
                exp = R1(tmp.numerator())/R1(tmp.denominator())
                new = self.evaluate(e).evaluate(exp)
                value += new
        else:
            print 'The available methods are either "moments" or "riemann_sum". The latter is only provided for consistency check, and should never be used.'
            return False
        return value

    def modular_form(self,z = None,level = 0,method = 'moments'):
        r"""
        Returns the modular form corresponding to ``self``.

        INPUT:

        - ``z`` - (default: None). If specified, returns the value of
          the form at the point ``zz`` in the `p`-adic upper half
          plane.

        - ``level`` - integer (default: 0). If ``method`` is
          'riemann_sum', will use a covering of `\PP^1(\QQ_p)` with
          balls of size `p^-\mbox{level]`.

        - ``method`` - string (default: ``moments``). It must be
          either ``moments`` or ``riemann_sum``.

        OUTPUT:

        - A function from the `p`-adic upper half plane to `\CC_p`. If
          an argument ``z`` was passed, returns instead the value at
          that point.

        """
        return self.derivative(z,level,method,order = 0)

    def derivative(self,z = None,level = 0,method = 'moments',order = 1):
        r"""
        Returns the derivative of the modular form corresponding to
        ``self``.

        INPUT:

        - ``z`` - (Default: None). If specified, evaluates the derivative
           at the point ``z`` in the `p`-adic upper half plane.

        - ``level`` - integer (default: 0). If ``method`` is
          'riemann_sum', will use a covering of `\PP^1(\QQ_p)` with
          balls of size `p^-\mbox{level]`.

        - ``method`` - string (default: ``moments``). It must be
          either ``moments`` or ``riemann_sum``.

        - ``order`` - integer (Default: 1). The order of the
          derivative to be computed.

        OUTPUT:

        - A function from the `p`-adic upper half plane to `\CC_p`. If
          an argument ``z`` was passed, returns instead the value of
          the derivative at that point.

        EXAMPLES::

        """
        def F(z,level = 0,method = 'moments'):
            R = PolynomialRing(z.parent(),'x,y').fraction_field()
            Rx = PolynomialRing(z.parent(),'x1').fraction_field()
            x1 = Rx.gen()
            subst = R.hom([x1,z],codomain = Rx)
            x,y = R.gens()
            center = self.parent()._source._BT.find_containing_affinoid(z)
            zbar = z.trace()-z
            f = R(1)/(x-y)
            k = self.parent()._n+2
            V = [f]
            for ii in range(order):
                V = [v.derivative(y) for v in V]+[k/(y-zbar)*v for v in V]
                k += 2
            return sum([self.integrate(subst(v),center,level,method) for v in V])
        if z is None:
            return F

        try:
            try: tmp = self._cached_values[z,level,method,order,self.precision()]
            except KeyError:
                tmp = F(z,level,method)
                self._cached_values[z,level,method,order,self.precision()] = tmp
            return tmp
        except TypeError:
            return F(z,level,method)


    # So far we can't break it into two integrals because of the pole at infinity.
    def coleman(self,t1,t2,E = None,method = 'moments',mult = False,delta = -1,level = 0):
        r"""
        If ``self`` is a `p`-adic automorphic form that
        corresponds to a rigid modular form, then this computes the
        coleman integral of this form between two points on the
        boundary `\PP^1(\QQ_p)` of the `p`-adic upper half plane.

        INPUT:

        - ``t1``, ``t2`` - elements of `\PP^1(\QQ_p)` (the endpoints
          of integration)

        - ``E`` - (Default: None). If specified, will not compute the
           covering adapted to ``t1`` and ``t2`` and instead use the
           given one. In that case, ``E`` should be a list of matrices
           corresponding to edges describing the open balls to be
           considered.

        - ``method`` - string (Default: 'moments'). Tells which
          algorithm to use (alternative is 'riemann_sum', which is
          unsuitable for computations requiring high precision)

        - ``mult`` - boolean (Default: False). Whether to use the
          multiplicative version.

        - ``delta`` - integer (Default: -1)

        - ``level`` - integer (Default: 0)

        OUTPUT:

          The result of the coleman integral

        EXAMPLES::

            sage: p = 7
            sage: lev = 2
            sage: prec = 20
            sage: X = BTQuotient(p,lev, use_magma = True) # optional - magma
            sage: k = 2 # optional - magma
            sage: M = HarmonicCocycles(X,k,prec) # optional - magma
            sage: B = M.basis() # optional - magma
            sage: f = 3*B[0] # optional - magma
            sage: MM = pAutomorphicForms(X,k,prec,overconvergent = True) # optional - magma
            sage: D = -11 # optional - magma
            sage: X.is_admissible(D) # optional - magma
            True
            sage: K.<a> = QuadraticField(D) # optional - magma
            sage: CM = X.get_CM_points(D,prec = prec) # optional - magma
            sage: Kp = CM[0].parent() # optional - magma
            sage: P = CM[0] # optional - magma
            sage: Q = P.trace()-P # optional - magma
            sage: F = MM.lift(f) # long time optional - magma
            sage: J0 = F.coleman(P,Q,mult = True) # long time optional - magma
            sage: E = EllipticCurve([1,0,1,4,-6]) # optional - magma
            sage: T = E.tate_curve(p) # optional - magma
            sage: xx,yy = getcoords(T,J0,prec) # long time optional -magma
            sage: P = E.base_extend(K).lift_x(algdep(xx,1).roots(QQ)[0][0]); P # long time optional - magma
            (7/11 : 58/121*a - 9/11 : 1)

            sage: p = 13 # optional - magma
            sage: lev = 2 # optional - magma
            sage: prec = 20 # optional - magma
            sage: Y = BTQuotient(p,lev, use_magma = True) # optional - magma
            sage: k = 2 # optional - magma
            sage: M = HarmonicCocycles(Y,k,prec) # optional - magma
            sage: B = M.basis() # optional - magma

            sage: f = B[1] # optional - magma
            sage: g = -4*B[0]+3*B[1] # optional - magma
            sage: MM = pAutomorphicForms(Y,k,prec,overconvergent = True) # optional - magma
            sage: D = -11 # optional - magma
            sage: Y.is_admissible(D) # optional - magma
            True
            sage: K.<a> = QuadraticField(D) # optional - magma
            sage: CM = Y.get_CM_points(D,prec = prec) # optional - magma
            sage: Kp = parent(CM[0]) # optional - magma
            sage: P = CM[0] # optional - magma
            sage: Q = P.trace()-P # optional - magma
            sage: F = MM.lift(f) # long time optional - magma
            sage: J11 = F.coleman(P,Q,mult = True) # long time optional - magma
            sage: E = EllipticCurve('26a2') # optional - magma
            sage: T = E.tate_curve(p) # optional - magma
            sage: xx,yy = getcoords(T,J11,prec) # long time optional - magma
            sage: HP = E.base_extend(K).lift_x(algdep(xx,1).roots(QQ)[0][0]); HP # long time optional - magma
            (-137/11 : 2/121*a + 63/11 : 1)

        AUTHORS:

        - Cameron Franc (2012-02-20)
        - Marc Masdeu (2012-02-20)
        """
        if(mult and delta >= 0):
            raise NotImplementedError, "Need to figure out how to implement the multiplicative part."
        p = self.parent().prime()
        K = t1.parent()
        R = PolynomialRing(K,'x')
        x = R.gen()
        R1 = LaurentSeriesRing(K,'r1')
        r1 = R1.gen()
        if(E is None):
            E = self.parent()._source._BT.find_covering(t1,t2)
            # print 'Got %s open balls.'%len(E)
        value = 0
        ii = 0
        value_exp = K(1)
        if(method == 'riemann_sum'):
            R1.set_default_prec(self.parent()._U.weight()+1)
            for e in E:
                ii += 1
                b = e[0,1]
                d = e[1,1]
                y = (b-d*t1)/(b-d*t2)
                poly = R1(y.log()) #R1(our_log(y))
                c_e = self.evaluate(e)
                new = c_e.evaluate(poly)
                value += new
                if mult:
                    value_exp  *=  K.teichmuller(y)**Integer(c_e[0].rational_reconstruction())

        elif(method == 'moments'):
            R1.set_default_prec(self.parent()._U.precision_cap())
            for e in E:
                ii += 1
                f = (x-t1)/(x-t2)
                a,b,c,d = e.list()
                y0 = f(R1([b,a])/R1([d,c])) #f( (ax+b)/(cx+d) )
                y0 = p**(-y0(0).valuation())*y0
                mu = K.teichmuller(y0(0))
                y = y0/mu-1
                poly = R1(0)
                ypow = y
                for jj in range(1,R1.default_prec()+10):
                    poly += (-1)**(jj+1)*ypow/jj
                    ypow *= y
                if(delta >= 0):
                    poly *= ((r1-t1)**delta*(r1-t2)**(self.parent()._n-delta))
                c_e = self.evaluate(e)
                new = c_e.evaluate(poly)
                value += new
                if mult:
                    value_exp  *=  K.teichmuller(((b-d*t1)/(b-d*t2)))**Integer(c_e[0].rational_reconstruction())

        else:
            print 'The available methods are either "moments" or "riemann_sum". The latter is only provided for consistency check, and should not be used in practice.'
            return False
        if mult:
            return K.teichmuller(value_exp) * value.exp()
        return value


class pAutomorphicForms(Module):
    Element = pAutomorphicFormElement
    r"""
    The module of (quaternionic) `p`-adic automorphic forms.

    INPUT:

    - `domain` - A BTQuotient.

    - `U` - A coefficient module or an integer. If U is a coefficient module then this creates the relevant space of automorphic forms. If U is an integer then the coefficients are the (`U-2`)nd power of the symmetric representation of  `\GL_2(\Qp)`.

    - `prec` - A precision (Default = None). If not None should be a
      positive integer

    - `t` - (Default = None).

    - `R` - (Default = None).

    - `overconvergent` - Boolean (Default = False).

    EXAMPLES::

        sage: X = BTQuotient(11,5)
        sage: H = HarmonicCocycles(X,2,prec=10)
        sage: A = pAutomorphicForms(X,H)

    AUTHORS:

    - Cameron Franc (2012-02-20)
    - Marc Masdeu (2012-02-20)
    """
    def __init__(self,domain,U,prec = None,t = None,R = None,overconvergent = False):
        if(R is None):
            if not isinstance(U,Integer):
                self._R = U.base_ring()
            else:
                if(prec is None):
                    prec = 100
                self._R = Qp(domain._p,prec)
        else:
            self._R = R
        #U is a CoefficientModuleSpace
        if isinstance(U,Integer):
            if t is None:
                if overconvergent:
                    t = prec-U+1
                else:
                    t = 0
            # self._U = OCVn(U-2,self._R,U-1+t)
            self._U = Distributions(U-2,base = self._R,prec_cap = U - 1 + t ,act_on_left = True,tuplegen = _btquot_tuplegen(),character = (None,-ZZ((U-2)/2))) #monoid = MatrixSpace(self._R,2,2))
        else:
            self._U = U
        self._source = domain
        self._list = self._source.get_list() # Contains also the opposite edges
        self._prec = self._R.precision_cap()
        self._n = self._U.weight()
        self._p = self._source._p
        Module.__init__(self,base = self._R)
        self._populate_coercion_lists_()

    def prime(self):
        return self._p

    def _repr_(self):
        r"""
        This returns the representation of self as a string.

        EXAMPLES::
        """
        s = 'Space of automorphic forms on '+str(self._source)+' with values in '+str(self._U)
        return s

    def _coerce_map_from_(self, S):
        r"""
        Can coerce from other HarmonicCocycles or from pAutomorphicForms
        """
        if isinstance(S,HarmonicCocycles):
            if S.weight()-2 != self._n:
                return False
            if S._source != self._source:
                return False
            return True
        if isinstance(S,pAutomorphicForms):
            if S._n != self._n:
                return False
            if S._source != self._source:
                return False
            return True
        return False

    def _element_constructor_(self,x):
        r"""
        """
        #Code how to coherce x into the space
        #Admissible values of x?
        if isinstance(x,(HarmonicCocycleElement,pAutomorphicFormElement)):
            return pAutomorphicFormElement(self,x)

    def _an_element_(self):
        r"""
        Returns an element of the module.
        """
        return pAutomorphicFormElement(self,1)

    def precision_cap(self):
        return self._prec

    def lift(self,f):
        r"""
        Lifts the harmonic cocycle ``f`` to an
        overconvegent automorphic form, thus computing
        all the moments.
        """
        F = self.element_class(self,f)
        F.improve()
        return F

    def _make_invariant(self, F):
        r"""
        EXAMPLES::
        """
        S = self._source.get_stabilizers()
        M  = [e.rep for e in self._list]
        newF = []
        for ii in range(len(S)):
            Si = S[ii]
            x = self._U(F[ii].moments) # was self._U(F[ii])
            if(any([v[2] for v in Si])):
                newFi = self._U(0)
                s = QQ(0)
                m = M[ii]
                for v in Si:
                    s += 1
                    # newFi  +=  x.r_act_by(m.adjoint()*self._source.embed(v[0],prec = self._prec)*m)  #TBC
                    newFi  +=  act_left( (m.adjoint() * self._source.embed(v[0],prec = self._prec)*m).adjoint(), x)
                newF.append((1/s)*newFi)
            else:
                newF.append(x)
        return newF

    def _apply_Up_operator(self,f,scale = False, fix_lowdeg_terms = True):
        r"""
        Apply the Up operator to ``f``.

        EXAMPLES::

            sage: X = BTQuotient(3,11)
            sage: M = HarmonicCocycles(X,4,30)
            sage: A = pAutomorphicForms(X,4,10, overconvergent = True)
            sage: F = A.lift(M.basis()[0]); F
            p-adic automorphic form on Space of automorphic forms on Quotient of the Bruhat Tits tree of GL_2(QQ_3) with discriminant 11 and level 1 with values in Overconvergent coefficient module of weight n = 2 over the ring 3-adic Field with capped relative precision 10 and depth 10:
            e   |   c(e)
            0 | 3^2 + O(3^12) + O(3^32)*z + O(3^26)*z^2 + (2*3^2 + 3^3 + 2*3^5 + 3^7 + 3^8 + 2*3^9 + O(3^10))*z^3 + (2*3^5 + 2*3^6 + 2*3^7 + 3^9 + O(3^10))*z^4 + (3^2 + 3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 3^7 + 2*3^8 + 3^9 + O(3^10))*z^5 + (3^2 + 2*3^3 + 3^4 + 2*3^6 + O(3^10))*z^6 + (2*3^3 + 3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 3^8 + 3^9 + O(3^10))*z^7 + (3^2 + 3^3 + 2*3^6 + 3^7 + 3^8 + 3^9 + O(3^10))*z^8 + (2*3^2 + 2*3^3 + 2*3^5 + 2*3^7 + 3^8 + 2*3^9 + O(3^10))*z^9
            1 | 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + 2*3^10 + 2*3^11 + O(3^12) + (3^2 + O(3^12))*z + (2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + 2*3^10 + 2*3^11 + O(3^12))*z^2 + (2*3^2 + 2*3^3 + 3^4 + 2*3^5 + 3^6 + 2*3^7 + 2*3^8 + O(3^10))*z^3 + (2*3^3 + 3^5 + 3^7 + 3^8 + O(3^10))*z^4 + (2*3^3 + 3^6 + 3^7 + 3^9 + O(3^10))*z^5 + (3^3 + 2*3^4 + 2*3^5 + 2*3^7 + 3^8 + 3^9 + O(3^10))*z^6 + (3^7 + 2*3^8 + 2*3^9 + O(3^10))*z^7 + (3^3 + 2*3^4 + 3^7 + O(3^10))*z^8 + (2*3^2 + 3^4 + 3^6 + 2*3^7 + 3^8 + 2*3^9 + O(3^10))*z^9
            2 | 3^2 + 2*3^3 + 2*3^6 + 3^7 + 2*3^8 + O(3^12) + (3 + 2*3^2 + 2*3^3 + 3^5 + 2*3^6 + 3^7 + 3^10 + O(3^11))*z + (2*3 + 2*3^2 + 3^4 + 2*3^5 + 2*3^6 + 2*3^8 + 3^10 + O(3^11))*z^2 + (2*3 + 3^2 + 2*3^7 + 3^9 + O(3^10))*z^3 + (3 + 2*3^2 + 2*3^4 + 3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + O(3^10))*z^4 + (3 + 3^2 + 3^4 + 2*3^9 + O(3^10))*z^5 + (3^3 + 2*3^5 + 3^6 + 3^8 + 2*3^9 + O(3^10))*z^6 + (3^5 + 2*3^7 + 3^9 + O(3^10))*z^7 + (2*3 + 3^3 + 3^4 + 2*3^6 + O(3^10))*z^8 + (2*3 + 2*3^3 + 2*3^4 + 2*3^6 + O(3^10))*z^9
            3 | 3^2 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + 2*3^10 + 2*3^11 + O(3^12) + (3^3 + 2*3^4 + 2*3^8 + O(3^13))*z + (3 + 3^2 + 3^3 + 2*3^4 + 2*3^5 + 3^6 + 3^7 + 2*3^8 + 2*3^9 + 2*3^10 + O(3^11))*z^2 + (3^2 + 2*3^3 + 3^4 + 3^7 + 3^8 + 2*3^9 + O(3^10))*z^3 + (3 + 2*3^2 + 2*3^3 + 3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 3^8 + O(3^10))*z^4 + (3 + 3^3 + 3^4 + 2*3^5 + 2*3^6 + 3^7 + 2*3^8 + 2*3^9 + O(3^10))*z^5 + (3 + 3^4 + 3^5 + 3^6 + 2*3^7 + O(3^10))*z^6 + (2*3 + 3^2 + 2*3^3 + 3^4 + 2*3^6 + 3^8 + 3^9 + O(3^10))*z^7 + (3 + 3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 3^7 + 2*3^9 + O(3^10))*z^8 + (2*3^2 + 3^4 + 3^5 + 3^8 + 3^9 + O(3^10))*z^9
        """
        HeckeData = self._source._get_Up_data()
        if scale == False:
            factor = (self._p)**(self._U.weight()/2)
        else:
            factor = 1
        Tf = []
        for jj in range(len(self._list)):
            tmp = self._U(0)
            for d in HeckeData:
                gg = d[0] # acter
                u = d[1][jj] # edge_list[jj]
                r = self._p**(-(u.power)) * (u.t(self.precision_cap()+1)*gg)
                # tmp += f._value[u.label].r_act_by(r) #TBC
                tmp += act_left(r.adjoint(), f._value[u.label]) #TBC
            tmp  *=  factor
            for ii in range(self._n+1):
                tmp.moments[ii] = f._value[jj].moments[ii]
            Tf.append(tmp)
        return pAutomorphicFormElement(self,Tf,quick = True)


