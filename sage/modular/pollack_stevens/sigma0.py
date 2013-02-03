r"""
The monoid `\Sigma_0(N)`.

This stands for a monoid of matrices over `\ZZ`, `\QQ`, `\ZZ_p`, or `\QQ_p`,
depending on a parameter `N`.

Over `\QQ` or `\ZZ`, it is the monoid of matrices of nonzero determinant if `N
= 0`, or if `N \ge 1` the monoid of matrices `[a, b; c, d]` with `a` a unit at
the primes dividing `N`, and `c` integral and congruent to 0 at the primes
dividing `N`.

Over `\Qp` or `\Zp`, the definition is the same, but we only allow `N = 0` or `N` a power of `p`.
"""

from sage.matrix.matrix_integer_2x2 import MatrixSpace_ZZ_2x2
from sage.matrix.matrix_space import MatrixSpace
from sage.structure.factory import UniqueFactory
from sage.structure.element import MonoidElement
from sage.categories.monoids import Monoids
from sage.categories.morphism import Morphism
from sage.structure.parent import Parent
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.infinity import Infinity
from sage.structure.unique_representation import UniqueRepresentation

# Need this to be pickleable
class _default_tuplegen(UniqueRepresentation):
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
            (1, 2, 3, 4)
        """
        return g[0,0], g[0,1], g[1,0], g[1,1]

class Sigma0_factory(UniqueFactory):

    def create_key(self, N, base_ring=ZZ, tuplegen=None):
        N = ZZ(N)
        if N <= 0:
            raise ValueError("Modulus should be > 0")
        if tuplegen is None:
            tuplegen = _default_tuplegen()

        if base_ring not in (QQ, ZZ):
            try:
                if not N.is_power_of(base_ring.prime()):
                    raise ValueError("Modulus must equal base ring prime")
            except AttributeError:
                raise ValueError("Base ring must be QQ, ZZ or a p-adic field")
        return (N, base_ring, tuplegen)

    def create_object(self, version, key):
        return Sigma0_class(*key)

Sigma0 = Sigma0_factory('sage.modular.pollack_stevens.sigma0.Sigma0')

class Sigma0Element(MonoidElement):

    def __init__(self, parent, mat):
        self._mat = mat
        MonoidElement.__init__(self, parent)

    def __hash__(self):
        return hash(self.matrix())

    def det(self):
        return self.matrix().det()

    def _mul_(self, other):
        return self.parent()(self._mat * other._mat, check=False)

    def __cmp__(self, other):
        return cmp(self._mat, other._mat)

    def _repr_(self):
        return self.matrix().__repr__()

    def matrix(self):
        return self._mat
    
    def __getitem__(self, *args):
        return self._mat.__getitem__(*args)

    def inverse(self):
        return self.parent()(self._mat._invert_unit())

class Sigma0Embedding(Morphism):
    r"""
    The embedding of `\Sigma_0` into the appropriate matrix space. This snippet
    of code is used so "x * y" will work if x is a matrix and y is a Sigma0
    element (and will return a matrix).
    """
    def __init__(self, domain):
        Morphism.__init__(self, domain, domain._matrix_space)

    def _call_(self, x):
        return x.matrix()

class Sigma0_class(Parent):

    Element = Sigma0Element

    def __init__(self, N, base_ring,tuplegen):
        self._N = N
        self._primes = list(N.factor())
        self._base_ring = base_ring
        self._tuplegen = tuplegen
        if base_ring == ZZ:
            self._matrix_space = MatrixSpace_ZZ_2x2()
        else:
            self._matrix_space = MatrixSpace(base_ring, 2)
        Parent.__init__(self, category=Monoids())
        self.register_embedding(Sigma0Embedding(self))

    def _an_element_(self):
        return self([1,0,0,1])

    def __cmp__(self, other):
        return cmp(type(self), type(other)) \
            or cmp(self._N, other._N)

    def level(self):
        return self._N

    def base_ring(self):
        return self._base_ring

    def _coerce_map_from_(self, other):
        r"""
        The *only* thing that coerces canonically into `\Sigma_0` is another
        `\Sigma_0`. It is *very bad* if integers are allowed to coerce in, as
        this leads to a noncommutative coercion diagram.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.sigma0 import Sigma0
            sage: Sigma0(1, QQ).has_coerce_map_from(Sigma0(3, ZZ))
            True
            sage: Sigma0(1, ZZ).has_coerce_map_from(ZZ)
            False
        """
        if isinstance(other, Sigma0_class) \
            and self.level().divides(other.level()) \
            and self.base_ring().has_coerce_map_from(other.base_ring()): 
            return True
        else:
            return False

    def _element_constructor_(self, x, check=True):
        if isinstance(x, Sigma0Element):
            x = x.matrix()
        if check:
            x = self._matrix_space(x)
            a,b,c,d = self._tuplegen(x)
            for (p, e) in self._primes:
                if c.valuation(p) < e:
                    raise ValueError("level %s^%s does not divide %s" % (p, e, c))
                if a.valuation(p) != 0:
                    raise ValueError("%s is not a unit at %s" % (a, p))
            if x.det() == 0:
                raise ValueError("matrix must be nonsingular")
        x.set_immutable()
        return self.element_class(self, x)

    def _repr_(self):
        return 'Monoid Sigma0(%s) with coefficients in %s' % (self.level(), self.base_ring())
