from sage.all import *

class pAdicLseries(SageObject):
    r"""
    The `p`-adic `L`-series associated to an overconvergent eigensymbol.
    """

    def __init__(self, symb):
        r"""

        INPUT:
            - ``symb`` -- overconvergent eigensymbol
        """
        self._symb = symb

    def symb(self):
        r"""
        """
        return self._symb

    def prime(self):
        r"""
        """
        return self._symb.parent().prime()

    if self.prime() == none:
        raise ValueError ("Not a p-adic overconvergent modular symbol.")

    def _repr_(self):
        r"""
        Return print representation.
        """
        s = "%s-adic L-series of $s"%(self.prime(), self.symb())
        return s

    def series(self, n, quadratic_twist, prec):
        r"""
        """
        pass

    def twisted_symbol_on_Da(self, a, D): # rename! should this be in modsym?
        """
        Returns `\Phi_{\chi}(\{a/p}-{\infty})` where `Phi` is the OMS
        corresponding to self and `\chi` is a character of conductor `D`


        INPUT:
            - ``a`` -- integer in [0..p-1]
            - ``D`` -- conductor of the quadratic twist `\chi`

        OUTPUT:

        `\Phi_{\chi}(\{a/p\}-\{\infty\})`

        EXAMPLES:

        """
        symb = self.symb()
        p = symb.parent().prime()
        twisted_dist = symb.parent().zero_element()
        m_map = symb._map
        for b in range(1, abs(D) + 1):
            if gcd(b, D) == 1:
                M1 = M2Z([1, b / abs(D), 0, 1])
                new_dist = m_map.__call__(M1 * M2Z[a, 1, p, 0])._act_right(M1)
                new_dist = new_dist.scale(kronecker(D, b)).normalize()
                twisted_dist = twisted_dist._add(new_dist)
                #ans = ans + self.eval(M1 * M2Z[a, 1, p, 0])._right_action(M1)._lmul_(kronecker(D, b)).normalize()
        return PSModularSybmolElement(twisted_dist.normalize(), symb.parent())

    def basic_integral(self, a, j, D=None):
        r"""
        Returns `\int_{a+pZ_p} (z-{a})^j d\Phi(0-infty)`
        -- see formula [Pollack-Stevens, sec 9.2]

        """
        if D == None:
            D = 1
        #check that a is between 0 and p - 1
        symb = self.symb()
        M = symb.precision_cap()
        if j > M:
            raise PrecisionError ("Too many moments for %s."%s(symb))
        p = self.prime()
        ap = symb.ap(p)
        ap = ap * kronecker(D, p)
        K = Qp(p, M)
        symb_twisted = twisted_symbol_on_Da(symb, a, D)
        return sum(binomial(j, r) * ((a - ZZ(K.teichmuller(a)))**(j - r)) *
                (p**r) * self.phi_on_Da(a, D).moment(r) for r in range(j+1)) / ap

def log_gamma_binomial(p,gamma,z,n,M):
    r"""
    Returns the list of coefficients in the power series
    expansion (up to precision `M`) of `{\log_p(z)/\log_p(\gamma) \choose n}`

    INPUT:

        - ``p`` --  prime
        - ``gam`` -- topological generator e.g., `1+p`
        - ``z`` -- variable
        - ``n`` -- nonnegative integer
        - ``M`` -- precision

    OUTPUT:

    The list of coefficients in the power series expansion of
    `{\log_p(z)/\log_p(\gamma) \choose n}`

    EXAMPLES:

        sage: R.<z> = QQ['z']
        sage: loggam_binom(5,1+5,z,2,4)
        [0, -3/205, 651/84050, -223/42025]
        sage: loggam_binom(5,1+5,z,3,4)
        [0, 2/205, -223/42025, 95228/25845375]
    """
    L = sum([ZZ(-1)**j/j*(gamma-1)**j for j in range (1,M)]) #log_p(1+z)
    loggam = L/(L(gamma-1))                  #log_{gamma}(1+z)= log_p(1+z)/log_p(gamma)
    return z.parent()(binomial(loggam,n)).truncate(M).list()
