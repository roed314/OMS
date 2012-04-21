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

    def basic_integral(self,a,j,ap,D=None):
        r"""
        Returns `\int_{a+pZ_p} (z-{a})^j d\Phi(0-infty)`
        -- see formula [Pollack-Stevens, sec 9.2]

        """
        symb = self.symb()
        M = symb.num_moments()
        p = self.prime()
        K = Qp(p,M)
        ap = ap*kronecker(D,p)
        return sum(binomial(j,r)*((a- ZZ(K.teichmuller(a)))**(j-r))*(p**r)*self.phi_on_Da(symb,a,D).moment(r) for r in range(j+1))/ap


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
