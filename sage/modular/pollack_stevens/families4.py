class WeightKAction_fam(WeightKAction_vector):
    """
    """
    def __init__(self, Dk, character, tuplegen, on_left):
        #Only difference is that it adds a cache for automorphy factors.
        self._autfactors = {}
        WeightKAction_vector.__init___(self, Dk, character, tuplegen, on_left)
    
    def clear_cache(self):
        #Only difference is that it clears the cache for automorphy factors.
        self._actmat = {}
        self._maxprecs = {}
        self._autfactors = {}
    
    def _compute_aut_factor_matrix(self, g):
        #compute the power series
        D = self.underlying_set()
        p_prec, var_prec = D.precision_cap()
        return automorphy_factor_matrix(D.prime(), g[0, 0], g[1,0], self._k, self._character, p_prec, var_prec, D.base_ring())
    
    def aut_factor_matrix(self, g):
        pass
        #manage cache
    
    def _call_(self, v, g):
        return v.parent()(v.moments * self.aut_factor_matrix(g) * self.acting_matrix(g, len(v.moments)))
