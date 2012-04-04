class modsym_dist(modsym):
	def ms(self):
		"""demotes to a regular modular symbol"""
		return modsym(self.level,self.data,self.manin)

	def p(self):
		"""returns the underlying prime"""
		return self.data[0].p

	def weight(self):
		"""returns the underlying weight"""		
		return self.data[0].weight

	def num_moments(self):
		"""returns the number of moments of each value of self"""
		return self.data[0].num_moments()

	def eval(self,A):
		"""here A is a 2x2 matrix and this function returns self evaluated and the divisor attached to A = A(\infty)-A(0)"""
		ans=self.ms().eval(A)
		return ans.normalize()

	def specialize(self):
		"""returns the underlying classical symbol of weight k -- i.e. applies the canonical map D_k --> Sym^k to all values of self"""
		v=[]
		for j in range(0,len(self.data)):
			v=v+[self.data[j].specialize()]
		return modsym_symk(self.level,v,self.manin)
	
	def valuation(self):
		"""returns the exponent of the highest power of p which divides all moments of all values of self"""
		return min([self.data[j].valuation() for j in range(0,len(self.data))])

	def normalize(self):
		"""normalized every values of self -- e.g. reduces each values j-th moment modulo p^(N-j)"""
		assert self.valuation()>=0, "moments are not integral"
		v=[]
		for j in range(0,len(self.data)):
			v=v+[self.data[j].normalize()]
		ans=modsym_dist(self.level,v,self.manin,self.full_data)
		ans.normalize_full_data()
		
		return ans

	def change_precision(self,M):
		"""only holds on to M moments of each value of self"""
		v=[]
		for j in range(0,len(self.data)):
			v=v+[self.data[j].change_precision(M)]
		return modsym_dist(self.level,v,self.manin)

	def lift(self,phi,ap):
		"""Greenberg trick of lifting and applying U_p --- phi is the (exact) classical symbol"""
		v=[]
		for a in range(self.ngens()):
			v=v+[self.data[a].lift()]
		Phi=modsym_dist(self.level,v,self.manin)
		k=self.weight()
		for a in range(self.ngens()):
			for j in range(k+1):
				Phi.data[a].moments[j]=(phi.data[a].coef(j))%(p^(Phi.num_moments()))
		return Phi.hecke(self.p()).scale(1/ap).normalize()

	def lift_to_modsym_dist_fam(self,w):
		p=self.p()
		deg=floor((self.num_moments()+1)*(p-2)/(p-1))
		v=[]
		for j in range(1,len(self.manin.gens)):
			v=v+[self.data[j].lift_to_dist_fam(deg,w)]
		#does this deal with two-torsion?
		t=v[0].zero()
		list=self.grab_relations()
		for j in range(2,len(list)):
			R=list[j]
			index=R[0][2]
			rj=self.manin.gens.index(index)
			t=t+self.data[rj].lift_to_dist_fam(deg,w).act_right(R[0][1]).scale(R[0][0])
		s=2
		gam=list[s][0][1]
		while (gam[0,0],gam[1,0]) == (1,0):
			s=s+1
			gam=list[s][0][1]
		K1=aut(self.data[0].p,deg,self.data[0].num_moments(),gam[0,0],gam[1,0],w).series[1]
		t0=t.moment(0)
		S=t0.parent()
		R=PowerSeriesRing(QQ,'ww')
		ww=R.gen()
		K1=K1.substitute({w:ww})	
		t0=t0.substitute({w:ww})	
		err=self.zero_elt().lift_to_dist_fam(deg,w)
		err.moments[1]=((t0/(K1)).truncate(deg)).substitute({ww:w})
		alt=(t0/K1).truncate(deg)
		print "The result will be a lifting modulo p^",valuation(alt.constant_coefficient(),p)
		v[0]=v[0]+err
		t=t+err-err.act_right(gam)
		#print "***",t.moment(0)
		mu=t.solve_diff_eqn()
		v=[mu.scale(-1)]+v

		Phis=modsym_dist_fam(self.level,v,self.manin)

		return Phis

	def is_Tq_eigen(Phi,q):
		Phiq=Phi.hecke(q)
		c=Phiq.data[0].moment(0)/Phi.data[0].moment(0)
		M=Phi.data[0].num_moments()
		p=Phi.data[0].p
		print Phiq-Phi.scale(c)

		return c%(p^M)
		
	
def random_OMS(N,p,k,M):
	"""Returns a random OMS with tame level N, prime p, weight k, and M moments --- requires no 2 or 3-torsion"""
	manin=manin_relations(N*p)
	v=[]
	for j in range(1,len(manin.gens)):
		g=manin.gens[j]
		if manin.twotor.count(g)==0:
			v=v+[random_dist(p,k,M)]
		else:
			rj=manin.twotor.index(g)
			gam=manin.twotorrels[rj]
			mu=random_dist(p,k,M)
			v=v+[(mu.act_right(gam)-mu).scale(1/2)]
	t=v[0].zero()
	for j in range(2,len(manin.rels)):
		R=manin.rels[j]
		if len(R)==1:
			if R[0][0]==1:
				rj=manin.gens.index(j)
				t=t+v[rj-1]
			else:
				index=R[0][2]
				rj=manin.gens.index(index)
				mu=v[rj-1]
				t=t+mu.act_right(R[0][1]).scale(R[0][0])
	v=[mu.scale(-1)]+v
	Phi=modsym_dist(N*p,v,manin)	
	return Phi

def random_OMS_char(N,p,k,chi,M):
	"""Returns a random OMS with tame level N, prime p, weight k, and M moments --- requires no 2 or 3-torsion"""
	manin=manin_relations(N*p)
	v=[]
	for j in range(1,len(manin.gens)):
		g=manin.gens[j]
		print g
		print manin.glue[g][1]
		gam=manin.glue[g][1]
		a=gam[0,0]
		c=gam[1,0]
		if manin.twotor.count(g)==0:
			mu=random_dist_char(p,k,chi,M).zero()
			mu.moments[0]=c
			mu.moments[1]=(chi(a)-a)*a/chi(a)  #formula to force difference equation to work 
			v=v+[mu]
		else:
			v=v+[random_dist_char(p,k,chi,M).zero()]
#			rj=manin.twotor.index(g)
#			gam=manin.twotorrels[rj]
#			mu=random_dist_char(p,k,chi,M)
#			v=v+[(mu.act_right(gam)-mu).scale(1/2)]
	t=v[0].zero()
	for j in range(2,len(manin.rels)):
		R=manin.rels[j]
		if len(R)==1:
			if R[0][0]==1:
				rj=manin.gens.index(j)
				t=t+v[rj-1]
				print v[rj-1]
			else:
				index=R[0][2]
				rj=manin.gens.index(index)
				mu=v[rj-1]
				print mu.act_right(R[0][1]).scale(R[0][0]).normalize()
				t=t+mu.act_right(R[0][1]).scale(R[0][0])
	t=t.normalize()
	mu=t.solve_diff_eqn()
	v=[mu.scale(-1)]+v

	return modsym_dist(N*p,v,manin)	

