from sage.structure.sage_object import SageObject

"""A family of distributions -- i.e. an element of R \hat{\otimes} D -- of the form \sum_i f_i(w) \otimes \mu_i is represented by a vector whose i-th entry is f_i(w) --- here \mu_i is the distribution which sends z^i to 1 and all other z^j's to 0."""

class dist_fam(SageObject):
	def __init__(self,p,deg,moments):
		"""deg denotes the max degree in the weight variable"""
		self.p=p
		self.deg=deg
		self.moments=moments
			
	def __repr__(self):
		return repr(self.moments)

	def moment(self,n):
		return self.moments[n]

	def num_moments(self):
		return len(self.moments)

	def change_deg(self,new_deg):
		assert new_deg<=self.deg, "can only lower degree"
		v=[self.moments[a].truncate(new_deg) for a in range(self.num_moments())]
		return dist_fam(self.p,new_deg,vector(v))

	def truncate(self):
		v=self.moments
		w=[v[j].truncate(self.deg) for j in range(self.num_moments())]
		return dist_fam(self.p,self.deg,vector(w))

	def __add__(self,right):
		assert self.num_moments()==right.num_moments(), "the accuracies are different"
		assert self.deg==right.deg, "the degrees in w are different"
		return dist_fam(self.p,self.deg,self.moments+right.moments)

	def scale(self,left):
		return dist_fam(self.p,self.deg,(left*self.moments)).truncate()

	def __sub__(self,right):
		return self+right.scale(-1)

	def __cmp__(self,right):
		return cmp((self.p,self.deg,self.moments),(right.p,right.deg,right.moments))

	def zero(self):
		return dist_fam(self.p,self.deg,vector([self.moment(0)*0 for i in range(0,len(self.moments))]))

	def gen(self):
		return self.moment(0).parent().gen()

	def specialize(self,k):
		"""evaluates at ((1+p)^k-1)/p"""
		w=self.gen()
		v=[]
		for j in range(0,self.num_moments()):
			v=v+[Rational(self.moment(j).substitute(w=((1+self.p)^k-1)/self.p))]
		return dist(self.p,k,vector(v))

	def valuation(self):
		return min([val(self.moment(j),self.p) for j in range(self.num_moments())])

	def normalize(self):
		N=self.num_moments()
		v=[]
		for j in range(0,N):
			v=v+[normalize(self.moment(j),self.p,j,N)]
		return dist_fam(self.p,self.deg,vector(v))

	def change_precision(self,M):
		"""only hangs onto M moments"""
		assert M<=self.num_moments(),"not enough moments"

		v=[self.moment(i) for i in range(M)]
		mu=dist_fam(self.p,self.deg,vector(v))
		return mu

	def act_by_ps_fam(self,F):
		gam=form_acting_matrix_on_dist_fam(F)
		v=(Matrix(self.moments)*gam)[0]
		w=[v[j].truncate(self.deg) for j in range(self.num_moments())]
		return dist_fam(self.p,self.deg,vector(w))

	def act_right_weight_zero(self,gam):
		a=gam[0,0]
		b=gam[0,1]
		c=gam[1,0]
		d=gam[1,1]
		G=form_acting_matrix_on_dist(self.p,self.num_moments(),0,a,b,c,d)

		return dist_fam(self.p,self.deg,(Matrix(self.moments)*G)[0])

	def act_right(self,gam):
		w=self.gen()
		K=aut(self.p,self.deg,self.num_moments(),gam[0,0],gam[1,0],w)
		return self.act_by_ps_fam(K).act_right_weight_zero(gam)
		
	def solve_diff_eqn(self):
		w=self.gen()
		mus=self.zero()
		for j in range(1,self.num_moments()):		
			v=[Rational(0) for i in range(self.num_moments())]
			v[j]=Rational(1)
			mu=dist(self.p,0,v)
			nu=mu.solve_diff_eqn()
			mus=mus+nu.lift_to_dist_fam(self.deg,w).scale(self.moment(j))
		return mus

@cached_function
def form_acting_matrix_on_dist_fam(F):
	"""first row is just F, then F shifted over 1, etc."""
	list=copy(F.series)
	v=[]
	for j in range(0,len(list)):
		v=v+[copy(list)]
		list.insert(0,0)
		list.pop()
	return Matrix(v).transpose()

def normalize(F,p,r,N):
	v=F.list()
	M=floor((N+1-r)*(p-2)/(p-1))
	v=[v[a]%(p^M) for a in range(len(v))]
	S=F.parent()
	return S(v)

def val(F,p):
	v=F.list()
	if v==[]:
		return Infinity
	else:
		return min([valuation(v[j],p) for j in range(len(v))])
