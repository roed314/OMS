class modsym_dist_fam(modsym):
	def ms(self):
		"""demotes to a regular modular symbol"""
		return modsym(self.level,self.data,self.manin)

	def num_moments(self):
		return self.data[0].num_moments()

	def p(self):
		return self.data[0].p

	def deg(self):
		return self.data[0].deg

	def change_deg(self,new_deg):
		v=[self.data[r].change_deg(new_deg) for r in range(self.ngens())]
		return modsym_dist_fam(self.level,v,self.manin)

	def specialize(self,k):
		v=[]
		for j in range(0,len(self.data)):
			v=v+[self.data[j].specialize(k)]
		return modsym_dist_aws(self.level,v,self.manin)
	
	def valuation(self):
		#print [self.data[j].valuation() for j in range(len(self.data))]
		return min([self.data[j].valuation() for j in range(len(self.data))])

	def normalize(self):
		v=[]
		for j in range(0,len(self.data)):
			v=v+[self.data[j].normalize()]
		return modsym_dist_fam(self.level,v,self.manin)
	
#	def normalize_aws(self):
#		v=[]
#		for j in range(0,len(self.data)):
#			v=v+[self.data[j].normalize_aws()]
#		return modsym_dist_fam(self.level,v,self.manin)

#	def normalize_aws2(self):
#		v=[]
#		for j in range(0,len(self.data)):
#			v=v+[self.data[j].normalize_aws2()]
#		return modsym_dist_fam(self.level,v,self.manin)

	def change_precision(self,M):
		v=[]
		for j in range(0,len(self.data)):
			v=v+[self.data[j].change_precision(M)]
		return modsym_dist_fam(self.level,v,self.manin)


		
	
