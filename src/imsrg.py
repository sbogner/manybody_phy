import numpy as np
from numpy import array, dot, diag, reshape, transpose
from scipy.linalg import eigvalsh
from scipy.integrate import odeint, ode
from sys import argv

class IMSRG:

	def __init__(self, d, g, smax, ds, holes, particles):

		self.d = d
		self.g = g
		self.smax = smax
		self.ds = ds
		self.holes = holes
		self.parts = particles
		self.dim1B = len(holes)+len(particles)
		self.dim2B = self.dim1B**2

		# bases and indices
		self.bas1B = [i for i in range(self.dim1B)]
		self.build_bas2B()
		self.build_ph_bas2B()

		# occupation number matrices
		self.build_occ1B()
		self.build_occ2B()
		self.build_ph_occ2B_A()

		# pairing hamiltonian
		self.normal_order()

		self.calc_eta_wegner()



	def build_bas2B(self):

		# construct list of states
		self.bas2B = []

		for i in self.holes:
			for j in self.holes:
				self.bas2B.append((i,j))

		for i in self.holes:
			for a in self.parts:
				self.bas2B.append((i,a))

		for a in self.parts:
			for i in self.holes:
				self.bas2B.append((a,i))

		for a in self.parts:
			for b in self.parts:
				self.bas2B.append((a,b))


		# store indices of states in dictionary
		self.idx2B = {}

		for i, state in enumerate(self.bas2B):
			self.idx2B[state] = i;


	def build_ph_bas2B(self):

		# construct list of states
		self.ph_bas2B = []

		for i in self.holes:
			for j in self.holes:
				self.ph_bas2B.append((i,j))

		for i in self.holes:
			for a in self.parts:
				self.ph_bas2B.append((i,a))

		for a in self.parts:
			for i in self.holes:
				self.ph_bas2B.append((a,i))

		for a in self.parts:
			for b in self.parts:
				self.ph_bas2B.append((a,b))


		# store indices of states in dictionary
		self.ph_idx2B = {}

		for i, state in enumerate(self.ph_bas2B):
			self.ph_idx2B[state] = i;


	def build_occ1B(self):

		self.occ1B = np.zeros(self.dim1B)

		for i in self.holes:
			self.occ1B[i] = 1


	def build_occ2B(self):

		# n_p - n_q
		self.occ2B_A = np.zeros((self.dim2B,self.dim2B))

		# 1 - n_p - n_q
		self.occ2B_B = np.zeros((self.dim2B,self.dim2B))

		# n_p * n_q
		self.occ2B_C = np.zeros((self.dim2B,self.dim2B))


		for i, (p,q) in enumerate(self.bas2B):
			self.occ2B_A[i,i] = self.occ1B[p]-self.occ1B[q]
			self.occ2B_B[i,i] = 1-self.occ1B[p]-self.occ1B[q]
			self.occ2B_C[i,i] = self.occ1B[p]*self.occ1B[q]


	def build_ph_occ2B_A(self):

		# n_p - n_q
		self.ph_occ2B_A = np.zeros((self.dim2B,self.dim2B))

		for i, (p,q) in enumerate(self.ph_bas2B):
			self.ph_occ2B_A[i,i] = self.occ1B[p]-self.occ1B[q]


	def normal_order(self):

		# construct pairing hamiltonian
		H1B = np.zeros((self.dim1B,self.dim1B))
		H2B = np.zeros((self.dim2B,self.dim2B))

		for i in self.bas1B:
			H1B[i,i] = self.d*np.floor_divide(i,2)

		for (i,j) in self.bas2B:
			if (i%2==0 and j==i+1):
				for (k,l) in self.bas2B:
					if (k%2==0 and l==k+1):
						H2B[self.idx2B[(i,j)],self.idx2B[(k,l)]] = -0.5*self.g
						H2B[self.idx2B[(i,j)],self.idx2B[(l,k)]] = 0.5*self.g
						H2B[self.idx2B[(j,i)],self.idx2B[(k,l)]] = 0.5*self.g
						H2B[self.idx2B[(j,i)],self.idx2B[(l,k)]] = -0.5*self.g


		# normal order hamiltonian
		self.E = 0.0
		for i in self.holes:
			self.E += H1B[i,i]
			for j in self.holes:
				self.E += 0.5*H2B[self.idx2B[(i,j)],self.idx2B[(i,j)]]

		self.f = H1B
		for p in self.bas1B:
			for q in self.bas1B:
				for i in self.holes:
					self.f[p,q] += H2B[self.idx2B[(p,i)],self.idx2B[(q,i)]]

		self.Gamma = H2B


	def ph_transform2B(self, matrix2B):

		ph_matrix2B = np.zeros(self.dim2B,self.dim2B)

		for i, (p,q) in enumerate(self.ph_bas2B):
			for j, (r,s) in enumerate(self.ph_bas2B):
				ph_matrix2B[i,j] -= matrix2B[idx2B[(p,s)],idx2B[(r,q)]]

		return ph_matrix2B


	def inverse_ph_transform2B(self, ph_matrix2B):

		matrix2B = np.zeros(self.dim2B,self.dim2B)

		for i, (p,q) in enumerate(self.ph_bas2B):
			for j, (r,s) in enumerate(self.ph_bas2B):
				matrix2B[i,j] -= ph_matrix2B[idx2B[(p,s)],idx2B[(r,q)]]

		return matrix2B


	def commutator(self, A, B):

		return dot(A,B)-dot(B,A)


	def commutator2B(self, A1B, A2B, B1B, B2B):

		# zero-body part
		C0B = 0.0

		# 1B-1B
		for i in self.holes:
			for a in self.parts:
				C0B += (A1B[i,a]*B1B[a,i]-A1B[a,i]*B1B[i,a])

		# 2B-2B
		for i in self.holes:
			for j in self.holes:
				for a in self.parts:
					for b in self.parts:
						ij = self.idx2B[(i,j)]
						ab = self.idx2B[(a,b)]
						C0B += 0.25*(A2B[ij,ab]*B2B[ab,ij]-B2B[ij,ab]*A2B[ab,ij])


		# one-body part
		# 1B-1B
		C1B = self.commutator(A1B,B1B)

		# 1B-2B
		for p in range(self.dim1B):
			for q in range(self.dim1B):
				for i in self.holes:
					for a in self.parts:
						ap = self.idx2B[(a,p)]
						iq = self.idx2B[(i,q)]
						ip = self.idx2B[(i,p)]
						aq = self.idx2B[(a,q)]
						C1B[p,q] += (A1B[i,a]*B2B[ap,iq]-B1B[i,a]*A2B[ap,iq]+B[a,i]*A2B[ip,aq]-A[a,i]*B2B[ip,aq])

		# 2B-2B
		AB = dot(A2B,dot(self.occ2B_B,B2B))
		ABT = transpose(AB)
		for p in range(self.dim1B):
			for q in range(self.dim1B):
				for i in self.holes:
					ip = self.idx2B[(i,p)]
					iq = self.idx2B[(i,q)]
					C1B[p,q] += 0.5*(AB[ip,iq]+ABT[ip,iq])

		AB = dot(A2B,dot(self.occ2B_C,B2B))
		ABT = transpose(AB)
		for p in range(self.dim1B):
			for q in range(self.dim1B):
				for r in range(self.dim1B):
					rp = self.idx2B[(r,p)]
					rq = self.idx2B[(r,q)]
					C1B[p,q] += 0.5*(AB[rp,rq]+ABT[rp,rq])

		# two-body part
		C2B = np.zeros_like((self.dim2B,self.dim2B))

		# 1B-2B
		for p in range(self.dim1B):
			for q in range(self.dim1B):
				for r in range(self.dim1B):
					for s in range(self.dim1B):
						pq = self.idx2B[(p,q)]
						rs = self.idx2B[(r,s)]
						for t in range(self.dim1B):
							tq = self.idx2B[(t,q)]
							tp = self.idx2B[(t,p)]
							ts = self.idx2B[(t,s)]
							tr = self.idx2B[(t,r)]
							C2B[pq,rs] += 0.25*(A1B[p,t]*B2B[tq,rs]-B1B[p,t]*A2B[tq,rs]
								               -A1B[q,t]*B2B[tp,rs]+B1B[q,t]*A2B[tp,rs]
								               -A1B[t,r]*B2B[pq,ts]+B1B[t,r]*A2B[pq,ts]
								               +A1B[t,s]*B2B[pq,tr]-B1B[t,s]*A2B[pq,tr])

		# 2B-2B
		AB = dot(A2B,dot(self.occ2B_B,B2B))
		ABT = transpose(AB)
		C2B += 0.125*(AB+ABT)
		


		#return C0B, C1B, C2B


	def calc_eta_wegner(self):

		# split into diag and off-diag parts
		fod = np.zeros_like(self.f)
		for i in self.holes:
			for a in self.parts:
				fod[i,a] = self.f[i,a]
				fod[a,i] = self.f[a,i]
		fd = self.f-fod

		Gammaod = np.zeros_like(self.Gamma)
		for i in self.holes:
			for j in self.holes:
				for a in self.parts:
					for b in self.parts:
						ij = self.idx2B[(i,j)]
						ab = self.idx2B[(a,b)]
						Gammaod[ij,ab] = self.Gamma[ij,ab]
						Gammaod[ab,ij] = self.Gamma[ab,ij]
		Gammad = self.Gamma-Gammaod

		self.commutator2B(fd,Gammad,fod,Gammaod)



holes = [0,1,2,3]
particles = [4,5,6,7]

PairingModel = IMSRG(1.0,0.5,10.0,0.01,holes,particles)


