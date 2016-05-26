import os,sys
import numpy as np
import struct
import pdb
class plink:
	def __init__(self,fbase,phenoFile=None,f=0.5,normGenotype=True):
		self.fbase = fbase
		self.normGenotype = normGenotype
		self.indivs = self.getIndivs(self.fbase)
		self.f = f
	def __del__(self):
		if self.fhandle:
			self.fhandle.close()
		if self.snpFileHandle:
			self.snpFileHandle.close()
	def getSNPiterator_tped(self):
		file = self.fbase + '.map'
		self.numSNPs = 0
		with open(file) as f:
			for line in f:
				self.numSNPs += 1
		self.have_read = 0
		self.snpFileHandle = open(file,'r')
		file = self.fbase + '.tped'
		self.fhandle = open(file,'r')
		return self
	def __iter__(self):
		return self.getSNPiterator_tped()
	def next(self):
		if self.have_read == self.numSNPs:
			raise StopIteration
		self.have_read += 1
		X = self.fhandle.readline()
		XX = X.strip().split()
		chrm, rsid, pos1, pos2 = tuple(XX[:4])
		XX = XX[4:]
		G = self.getGenos_tped(XX)
		if self.normGenotype:
			G = self.normalizeGenotype(G)
		return G,self.snpFileHandle.readline().strip().split()[1]
	def getGenos_tped(self,X):
		G = []
		for i in range(0,len(X)-1,2):
			a = X[i]
			b = X[i+1]
			if a == b == '0':
				g = np.nan
			if a == b == '1':
				g = 0 
			if a == b == '2':
				g = 1
			if a != b:
				g = 0.5
			try:
				G.append(g)
			except:
				print a,b
				raw_input('pause')
		return np.array(G)
	def normalizeGenotype(self,G):
		G = G[self.keep]
		x = True - np.isnan(G)
		if not len(G[x]):
			return G[x]
		m = G[x].mean()
		if G[x].var == 0:
			s = 1.0
		else:
			s = np.sqrt(G[x].var())
		G[np.isnan(G)] = m
		if s == 0:
			G = G - m
		else:
			G = (G -m) / s
		return G
	def getIndivs(self,base):
		famFile = '%s.tfam' % base 
		keys = []
		i = 0
		f = open(famFile,'r')
		for line in f:
			v = line.strip().split()
			famId = v[0]
			indivId = v[1]
			k = (famId.strip(),indivId.strip())
			keys.append(k)
			i += 1
		f.close()
		self.N = len(keys)
		return keys
	def getPhenos(self,phenoFile=None):
		f = open(phenoFile,'r')
		keys = []
		P = []
		head = f.readline().strip().split()
		for line in f:
			v = line.strip().split()
			keys.append((v[0],v[1]))
			P.append([(x.strip() == 'NA') or x.strip() == '-9' and np.nan or float(x) for x in v[2:]])
		f.close()
		P = np.array(P)
		D = {}
		L = []
		keep = []
		for i in range(len(keys)):
			D[keys[i]] = i
		for i in range(len(self.indivs)):
			if not D.has_key(self.indivs[i]):
				continue
			keep.append(i)
			L.append(D[self.indivs[i]])
		P = P[L,:]
		self.phenos = P
		Y2d = self.phenos[:,1]
		self.phenos[:,0] = self.phenos[:,0]
		self.keep = keep
		return P


		