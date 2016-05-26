import sys, time, pdb, sys
import numpy as np
from scipy import linalg, stats, optimize, stats
from scipy.stats import norm,chi2, pearsonr
class SLM:
	def __init__(self,Y1,Y2d,X,f,verbose=False):
		self.Y1 = Y1 
		self.Y2d = Y2d
		self.X = X
		self.n = self.Y1.shape[0]
		self.f = f
		self.initialize()
		self.verbose = verbose
	def initialize(self):
		self.pvalue = 1
		self.beta1 = 0
		self.beta2 = 0
		self.t = norm.ppf(1-self.f)
		self.Y2 = np.copy(self.Y2d)
		for i in xrange(self.Y2d.shape[0]):
			if self.Y2d[i,0] == 1:
				self.Y2[i,0] = norm.pdf(self.t) / ( 1- norm.cdf(self.t))
			else:
				self.Y2[i,0] = - norm.pdf(self.t) / norm.cdf(self.t)
		self.N = np.dot(self.X.T,self.X)
		self.residue1 = self.Y1 - self.Y1.mean()
		self.residue2 = self.Y2 - self.Y2.mean()
		self.rho = float(np.dot(self.residue1.T,self.residue2) / np.sqrt(np.dot(self.residue1.T,self.residue1)*np.dot(self.residue2.T,self.residue2)))		
		self.sigma1 = np.sqrt(np.dot(self.residue1.T, self.residue1)/(self.n - 1))
		self.sigma2 = np.sqrt(np.dot(self.residue2.T,self.residue2)/(self.n-1))
		self.LL_control = -self.n * np.log(4*np.pi**2*self.sigma2*self.sigma1*(1-self.rho**2)) - 1/2/(1-self.rho**2)*(np.dot(self.residue1.T,self.residue1)/self.sigma1**2 - 2 * self.rho * np.dot(self.residue1.T,self.residue2) / self.sigma1 / self.sigma2 + np.dot(self.residue2.T,self.residue2) / self.sigma2 ** 2)			
		self.beta1 = float(np.dot(self.X.T,self.Y1) / self.N)
		self.residue1 = self.Y1 - self.Y1.mean() - self.X * self.beta1
		self.sigma1 = np.sqrt(float(np.dot(self.residue1.T,self.residue1) / (self.n-2)))
	def fit(self):
		if self.verbose:
			print 'rho\tbeta1\tbeta2\tsigma2'
		self.em()
		self.residue1 = self.Y1 - self.Y1.mean() - self.X * self.beta1 
		self.residue2 = self.Y2 - self.Y2.mean() - self.X * self.beta2
		self.sigma1 = np.sqrt(np.dot(self.residue1.T, self.residue1)/(self.n - 2))
		self.sigma2 = np.sqrt(np.dot(self.residue2.T,self.residue2)/(self.n-2))
		self.LL = -self.n * np.log(4*np.pi**2*self.sigma2*self.sigma1*(1-self.rho**2)) - 1/2/(1-self.rho**2)*(np.dot(self.residue1.T,self.residue1)/self.sigma1**2 - 2 * self.rho * np.dot(self.residue1.T,self.residue2) / self.sigma1 / self.sigma2 + np.dot(self.residue2.T,self.residue2) / self.sigma2 ** 2)
		self.pvalue = chi2.sf(2*(self.LL - self.LL_control),2)
		self.beta2_norm = float(self.beta2 / self.sigma2)
		self.pvalue = 1 - norm.cdf(self.beta2_norm * np.sqrt(self.n))
		if self.verbose:
			print self.beta2_norm, self.pvalue
	def em(self):
		beta2 = float(np.dot(self.X.T,self.Y2) / self.N)
		self.residue2 = self.Y2 - self.Y2.mean() - self.X * beta2
		self.sigma2 = np.sqrt(np.dot(self.residue2.T,self.residue2)/(self.n-2))
		rho = float(np.dot(self.residue1.T,self.residue2) / np.sqrt(np.dot(self.residue1.T,self.residue1)*np.dot(self.residue2.T,self.residue2)))
		if self.beta2 != 0 and abs((beta2 - self.beta2) / self.beta2) < 1e-3 and abs((rho - self.rho) / self.rho) < 1e-3:
			return self
		self.rho, self.beta2 = rho, beta2
		mu_star = self.Y2.mean() + self.X * self.beta2 + self.sigma2 / self.sigma1 * self.residue1 * self.rho
		sigma_star = np.sqrt((1-self.rho**2)*self.sigma2)
		alpha_star = (self.t - mu_star) / sigma_star
		if self.verbose:
			print self.rho, self.beta1, self.beta2, float(self.sigma2)
		for i in xrange(self.Y2.shape[0]):
			if self.Y2d[i,0] == 1:
				self.Y2[i,0] = mu_star[i] + sigma_star * norm.pdf(alpha_star[i])/max((1-norm.cdf(alpha_star[i])),1e-6)
			elif self.Y2d[i,0] == 0:
				self.Y2[i,0] = mu_star[i] - sigma_star * norm.pdf(alpha_star[i]) / max(norm.cdf(alpha_star[i]),1e-6)
		self.em()



