import numpy as np 
import scipy

import Assignment1 #utilize the functions from last assignement
import numpy as np
import matplotlib.pyplot as plt

from datetime import datetime
from dateutil.parser import parse

import pandas as pd
import warnings

class adjustedMatrix:
	def __init__(self, Nt,Ns):
		self.Ns = Ns
		self.Nt = Nt
		self.matrix = np.zeros( (Nt + 1, (Ns * 2) + 1))

	def __getitem__(self, tup):
		try:
			t,s = tup
			return self.matrix[t,s + self.Ns]
		except:
			return self.matrix[tup, :]
	def __setitem__(self, tup, value):
		try:
			t,s = tup
			self.matrix[t,s+self.Ns] = value
		except:
			self.matrix[tup,:] = value
	def __str__(self):
		return self.matrix


def solveTriDiagonal(a,b,c,ans):
	for i in range(len(ans)):
		if i == 0:
			pass
		else:
			w = a[i] / b[i-1]
			b[i] = b[i] - (w * c[i-1])
			ans[i] = ans[i] - w * ans[i-1]

	x = np.zeros(len(ans))
	x[len(ans) -1] = ans[len(ans) -1] / b[len(ans) -1]
	for i in range(len(ans)-2, -1, -1):
		x[i] = (ans[i] - c[i] * x[i+1]) / b[i]

	return x

def ExplicitFiniteDifference(Nt,Np,dp,S,K,sigma,div, r,T, call = True, getGreek = False):
	'''
	Nt - Number of time steps to discritize
	Np - Number of prices to discritize on each size of the current point
	dp - Change in price per discrete price point
	S - Current Price
	K - Strike
	sigma - volitility
	div - dividend
	'''
	dt = T/N
	u = np.exp(sigma * np.sqrt(2*dt) ) 
	d = np.exp(-sigma * np.sqrt(2*dt) )
	m = 1

	pu = ( (np.exp( (r-div)*dt/2) - np.exp(-sigma * np.sqrt(dt/2) ) ) / (np.exp(sigma*np.sqrt(dt/2)) - np.exp(-sigma*np.sqrt(dt/2))) )**2
	pd = ( (np.exp(sigma* np.sqrt(dt/2)) - np.exp( (r-div)*dt/2 ) )/ (np.exp(sigma*np.sqrt(dt/2)) - np.exp( -sigma*np.sqrt(dt/2))) ) **2
	pm = 1 - (pu + pd)
	tTree = np.zeros((N+1, ((N+1) * 2)-1 ))
	#Create Tree of Prices
	tTree[0][0] = S
	for i in range(1,N+1):
		tTree[i][0] = tTree[i-1][0] * u
		for j in range(1, (i*2)):
			tTree[i][j] = tTree[i-1][j-1]
		tTree[i][(i*2)] = tTree[i-1][(i - 1)*2] * d
	
	optionTree = np.zeros((N+1, ((N+1) * 2)-1 ))
	#Fill in known values at expiration
	for i in range(N+1):
		if call:
			optionTree[N,i] = max(0, tTree[N][i]-K)
		else:
			optionTree[N,i] = max(0, K-tTree[N][i])

	#back solve values
	for i in range(N-1, -1, -1):
		for j in range( (i*2)+1 ):
			if call:
				optionTree[i][j] = max(0,tTree[i,j]-K, np.exp(-r*dt)*(pu*optionTree[i+1][j] + pm*optionTree[i+1][j+1] + pd*optionTree[i+1][j+2]) )
			else:
				optionTree[i][j] = max(0, K-tTree[i,j], np.exp(-r*dt)*(pu*optionTree[i+1][j] + pm*optionTree[i+1][j+1] + pd*optionTree[i+1][j+2]) )

	return optionTree[0][0]

def ImplicitFiniteDifference(Nt,Np,dp,S,K,sigma,div, r,T, call = True):
	'''
	Nt - Number of time steps to discritize
	Np - Number of prices to discritize on each size of the current point
	dp - Change in price per discrete price point
	S - Current Price
	K - Strike
	sigma - volitility
	div - dividend
	'''
	dt = T/Np
	nu = r - div - (sigma**2 / 2)
	#Probabilities
	pu = -0.5 * dt * ( (sigma / dp )**2+ (nu/dp) )
	pm = 1 + dt*(sigma / dp)**2 + r*dt 
	pd = -0.5 * dt *( (sigma/ dp )**2 - (nu/dp) )

	#Stock prices
	SVDown = [S*np.exp(-(x+1)*dp)  for x in range(Np)][::-1] #Stock Price Vector Down
	SVUp = [S*np.exp((x+1)*dp)  for x in range(Np)]
	SV = SVDown + [S] + SVUp

	CM = adjustedMatrix(Nt, Np) #Option Price Matrix

	#Set Option value at maturity and boundary conditions

	for i in range(Np + 1):
		CM[Nt, i] = max(0, SV[i + Np] - K)
		#print(i, SV[i+ Np], CM[Nt, i])
		CM[Nt, -i] = max(0, SV[Np - i] - K)
		#print(-i, SV[Np - i], CM[Nt, -i])
	lambdaU = SV[-1] - SV[-2]
	lambdaD = 0

	#Back solve through matrix using tri diagonal
	count = 0
	for i in range(Nt - 1, -1,-1):
		count += 1
		linearEqSolutions = CM[i + 1][::-1] #flip the array so the higer prices are near 0

		linearEqSolutions[0] = lambdaU
		linearEqSolutions[-1] = lambdaD
		#print(linearEqSolutions)
		
		a = [pu for _ in range(len(linearEqSolutions))]
		b = [pm for _ in range(len(linearEqSolutions))]
		c = [pd for _ in range(len(linearEqSolutions))]

		a[0] = 0
		b[0] = 1
		c[0] = -1

		a[-1] = 1
		b[-1] = -1
		c[-1] = 0
		
		solution = solveTriDiagonal(a,b,c,linearEqSolutions )
		
		CM[i] = solution[::-1]


	if call:
		return CM[0,0]
	else:
		return CM[0,0] + K*np.exp(-r*T) - S

def CNFiniteDifference(Nt,Np,dp,S,K,sigma,div, r,T, call = True):
	#Crank Nicholson
	'''
	Nt - Number of time steps to discritize
	Np - Number of prices to discritize on each size of the current point
	dp - Change in price per discrete price point
	S - Current Price
	K - Strike
	sigma - volitility
	div - dividend
	'''
	dt = T/Np
	nu = r - div - (sigma**2 / 2)
	#Probabilities
	pu = -0.25 * dt * ( (sigma / dp )**2 + (nu/dp) )
	pm = 1 + 0.5*dt*(sigma / dp)**2 + 0.5*r*dt 
	pd = -0.25 * dt *( (sigma/ dp )**2 - (nu/dp) )

	#Stock prices
	SVDown = [S*np.exp(-(x+1)*dp)  for x in range(Np)][::-1] #Stock Price Vector Down
	SVUp = [S*np.exp((x+1)*dp)  for x in range(Np)]
	SV = SVDown + [S] + SVUp

	CM = adjustedMatrix(Nt, Np) #Option Price Matrix

	#Set Option value at maturity and boundary conditions
	
	for i in range(Np + 1):
		CM[Nt, i] = max(0, SV[i + Np] - K)
		CM[Nt, -i] = max(0, SV[Np - i] - K)

	lambdaU = SV[-1] - SV[-2]
	lambdaD = 0
	


	#Back solve through matrix using tri diagonal

	for i in range(Nt - 1, -1,-1):

		forwardPrices = CM[i + 1][::-1] #flip the array so the higer prices are near 0
		linearEqSolutions = [lambdaU]
		
		for j in range(1,len(forwardPrices) - 1):
			linearEqSolutions.append(-pu*forwardPrices[j+1] - (pm - 2)*forwardPrices[j]-pd*forwardPrices[j-1])

		linearEqSolutions.append(lambdaD)

		#print(linearEqSolutions)
		a = [pu for _ in range(len(linearEqSolutions))]
		b = [pm for _ in range(len(linearEqSolutions))]
		c = [pd for _ in range(len(linearEqSolutions))]

		a[0] = 0
		b[0] = 1
		c[0] = -1

		a[-1] = 1
		b[-1] = -1
		c[-1] = 0

		solution = solveTriDiagonal(a,b,c,linearEqSolutions)
		CM[i] = solution[::-1]

	if call:
		return CM[0,0]
	else:
		return CM[0,0] + K*np.exp(-r*T) - S
