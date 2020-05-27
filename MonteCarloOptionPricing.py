import numpy as np
from scipy.stats import norm


'''
S = 50.00
K = 52.00
sigma = 0.20
r = 0.02
tau = 0.50
'''

dt = 1 / 252 #simulate on daily basis
trials = 10**4

	
def generatePath(S, K, sigma, r, tau):
	St = S
	t = 0
	path = [St]
	while t < tau:
		ds = St * (r*dt + sigma*np.sqrt(dt)*np.random.normal())
		St += ds
		path.append(St)
		t+= dt
	return path

def EuropeanOptionPathPrice(path, K, r, tau):
	return max(0, path[-1] - K)* np.exp(-r*tau)

def EuropeanOption(S, K, sigma, r, tau):
	valueSum = 0
	for _ in range(trials):
		path = generatePath(S, K, sigma, r, tau)
		value = EuropeanOptionPathPrice(path, K, r, tau)
		valueSum += value
	return valueSum / trials

def AsianOptionPathPrice(path, K, r, tau):
	return max(0, mean(path)-K)* np.exp(-r*tau)
	
def AsianOption(S, K, sigma, r, tau):
	valueSum = 0
	for _ in range(trials):
		path = generatePath(S, K, sigma, r, tau)
		value = AsianOptionPathPrice(path, K, r, tau)
		valueSum += value
	return valueSum / trials 

def KnockOutBarrierOptionPathPrice(path, K,r,tau, LBarrier):
	if min(path) < LBarrier:
		return 0
	return max(0, path[-1] - K) * np.exp(-r*tau)
	
def KnockOutBarrierOption(S, K, sigma, r, tau,LBarrier):
	valueSum = 0
	for _ in range(trials):
		path = generatePath(S, K, sigma, r, tau)
		value = KnockOutBarrierOptionPathPrice(path,K,r,tau,LBarrier)
		valueSum += value
	return valueSum / trials

def AmericanOptionPrice(paths, K,r,tau):
	#least square approach proposed by Longstaff and Schwartz (2001)
	paths = np.array(paths).T 
	payoffs = np.array([0.0] * trials)
	for t in range(-1, -len(paths), -1):
		if t == -1:
			for i in range(len(paths[t])):
				if paths[t][i] > K:
					payoffs[i] =  paths[t][i] - K
			
		else:
			payoffs *= np.exp(-r*dt)
			prices = []
			HVs = []
			toExercise = []
			for i in range(len(paths[t])):
				if paths[t][i] > K:
					prices.append(paths[t][i])
					HVs.append(payoffs[i])
					toExercise.append(i)
			#Regression
			if len(toExercise) > 0:
				EHV = np.polyfit(prices, HVs, 2)
				for i in toExercise:
					EH = EHV[0] * paths[t][i]**2 +  EHV[1] *paths[t][i] + EHV[2]
					if paths[t][i] - K > EH:
						payoffs[i] = paths[t][i] - K
	#get present values at time 0
	payoffs *= np.exp(-r*dt)
	return max(np.mean(payoffs), S-K)


def AmericanOption(S, K, sigma, r, tau):
	paths = []
	for _ in range(trials):
		paths.append(generatePath(S, K, sigma, r, tau))
	return AmericanOptionPrice(paths, K,r,tau)

