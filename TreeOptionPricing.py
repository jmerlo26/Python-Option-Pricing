import numpy as np


def AmericanBinomialTree(N,S,K,sigma,div, r,T, call = True):
	#solve values to create tree
	dt = T/N
	u = np.exp(sigma * np.sqrt(dt) ) 
	d = np.exp(-sigma * np.sqrt(dt) )

	P = (np.exp(r * dt) - d ) / (u - d)
	Q = 1.-P

	#Create Binomial Tree
	bTree = np.zeros( (N+1, N+1))

	#Fill in prices of Binomial Tree
	bTree[0,0] = S

	for i in range(1, N+1):
		bTree[i,0] = bTree[i-1,0] * u 
		for j in range(1,i+1):
			bTree[i,j] = bTree[i-1,j-1] * d

	#Fill in known values - at expiration
	optionTree = np.zeros( (N+1, N+1) )

	for j in range(N+1):
		if call: # Call
			optionTree[N,j] = max(0, bTree[N][j]-K)
		else: #Put
			optionTree[N,j] = max(0, K-bTree[N,j])
    
	#Back solve values
	for i in range(N-1, -1,-1):
		for j in range(i + 1):
			if call:
				optionTree[i,j] = max(0, bTree[i,j]-K, np.exp(-r*dt)*(P*optionTree[i+1,j]+Q*optionTree[i+1,j+1]))
			else:
				optionTree[i,j] =  max(0, K-bTree[i,j],  np.exp(-r*dt)*(P*optionTree[i+1,j]+Q*optionTree[i+1,j+1]) )

	return optionTree[0,0]

def EuropeanBinomialTree(N,S,K,sigma,div, r,T, call = True):
	#solve values to create tree
	dt = T/N
	u = np.exp(sigma * np.sqrt(dt) ) 
	d = np.exp(-sigma * np.sqrt(dt) )

	P = (np.exp(r * dt) - d ) / (u - d)
	Q = 1.-P

	#Create Binomial Tree
	bTree = np.zeros( (N+1, N+1))

	#Fill in prices of Binomial Tree
	bTree[0,0] = S

	for i in range(1, N+1):
		bTree[i,0] = bTree[i-1,0] * u 
		for j in range(1,i+1):
			bTree[i,j] = bTree[i-1,j-1] * d

	#Fill in known values - at expiration
	optionTree = np.zeros( (N+1, N+1) )

	for j in range(N+1):
		if call: # Call
			optionTree[N,j] = max(0, bTree[N][j]-K)
		else: #Put
			optionTree[N,j] = max(0, K-bTree[N,j])
    
	#Back solve values
	for i in range(N-1, -1,-1):
		for j in range(i + 1):
			if call:
				optionTree[i,j] = max(0, np.exp(-r*dt)*(P*optionTree[i+1,j]+Q*optionTree[i+1,j+1]))
			else:
				optionTree[i,j] =  max(0, np.exp(-r*dt)*(P*optionTree[i+1,j]+Q*optionTree[i+1,j+1]) )

	return optionTree[0,0]

def AmericanTrinomialTree(N,S,K,sigma,div, r,T, call = True):
	dt = T/N
	u = np.exp(sigma * np.sqrt(2*dt) ) 
	d = np.exp(-sigma * np.sqrt(2*dt) )
	m = 1

	pu =( (np.exp( (r-div)*dt/2) - np.exp(-sigma * np.sqrt(dt/2) ) ) / (np.exp(sigma*np.sqrt(dt/2)) - np.exp(-sigma*np.sqrt(dt/2))) )**2
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