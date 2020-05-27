import numpy as np 
import scipy

import warnings


def HestonExplicitPDE(SV, VV,Nt,kappa,theta,S,K,sigma,div, r,T):
	'''
	SV : Vector of prices
	VV : Vector of volatilities
	Nt : Number of time steps
	'''
	dt = T/Nt
	Nv = len(VV)
	Ns = len(SV)
	dv = np.abs(VV[-1] - VV [0]) / Nv
	ds = np.abs(SV[-1] - SV[0]) / Ns
	
	U = np.zeros((Ns, Nv))


	Smax = max(SV)
	#Terminal Conditions
	for s in range(Ns):
		for v in range(Nv):

			U[s,v] = max(SV[s] - K, 0);

	for t in range(Nt-1):
		for v in range(Nv-1):
			U[0,v] = 0
			U[Ns-1,v] = max(Smax - K,0)

		for s in range(Ns):
			U[s,Nv-1] = max(0, SV[s] - K);

		u = U[:] 
		for s in range(1,Ns-1):
			DerV = (u[s,1] - u[s,0]) / (VV[1] - VV[0])
			DerS = (u[s+1,0] - u[s-1,0])/(SV[s+1] - SV[s-1])

			if np.isnan(DerS):
				DerS = 0
			if np.isnan(DerV):
				DerV = 0
			dU = dt*(-r*u[s,0] + (r-div)*SV[s]*DerS + kappa*theta*DerV)
			U[s,0] = u[s,0] + dU 

		u = U[:]
		for s in range(1,Ns-1):
			for v in range(1, Nv-1):
				A = 1-dt*(s**2 * v * dv + (sigma**2 * v)/dv + r )
				B = s*dt/2*(s*v*dv - r + div) 
				C = s*dt/2*(s*v*dv + r - div) 
				D = dt/2/dv * (sigma**2 * v - kappa *(theta - v*dv))
				E = dt/2/dv * (sigma**2 * v + kappa *(theta - v*dv))
				F = dt*sigma*(s)*(v)/4
				U[s,v] = A*u[s,v] + B*u[s-1,v] + C*u[s+1,v] + D*u[s,v-1] + E*u[s,v+1] + F*(u[s+1,v+1] + u[s-1, v-1] - u[s-1, v+1] - u[s+1,v-1])
	return U

def HestonModelPrice(Nt,Ns,Nv,ds,kappa,theta,S,K,sigma,div, r,T):
	SVDown = [S*np.exp(-(x+1)*ds)  for x in range(Ns)][::-1] #Stock Price Vector Down
	SVUp = [S*np.exp((x+1)*ds)  for x in range(Ns)]
	SV = SVDown + [S] + SVUp

	kappa = 2
	theta = 0.09
	dv = 1 / Nv
	VV = [x*dv for x in range(Nv)]
	print(VV)

	HestonPDE = HestonExplicitPDE(SV=SV, VV = VV,Nt = Nt,kappa = kappa,theta = theta,S=S,K = K,sigma=sigma,div=div, r=r,T=T)
	return scipy.interpolate.griddata(SV,VV, U,s,sigma)