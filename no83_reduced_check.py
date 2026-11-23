import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
alpha=np.arcsin(1/np.sqrt(3))
y0=np.array([0.,alpha,np.pi-alpha,-alpha])
def rhs(g):
 def f(t,v):
  h,x,y,z=v
  H=(np.sin(x-h)+np.sin(y-h)+2*np.sin(z-h))/5
  if g==0:
   return [H,(np.sin(h-x)+2*np.sin(z-x))/5,(np.sin(h-y)+2*np.sin(z-y))/5,(np.sin(h-z)+np.sin(x-z)+np.sin(y-z))/5]
  return [H,(np.sin(h-x)+np.sin(y-x))/5,(np.sin(h-y)+np.sin(x-y))/5,np.sin(h-z)/5]
 return f
s0=solve_ivp(rhs(0),(0,3.12),y0,method='DOP853',rtol=2e-13,atol=2e-14,dense_output=True)
s1=solve_ivp(rhs(1),(0,3.12),y0,method='DOP853',rtol=2e-13,atol=2e-14,dense_output=True)
def D(t):return s1.sol(t)[0]-s0.sol(t)[0]
for t in [3.10,3.12]: print(f'D({t})={D(t):+.15f}')
print('root',brentq(D,3.10,3.12,xtol=1e-13))
print('NOTE: numerical check only, not validated integration')
