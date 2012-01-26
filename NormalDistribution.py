import math
import numpy as np
import matplotlib.pyplot as plt

def normal(mu,sigma):
    def f(x):
        z = 1.0*(x-mu)/sigma
        e = math.e**(-0.5*z**2)
        C = math.sqrt(2*math.pi)*sigma
        return 1.0*e/C
    return f

X = 16.
dx = 1.
R = np.arange(-X,X+dx,dx)

L = list()
sd = 6
f = normal(mu=0,sigma=sd)
L.append([f(x) * (X-1) for x in R])
colors = ['b']

for c,P in zip(colors,L):
    plt.plot(R,P,zorder=1,color='0.2',lw=1.5)
    #plt.scatter(R,P,zorder=2,s=50,color=c)
    
ax = plt.axes()
ax.set_xlim(-16.,16.)
#ax.set_ylim(-0.01,0.5)

y = (np.negative(np.power(R,2)) + 256.0) / 256.0
plt.plot(R,y)

plt.show()

#plt.savefig('example.png')
