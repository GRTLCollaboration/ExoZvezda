import numpy as np
import matplotlib.pyplot as plt

v = '30'

x = []
y = []
data = []
data2 = []
with open('StarTracking.dat','r') as f:
	d = f.readlines()
	for i in d:
		k = i.rstrip().split(" ")
		for element in k:
			if element != ' ' and element != '':
				try:
					data.append(float(element))
				except:
					pass

'''with open('inspiral/cont/Weyl_integral_22','r') as f:
	d = f.readlines()
	for i in d:
		k = i.rstrip().split(" ")
		for element in k:
			if element != ' ' and element != '':
				try:
					data2.append(float(element))
				except:
					pass'''


data = np.array(data[4:])
'''data = np.array(data)'''
data2 = np.array(data2)

for i in range(int(data.size/4)):
	x.append(data[4*i+1])
	y.append(data[4*i+3])

x = np.array(x)
y = np.array(y)
print(x.size)
print(y.size)


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

#ax.set_ylim([0,1.2])
ax.plot(x,y,'k',lw=0.7)
plt.rc('text', usetex=True)
plt.title(r'Gravitational Wave Signal in 22 Spherical Harmonic')
plt.xlabel(r'Time $t \cdot m$')
plt.ylabel(r'Newman-Penrose Scalar $\psi_4$')
plt.show()
fig.savefig('inspiral64v'+v+'/v'+v+'Weyl22.png')
