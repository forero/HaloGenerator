import numpy
import random
import pylab 
import timeit

start = timeit.default_timer()

data = numpy.loadtxt(open('parameters.dat', 'r'))
R_s = data[0]
rho_0 = data[1]
mass_element = 1.7
R_max = 2.5
n_points = 1000

def nfw_density(r,R_s,rho_0):

	return rho_0/((r/R_s)*(1+(r/R_s))**2.0)
 
def nfw_mass(r,R_s,rho_0):

	return 4.0*numpy.pi*rho_0*(R_s**3.0)*(numpy.log((R_s+r)/(R_s))-((r)/(R_s + r)))
norm = nfw_mass(R_max,R_s,rho_0)
n_particles = int(norm/mass_element)
dr = R_max/n_points

r = [i*dr for i in range(1,n_points+1)]
m = nfw_mass(r,R_s,rho_0)
density = nfw_density(r,R_s,rho_0)
delta_m = [4 * numpy.pi * r[i] * r[i] * density [i] for i in range(n_points)]

dist = [int(n_particles*(d_m/norm)*dr) for d_m in delta_m]

pylab.plot(r , dist,'k.')
pylab.title('Particle Distribution')
pylab.xlabel('Radius')
pylab.ylabel('Particles')
pylab.savefig('distribution.png')
pylab.close()

radius = []

for i in range(n_points):

	for j in range(dist[i]):

		radius.append(r[i])

phi = [2.0 * numpy.pi * random.random() for i in range(len(radius))]
theta = [numpy.arccos(2.0 * random.random() - 1.0) for i in range(len(radius))]

x = radius * numpy.sin(theta) * numpy.cos(phi)
y = radius * numpy.sin(theta) * numpy.sin(phi)
z = radius * numpy.cos(theta)

print 'Plotting halo...',
pylab.plot(x , y, 'k,')
pylab.xlabel('$x$')
pylab.ylabel('$y$')
pylab.title('x vs y')
pylab.savefig('halo_xy.png')
pylab.close()

pylab.plot(y , z, 'k,')
pylab.xlabel('$y$')
pylab.ylabel('$z$')
pylab.title('y vs z')
pylab.savefig('halo_yz.png')
pylab.close()

pylab.plot(z , x, 'k,')
pylab.xlabel('$z$')
pylab.ylabel('$x$')
pylab.title('z vs x')
pylab.savefig('halo_zx.png')
pylab.close()

print 'Done'

tabla = open('coordinates.dat', "w")
tabla.write('\n'.join('%lf,%lf,%lf' % (x[i],y[i],z[i]) for i in range(len(x))))

stop = timeit.default_timer()
print 'Done in ' + str(stop - start)+' seconds' 


