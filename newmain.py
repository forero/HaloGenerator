import numpy
import random
import pylab
import sys

# First it opens the file given as an argument and takes the parameters R_s and rho_0

parameters = numpy.loadtxt(open(sys.argv[1],'r'))
R_s = parameters[0]
rho_0 = parameters[1]

# Then it declares the value for the maximum radius of the halo and the value for the mass of each particle

R_max = 2.5
mass_element = 1.7

# It defines NFW profile for density

def NFW_density(r):

    return rho_0/((r/R_s)*(1+(r/R_s))**2.0)

# Also defines the mass per radius

def mpr(r):

    return NFW_density(r) * 4.0  * numpy.pi * r**2.0

# Here is where the Metropolis-Hastings Algorithm begins, first it declares an array for the values of r and appends a value

r_walk = numpy.empty((0))
r_0 = 1.0
r_walk = numpy.append(r_walk,r_0)

# Then the random walk begins, it does n_interations steps

n_iterations = 500000 
for i in range(n_iterations):
    print i
    r_prime = numpy.random.normal(r_walk[i], 0.1) 
    alpha = mpr(r_prime)/mpr(r_walk[i])

    if(alpha >= 1.0):

        r_walk  = numpy.append(r_walk,r_prime)
    
    else:

        beta = random.random()

        if(beta <= alpha):

            r_walk = numpy.append(r_walk,r_prime)

        else:

            r_walk = numpy.append(r_walk,r_walk[i])


# Now it makes an histogram that corresponds to the number of particles per each value of the radius and plots it

particles = [int(p/mass_element) for p in numpy.histogram(r_walk,10000)[0]]
r_values = numpy.delete(numpy.histogram(r_walk,10000)[1],1,0)

pylab.plot(r_values,particles,'k')
pylab.title('Particle Distribution')
pylab.xlabel('Radius')
pylab.ylabel('Number of Particles')
pylab.savefig('distribution.png')
pylab.close()

# Then it counts the total number of particles and makes an array for the radial position of each particle

n_points = sum(particles)
radius = []

for i in range(len(particles)):

	for j in range(particles[i]):

		radius.append(r_values[i])

# Generates random values for the polar and azimutal angles for each particle

phi = [2.0 * numpy.pi * random.random() for i in range(n_points)]
theta = [numpy.arccos(2.0 * random.random() - 1.0) for i in range(n_points)]

# Converts everything to cartesian coordinates and plots

x = radius * numpy.sin(theta) * numpy.cos(phi)
y = radius * numpy.sin(theta) * numpy.sin(phi)
z = radius * numpy.cos(theta)

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

# Exports

tabla = open('coordinates.dat', "w")
tabla.write('\n'.join('%lf,%lf,%lf' % (x[i],y[i],z[i]) for i in range(len(x))))
