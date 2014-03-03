import numpy as np, pylab, sys

n_haloes = int(sys.argv[1])+1

mass_element = 1.7
    
n_iterations = 50000 

a,b,c = 2.85739253e+02,9.83269496e-01,3.79708464e-02

def f(x):
    return a*x**(-b)+c
rho_array = np.arange(500.0,5000.0,(5000.0-500.0)/n_haloes)
rs_array = f(rho_array)
rad_array = 30*rs_array

np.savetxt('parameters.dat', zip(*[range(1,n_haloes),rho_array,rs_array]),delimiter=',')

def NFW_density(r,rho_0,R_s):

    return rho_0/((r/R_s)*(1+(r/R_s))**2.0)

def NFW_mass(r,rho_0,R_s):

	return 4.0*np.pi*rho_0*(R_s**3.0)*(np.log((R_s+r)/(R_s))-((r)/(R_s + r)))

def mpr(r,rho_0,R_s):

    return NFW_density(r,rho_0,R_s) * 4.0  * np.pi * r**2.0
 
for n in range(1,n_haloes):

    R_s = rs_array[n]
    rho_0 = rho_array[n]
    R_max = rad_array[n]

    r_walk = np.empty((0))
    r_0 = 0.000001
    r_walk = np.append(r_walk,r_0)

    for i in range(n_iterations):
        r_prime = np.random.normal(r_walk[i], 0.01) 
        while (r_prime > R_max):
            r_prime = np.random.normal(r_walk[i], 0.01) 
    
        alpha = mpr(r_prime,rho_0,R_s)/mpr(r_walk[i],rho_0,R_s)

        if(alpha >= 1.0):

            r_walk  = np.append(r_walk,r_prime)
    
        else:

            beta = np.random.random()

            if(beta <= alpha):

                r_walk = np.append(r_walk,r_prime)

            else:

                r_walk = np.append(r_walk,r_walk[i])


    
    histo = np.histogram(r_walk,int(1000 * R_max/rad_array[1]))


    r_values = np.delete(histo[1],0)
    
    norm = NFW_mass(R_max,rho_0,R_s)/(float(sum(histo[0])))
    particles = histo[0]*norm
    
    radius = []
    n_points = int(sum(np.rint(particles)))

    print sum(particles),sum(np.rint(particles))

    for i in range(len(particles)):
        for j in range(int(np.rint(particles[i]))):
            radius.append(r_values[i])
    phi = 2.0 * np.pi * np.random.random(n_points)
    theta = np.arccos(2.0 * np.random.random(n_points) - 1.0)

    x = radius * np.sin(theta) * np.cos(phi)
    y = radius * np.sin(theta) * np.sin(phi)
    z = radius * np.cos(theta)

    np.savetxt('halo_'+str(n)+'_fake.dat', zip(*[x,y,z]),delimiter=',')
    print 'Generated',n,'of',n_haloes-1,'halos'
