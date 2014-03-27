import numpy as np, scipy.optimize as opt, pylab

dm =1.7
n_haloes = 5
rmax = 4

points = np.arange(500,10000,(10000-500)/n_haloes)

mtot = dm*points
b_array = np.random.random(n_haloes)
a_array = mtot/(4*np.pi*(b_array**3)*(np.log(1+rmax/b_array)-rmax/(rmax+b_array)))

parameters = np.dstack((a_array,b_array))[0]

np.savetxt('parameters.dat',parameters,delimiter=',')

for i in range(n_haloes):

    a = a_array[i]
    b = b_array[i]

    n_points = points[i]

    r = np.empty(n_points)

    for n in range(n_points):
        def mass(r):
            return 4*np.pi*a*(b**3)*(np.log(1+r/b)-r/(r+b)) - n*dm
        r[n] = opt.bisect(mass,0,rmax)
    
    theta = np.arccos(2*np.random.random(n_points)-1)
    phi =  2.0 * np.pi * np.random.random(n_points)

    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    
    data = np.dstack((x,y,z))[0]
    np.savetxt('fake_'+str(i)+'_halo.dat',data,delimiter=',')
