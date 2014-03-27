import numpy as np, scipy.optimize as opt, pylab

a_array = (150000-100000)*np.random.random(10)+100000
b_array = np.random.random(10)
print a_array, b_array
for i in range(10):

    a = a_array[i]
    b = b_array[i]

    rmax = 10000000.0
    dm = 1.7

    n_points = 765+i*572

    r = np.empty(n_points)

    for n in range(n_points):
        def mass(r):
            return 4*np.pi*a*(b**3)*(np.log(1+r/b)-r/(r+b)) - n*dm
        r[n] = opt.bisect(mass,0,rmax)

    theta = np.arccos((2*np.random.random(n_points)-1))
    phi = 2*np.arccos((2*np.random.random(n_points)-1))
    x,y,z = r*np.sin(theta)*np.cos(phi),r*np.sin(theta)*np.sin(phi),r*np.cos(theta)
    
    data = np.dstack((x,y,z))[0]
    np.savetxt('fake_'+str(i)+'_halo.dat',data,delimiter=',')

