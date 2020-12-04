# Functions for plotting output of hydro code
# Created by Tom Hilder for PHS3350 fallback supernovae project

def plot_density(step=0, path='', ev=False):

    import numpy as np
    import matplotlib.pyplot as plt

    file = path + 'output' + str(step).zfill(8) + '.dat'

    u = np.loadtxt(file)

    plt.plot(u[:,0]/1E5,u[:,1])

    plt.xscale('log')
    plt.yscale('log')

    plt.xlabel(r'Radius [km]')
    plt.ylabel(r'Density [$g \, cm^{-3}$]')

    if ev == False:
        plt.savefig(path + 'density' + str(step).zfill(8) + '.png', dpi=150)
        plt.show()

def plot_velocity(step=0, path='', ev=False):

    import numpy as np
    import matplotlib.pyplot as plt

    file = path + 'output' + str(step).zfill(8) + '.dat'

    u = np.loadtxt(file)

    plt.plot(u[:,0]/1E5,u[:,2])

    plt.xscale('log')

    plt.xlabel(r'Radius [km]')
    plt.ylabel(r'Velocity [$cm \, s^{-1}$]')

    if ev == False:
        plt.savefig(path + 'velocity' + str(step).zfill(8) + '.png', dpi=150)
        plt.show()

def plot_energy(step=0, path='', ev=False):

    import numpy as np
    import matplotlib.pyplot as plt

    file = path + 'output' + str(step).zfill(8) + '.dat'

    u = np.loadtxt(file)

    plt.plot(u[:,0]/1E5,u[:,4])

    plt.xscale('log')
    plt.yscale('log')

    plt.xlabel(r'Radius [km]')
    plt.ylabel(r'Energy [erg]')

    if ev == False:
        plt.savefig(path + 'energy' + str(step).zfill(8) + '.png', dpi=150)
        plt.show()

def evolution_2D_plot(finalstep, stepsize, type='density', path=''):

    import numpy as np
    import matplotlib.pyplot as plt

    steps = np.arange(0,finalstep,stepsize)

    for i in steps:
        eval('plot_' + type)(i, path, ev=True)

    plt.savefig(path + type + '_2D_evolution' + '.png', dpi=150)
    plt.show()
    plt.close()

def evolution_3D_plot(finalstep, type, path=''):

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    from mpl_toolkits.mplot3d import Axes3D

    steps = np.arange(0,finalstep,1000)
    u = [0]*len(steps)

    for i in range(len(steps)):
        file = path + 'output' + str(steps[i]).zfill(8) + '.dat'
        u[i] = np.loadtxt(file)

    u = np.array(u)

    radii = u[0,:,0]
    density = u[:,:,1]
    velocity = u[:,:,2]
    energy = u[:,:,4]

    time = u[:,0,7]

    #radii, steps = np.meshgrid(radii, steps)
    radii, time = np.meshgrid(radii, time)

    if type == 'density':
        Z = np.log(density)
        Zax_label = r'$\log_{10} \left( \frac{Density}{1 \, g \, cm^{-3}} \right)$'
    elif type == 'velocity':
        Z = velocity
        Zax_label = r'$ Velocity \, [cm \, s^{-1}]$'
    elif type == 'energy':
        Z = np.log(energy)
        Zax_label = r'$\log_{10} \left( \frac{Energy}{1 \, erg} \right)$'
    else:
        return 'Please choose either density, velocity or energy for type'

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(np.log10(radii), time, Z, cmap=plt.cm.coolwarm, linewidth=0, antialiased=True)
    fig.colorbar(surf, shrink=0.5, aspect=5)

    ax.set_xlabel(r'$\log_{10} \left( \frac{Radius}{1 \, cm} \right)$', fontsize=10)
    ax.set_ylabel(r'Time [s]', fontsize=10)
    ax.set_zlabel(Zax_label, fontsize=10)

    plt.savefig(path + type + '_3D_evolution' + '.png', dpi=150)
    plt.show()

def animation(finalstep, type, path=''):

    import numpy as np
    import matplotlib.pyplot as plt
    from celluloid import Camera

    steps = np.arange(0,finalstep + 1000,1000)
    u = [0]*len(steps)

    for i in range(len(steps)):
        file = path + 'output' + str(steps[i]).zfill(8) + '.dat'
        u[i] = np.loadtxt(file)

    u = np.array(u)

    radii = u[0,:,0]
    density = u[:,:,1]
    velocity = u[:,:,2]
    energy = u[:,:,4]

    time = u[:,0,7]

    if type == 'density':
        Y = density
        Yax_label = r'Density [$g \, cm^{-3}$]'
    elif type == 'velocity':
        Y = velocity
        Yax_label = r'$ Velocity \, [cm \, s^{-1}]$'
    elif type == 'energy':
        Y = energy
        Yax_label = r'Energy [erg]'
    else:
        return 'Please choose either density, velocity or energy for type'

    fig = plt.figure()
    camera = Camera(fig)
    for i in range(len(steps)):
        t = plt.plot(radii/1E5, Y[i], color='green')
        plt.xlabel(r'Radius [km]')
        plt.ylabel(Yax_label)
        plt.xscale('log')
        if type == 'density' or type == 'energy':
            plt.yscale('log')
        plt.legend(t, [f'time = {time[i]} s'])
        camera.snap()

    animation = camera.animate()
    animation.save(path + type + 'animation.mp4')

def make_everything(finalstep, path='', animations=True):

    import matplotlib.pyplot as plt
    plt.close('all')

    evolution_2D_plot(finalstep, 5000, 'density', path)
    evolution_2D_plot(finalstep, 5000, 'velocity', path)
    evolution_2D_plot(finalstep, 5000, 'energy', path)

    evolution_3D_plot(finalstep, 'density', path)
    evolution_3D_plot(finalstep, 'velocity', path)
    evolution_3D_plot(finalstep, 'energy', path)

    if animations == True:
        animation(finalstep, 'density', path)
        animation(finalstep, 'velocity', path)
        animation(finalstep, 'energy', path)

















































#
#
#
