# Functions for plotting output of hydro code
# Created by Tom Hilder for PHS3350 fallback supernovae project

def units_dict():
    """A dictionary containing axis labels with units for all plotting paramters used in other functions.
    """

    units = {
    'radius':r'Radius [cm]',
    'density':r'Density [$g \, cm^{-3}$]',
    'velocity':r'Velocity [$cm \, s^{-1}$]',
    'energy':r'Energy [erg]',
    'phi':r'Gravitational Potential [erg]',
    'g_acc':r'Gravitational Acceleration [$cm \, s^{-2}$]',
    'time':r'Time [s]',
    'pressure':r'Pressure [Ba]',
    'c_s':r'Sound Speed [$cm \, s^{-1}$]',
    'mass_coord':r'Mass Coordinate [g]',
    'energy_tot':r'Total Energy [erg]'
    }
    return units

def plot(x='radius', y='density', step=0, path='', ev=False):
    """
    INPUTS:     x = parameter to be plotted on X-axis
                y = parameter to be plotted on Y-axis (choose x and y from: radius, density, velocity, energy, phi, g_acc, time, pressure, c_s, mass_coord, energy_tot)
                step = timestep to be plotted
                path = the relative path of the output files to be plotted, and the location where the plot will be saved
                ev = variable that allows you to suppress the output of the function. That is, the plot is still created, but not saved or shown. Useful for the case
                where you would like to call this function inside another function that uses the plots created further.

    OUTPUT:     2D plot of any two parameters outputted by the code against each other, for any time step also outputted.
    """

    import numpy as np
    import matplotlib.pyplot as plt

    file = path + 'output' + str(step).zfill(8) + '.dat' # convert specified step to appropriate output file name

    u = np.loadtxt(file) # load output file as a numpy array

    radius = u[:,0]          # [cm]
    density = u[:,1]        # [g/cm^3]
    velocity = u[:,2]       # [cm/s]
    energy = u[:,4]         # [erg]
    phi = u[:,5]            # [erg]
    g_acc = u[:,6]          # [cm/s^2]
    time = u[0,7]           # [s]
    pressure = u[:,8]       # [Ba]
    c_s = u[:,9]            # [cm/s]
    mass_coord = u[:,10]    # [g]
    energy_tot = u[:,11]    # [erg]

    X = eval(x)  # set x-axis plotting variable equal to that specified
    Y = eval(y)  # set y-axis plotting variable equal to that specified

    plt.plot(X,Y)

    plt.xscale('log')
    ylog = True
    for j in Y: # set y axis scale to log, provided that no values are negative
        if j <= 0:
            ylog = False
    if ylog == True:
        plt.yscale('log')

    units = units_dict()

    plt.xlabel(units[x]) # set axes labels using dictionary containing units for each parameter
    plt.ylabel(units[y])

    if ev == False: # save and show figure if ev ==  False
        plt.savefig(path + x + '-' + y + str(step).zfill(8) + '.png', dpi=150)
        plt.show()


def evolution_2D_plot(finalstep, stepsize, x='radius', y='density', path=''):
    """
    INPUTS:     finalstep = last time step outputted by code that you would like to plot to
                stepsize = number of steps between each result you would like plotted
                x = parameter to be plotted on X-axis
                y = parameter to be plotted on Y-axis (choose x and y from: radius, density, velocity, energy, phi, g_acc, time, pressure, c_s, mass_coord, energy_tot)
                path = the relative path of the output files to be plotted, and the location where the plot will be saved

    OUTPUT:     2D plot of any two parameters outputted by the code against each other, for all time steps up until specified.
    """

    import numpy as np
    import matplotlib.pyplot as plt

    if y == 'pressure' or y == 'c_s':
        steps = np.arange(1000,finalstep,stepsize)
    else:
        steps = np.arange(0,finalstep,stepsize)

    for i in steps:
        plot(x, y, i, path, ev=True)

    plt.savefig(path + x + y + '_2D_evolution' + '.png', dpi=150)
    plt.show()

    plt.close('all')

def evolution_3D_plot(finalstep, z='density', zlog=True, path=''):
    """
    INPUTS:     finalstep = last time step outputted by code that you would like to plot to
                z = parameter to be plotted on Z-axis. Choose from: density, velocity, energy, phi, g_acc, time, pressure, c_s, mass_coord, energy_tot
                path = the relative path of the output files to be plotted, and the location where the plot will be saved
                zlog = True for logarithmic Z-axis, False otherwise

    OUTPUTS:    3D surface plot of supernova evolution, with radius on the X-axis, time on the Y-axis and the specified parameter on the Z-axis. Z-axis scale
                will be logarithmic provided that no values are zero (ie. velocity will not be logarithmic).
    """

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

    radius = u[0,:,0]           # [cm]
    density = u[:,:,1]          # [g/cm^3]
    velocity = u[:,:,2]         # [cm/s]
    energy = u[:,:,4]           # [erg]
    phi = u[:,:,5]              # [erg]
    g_acc = u[:,:,6]            # [cm/s^2]
    time = u[:,0,7]             # [s]
    pressure = u[1:,:,8]         # [Ba]
    c_s = u[1:,:,9]              # [cm/s]
    mass_coord = u[:,:,10]      # [g]
    energy_tot = u[:,:,11]      # [erg]

    if z == 'pressure' or z == 'c_s': # first timestep has zero pressure and c_s, so we avoid it in this case
        time = u[1:,0,7]

    X, Y = np.meshgrid(radius, time)
    Z = eval(z)

    X = np.log10(X)
    if zlog == True:
        Z = np.log10(Z)

    units = units_dict()

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X, Y, Z, cmap=plt.cm.coolwarm, linewidth=0, antialiased=True)
    fig.colorbar(surf, shrink=0.5, aspect=8)

    ax.set_xlabel(r'$\log_{10}$ ' + units['radius'], fontsize=10)
    ax.set_ylabel(r'Time [s]', fontsize=10)
    if zlog == True:
        ax.set_zlabel(r'$\log_{10}$ ' + units[z], fontsize=10)
    else:
        ax.set_zlabel(units[z], fontsize=10)

    plt.savefig(path + z + '_3D_evolution' + '.png', dpi=150)
    plt.show()

    plt.close('all')

def animation(finalstep, y='density', path=''):
    """
    INPUTS:     finalstep = last time step outputted by code that you would like to plot to
                y = parameter to be plotted on Y-axis. Choose from: density, velocity, energy, phi, g_acc, time, pressure, c_s, mass_coord, energy_tot
                path = the relative path of the output files to be plotted, and the location where the plot will be saved

    OUTPUTS:    2D plot of any two parameters outputted by the code against each other, animated for all time steps up until specified.
    """

    import numpy as np
    import matplotlib.pyplot as plt
    from celluloid import Camera

    steps = np.arange(0,finalstep + 1000,1000)
    u = [0]*len(steps)

    for i in range(len(steps)):
        file = path + 'output' + str(steps[i]).zfill(8) + '.dat'
        u[i] = np.loadtxt(file)

    u = np.array(u)

    radius = u[:,:,0]          # [cm]
    density = u[:,:,1]        # [g/cm^3]
    velocity = u[:,:,2]       # [cm/s]
    energy = u[:,:,4]         # [erg]
    phi = u[:,:,5]            # [erg]
    g_acc = u[:,:,6]          # [cm/s^2]
    time = u[:,0,7]           # [s]
    pressure = u[:,:,8]       # [Ba]
    c_s = u[:,:,9]            # [cm/s]
    mass_coord = u[:,:,10]    # [g]
    energy_tot = u[:,:,11]    # [erg]

    X = radius
    Y = eval(y)  # set y-axis plotting variable equal to that specified

    ylog = True
    for j in Y: # set y axis scale to log, provided that no values are negative
        for k in j:
            if k < 0:
                ylog = False
            elif k == 0:
                k = 10**(-32)

    units = units_dict()

    if y == 'pressure' or y == 'c_s':
        start = 1
    else:
        start = 0

    fig = plt.figure()
    camera = Camera(fig)
    for i in range(start, len(steps)):
        t = plt.plot(X[i], Y[i], color='green')
        plt.xlabel(r'Radius [km]')
        plt.ylabel(units[y])
        plt.xscale('log')
        if ylog == True:
            plt.yscale('log')
        plt.legend(t, [f'time = {time[i]} s'])
        camera.snap()

    animation = camera.animate()
    animation.save(path + y + 'animation.mp4')

def make_everything(finalstep, path='', animations=True):
    """Function to create many plots and animations in one go.
    """

    # density, velocity, energy, phi, g_acc, time, pressure, c_s, mass_coord, energy_tot

    import matplotlib.pyplot as plt
    plt.close('all')

    evolution_2D_plot(finalstep, 5000, 'radius', 'density', path)
    evolution_2D_plot(finalstep, 5000, 'radius', 'velocity', path)
    evolution_2D_plot(finalstep, 5000, 'radius', 'energy', path)
    evolution_2D_plot(finalstep, 5000, 'radius', 'pressure', path)
    evolution_2D_plot(finalstep, 5000, 'radius', 'c_s', path)
    evolution_2D_plot(finalstep, 5000, 'radius', 'mass_coord', path)
    evolution_2D_plot(finalstep, 5000, 'radius', 'phi', path)

    evolution_3D_plot(finalstep, 'density', True, path)
    evolution_3D_plot(finalstep, 'velocity', False, path)
    evolution_3D_plot(finalstep, 'energy', True, path)
    evolution_3D_plot(finalstep, 'pressure', True, path)
    evolution_3D_plot(finalstep, 'c_s', False, path)
    evolution_3D_plot(finalstep, 'mass_coord', False, path)

    if animations == True:
        animation(finalstep, 'density', path)
        animation(finalstep, 'velocity', path)
        animation(finalstep, 'energy', path)
        animation(finalstep, 'pressure', path)
        animation(finalstep, 'c_s', path)
        animation(finalstep, 'mass_coord', path)
        animation(finalstep, 'phi', path)

















































#
#
#
