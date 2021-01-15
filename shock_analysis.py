# Functions for finding shock in output of hydro code
# Created by Tom Hilder for PHS3350 fallback supernovae project

def read_output(finalstep, path=''):

    import numpy as np

    steps = np.arange(0,finalstep,1000)
    u = [0]*len(steps)

    for i in range(len(steps)):
        file = path + 'output' + str(steps[i]).zfill(8) + '.dat'
        u[i] = np.loadtxt(file)

    return np.array(u)

def locate_shock(velocities):

    import numpy as np

    dv = np.gradient(velocities)
    mag_dv = np.absolute(dv)

    shock_ind = 0
    for i in range(len(mag_dv)-1):
        if mag_dv[i+1] >= 2E7:
            shock_ind = i

    return shock_ind

def extract_shocks(finalstep, path=''):

    import numpy as np

    u = read_output(finalstep, path)

    radius = u[0,:,0]           # [cm]
    velocity = u[:,:,2]         # [cm/s]
    pressure = u[1:,:,8]         # [Ba]

    time = u[:,0,7]             # [s]

    shock_indices = []
    shock_positions = []

    for i in range(len(time)):

        shock_ind = locate_shock(velocity[i,:])

        shock_indices.append(shock_ind)
        shock_positions.append(radius[shock_ind])

    return shock_positions, shock_indices

def velocity_jumps(finalstep, path=''):

    import numpy as np

    u = read_output(finalstep, path)

    time = u[:,0,7]             # [s]
    velocity = u[:,:,2]         # [cm/s]

    shock_indices = extract_shocks(finalstep, path)[1]

    v_jumps = []
    for i in range(len(time)):
        j = shock_indices[i]
        dv = velocity[i,j+3] - velocity[i,j-3]
        v_jumps.append(dv)

    return np.abs(np.array(v_jumps))

def downstream_quantity(finalstep, quantity, path=''):

    import numpy as np

    u = read_output(finalstep, path)

    time = u[:,0,7]         # [s]

    radius = u[:,:,0]         # [cm]
    density = u[:,:,1]        # [g/cm^3]
    velocity = u[:,:,2]       # [cm/s]
    energy = u[:,:,4]         # [erg]
    phi = u[:,:,5]            # [erg]
    g_acc = u[:,:,6]          # [cm/s^2]
    pressure = u[:,:,8]       # [Ba]
    c_s = u[:,:,9]            # [cm/s]
    mass_coord = u[:,:,10]    # [g]
    energy_tot = u[:,:,11]    # [erg]
    e_pos = u[:,0,12]         # [erg]

    shock_indices = extract_shocks(finalstep, path)[1]

    target = eval(quantity)

    downstream_quantity = []

    for i in range(len(time)):
        j = shock_indices[i]

        downstream_quantity_i = target[i,j-3]
        downstream_quantity.append(downstream_quantity_i)

    return downstream_quantity

def upstream_quantity(finalstep, quantity, path=''):

    import numpy as np

    u = read_output(finalstep, path)

    time = u[:,0,7]           # [s]

    radius = u[:,:,0]         # [cm]
    density = u[:,:,1]        # [g/cm^3]
    velocity = u[:,:,2]       # [cm/s]
    energy = u[:,:,4]         # [erg]
    phi = u[:,:,5]            # [erg]
    g_acc = u[:,:,6]          # [cm/s^2]
    pressure = u[:,:,8]       # [Ba]
    c_s = u[:,:,9]            # [cm/s]
    mass_coord = u[:,:,10]    # [g]
    energy_tot = u[:,:,11]    # [erg]
    e_pos = u[:,0,12]         # [erg]

    shock_indices = extract_shocks(finalstep, path)[1]

    target = eval(quantity)

    upstream_quantity = []

    for i in range(len(time)):
        j = shock_indices[i]

        upstream_quantity_i = target[i,j+2]
        upstream_quantity.append(upstream_quantity_i)

    return upstream_quantity

def shock_animation(finalstep, path=''):

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
    velocity = u[:,:,2]       # [cm/s]
    time = u[:,0,7]           # [s]

    shock_positions = extract_shocks(finalstep, path)[0]
    v_jumps = velocity_jumps(finalstep, path)

    fig = plt.figure()
    camera = Camera(fig)
    for i in range(len(steps)-1):
        t = plt.plot(radius[i], velocity[i], color='green')
        t = plt.vlines(shock_positions[i], -1E9, 1E9, color='r')
        #t = plt.vlines(shock_positions[i], 0, -1*v_jumps[i], color='r')
        plt.xlabel(r'Radius [cm]')
        plt.ylabel(r'Velcity [cm/s]')
        plt.xscale('log')
        #plt.legend(t, [f'time = {time[i]} s'])
        camera.snap()

    animation = camera.animate()
    animation.save(path + 'shocks_animation.mp4')

def cross_sec_area(radius):

    import numpy as np

    return 4*np.pi*radius**2

def shock_eqn(gamma, A, rho0, c0, dv):

     import numpy as np

     return ((gamma + 1) * A * rho0 * dv**3) / (12 * c0)

def eval_shock_eqn_RHS(finalstep, path=''):

     import numpy as np

     gamma = 5/3

     shock_radius = np.array(extract_shocks(finalstep, path)[0])
     A = cross_sec_area(shock_radius)

     #rho0 = np.array(upstream_quantity(finalstep, 'density', path))
     #c0 = np.array(upstream_quantity(finalstep, 'c_s', path))

     rho0 = np.array(downstream_quantity(finalstep, 'density', path))
     c0 = np.array(downstream_quantity(finalstep, 'c_s', path))

     dv = velocity_jumps(finalstep, path)

     return shock_eqn(gamma, A, rho0, c0, dv)

def eval_shock_eqn_LHS(finalstep, path=''):

    import numpy as np

    u = read_output(finalstep, path)
    e_pos = u[:,0,12]
    e_kinetic = u[:,0,13]
    wave_lum = wave_luminosity(finalstep, path)

    wave_energy = e_pos # remove -1 factor from return line if you use this
    #wave_energy = e_kinetic
    #wave_energy = wave_lum
    shock_radius = np.array(extract_shocks(finalstep, path)[0])

    return (np.gradient(wave_energy, shock_radius))
    #return wave_energy, shock_radius

def shock_eqn_plot(finalstep, regime_line=False, path=''):

    import numpy as np
    import matplotlib.pyplot as plt

    u = read_output(finalstep, path)
    time = u[:,0,7]

    LHS = eval_shock_eqn_LHS(finalstep, path)
    RHS = 0.5*eval_shock_eqn_RHS(finalstep, path)

    time2 = time[LHS > 0]
    LHS2 = LHS[LHS > 0]
    RHS2 = RHS[LHS > 0]

    dif = (RHS - LHS)/RHS *100
    dif2 = (RHS2 - LHS2)/RHS2 *100

    if regime_line == True:

        time = u[:,0,7]
        e_kinetic = u[:,0,13]
        e_binding_ahead = np.abs(total_binding_ahead_shock(finalstep, path))

        index = index_array2_larger(e_binding_ahead, e_kinetic)
        regime_time = time[index]

    plt.close('all')

    fig, (ax1, ax2) = plt.subplots(2, sharex=True, gridspec_kw={'height_ratios': [3, 1]})

    #ax1.plot(time, LHS, color='red', label=r'$- \frac{dE_w}{dr}$')
    ax1.plot(time2, LHS2, color='red', label=r'$- \frac{dE_w}{dr}$')

    ax1.plot(time, RHS, color='blue', label=r'$\frac{1}{2} \, \frac{(\gamma + 1) A \rho_0 (\Delta v)^3}{12 c_0}$')

    if regime_line == True:
        ax1.vlines(regime_time,  1E37, 1E42, color='green', label=r'Time: kinetic_E > binding_E')

    ax1.legend(loc='best', fontsize='10')
    ax1.set_yscale('symlog')

    #ax2.plot(time, dif, color='orange')
    ax2.plot(time2, dif2, color='orange')

    ax2.set_ylabel('% Difference', fontsize='10')
    ax2.set_xlabel('Time [s]', fontsize='10')

    plt.savefig(path + 'shock_eqn_epos' + '.png', dpi=200)
    plt.show()

def initial_final_explosion_energy(finalstep):

    # ADD PATH AS AN INPUT

    import numpy as np
    import matplotlib.pyplot as plt

    u0_50 = read_output(finalstep, 'Mach 0.5/')
    u0_75 = read_output(finalstep, 'Mach 0.75/')
    u1_00 = read_output(finalstep, 'Mach 1.0/')
    u1_25 = read_output(finalstep, 'Mach 1.25/')
    u1_50 = read_output(finalstep, 'Mach 1.5/')
    u1_75 = read_output(finalstep, 'Mach 1.75/')
    u2_00 = read_output(finalstep, 'Mach 2.0/')


    #energies = [u0_50[:,0,12],u0_75[:,0,12],u1_00[:,0,12],u1_25[:,0,12],u1_50[:,0,12]]
    #energies = [u1_00[:,0,12],u1_25[:,0,12],u1_50[:,0,12],u1_75[:,0,12],u2_00[:,0,12]] # e_pos

    energies = [u1_00[:,0,13],u1_25[:,0,13],u1_50[:,0,13],u1_75[:,0,13],u2_00[:,0,13]] # e_kinetic

    init_energies = []
    final_energies = []

    for j in energies:
        init_energies.append(np.amax(j[0:8]))
        final_energies.append(j[-1])

    #init_mach = [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]
    init_mach = [1.0, 1.25, 1.5, 1.75, 2.0]

    plt.close('all')

    plt.plot(init_mach, init_energies, c='blue', label='Initial Explosion Energy')
    plt.plot(init_mach, final_energies, c='red', label='Final Explosion Energy')

    plt.xlabel('Initial Mach No.')
    plt.ylabel('Energy [erg]')

    #plt.yscale('log')

    plt.legend(loc='best')

    plt.savefig('init_final_energy_ekin.png')
    plt.show()

    return init_energies, final_energies

def binding_ahead_shock(time_index, u, shock_indices):

    import numpy as np

    density = u[time_index,:,1]
    energy_dens = u[time_index,:,4]
    phi = u[time_index,:,5]
    dv = u[time_index,:,14]

    shock_index = shock_indices[time_index]

    binding_energies_array = (energy_dens + density*phi)*dv

    binding_energy_ahead = np.sum(binding_energies_array[shock_index:])

    return binding_energy_ahead

def total_binding_ahead_shock(finalstep, path=''):

    import numpy as np

    u = read_output(finalstep, path)
    time = u[:,0,7]

    shock_indices = extract_shocks(finalstep, path)[1]

    energies = np.zeros(len(time))
    for i in range(0, len(time)):
        energies[i] = binding_ahead_shock(i, u, shock_indices)

    return energies

def shock_velocity(finalstep, path=''):

    import numpy as np

    u = read_output(finalstep, path)

    time = u[:,0,7]             # [s]
    velocity = u[:,:,2]         # [cm/s]

    shock_indices = extract_shocks(finalstep, path)[1]

    shock_v = np.zeros(len(time))
    for i in range(len(time)):
        j = shock_indices[i]
        if j-4 >= 0 and j+2 < len(velocity[0,:]):
            shock_v[i] = np.amax(velocity[i,j-4:j+2])
        elif j-4 < 0:
            shock_v[i] = np.amax(velocity[i,0:j+2])
        elif j+2 > len(velocity[0,:]):
            shock_v[i] = np.amax(velocity[i,j-4:len(velocity[0,:])])
        else:
            print('Error, looking for shock outside of radii')
            return 'Error'

    return shock_v

def escape_velocity(solar_masses, radius):

    import numpy as np

    G = 6.674E-8 # Gravitational constant in cgs
    M = 1.988E33*solar_masses # Mass in grams

    return np.sqrt(2*G*M/radius)

def shock_velocity_plot(finalstep, path=''):

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import (AutoLocator, AutoMinorLocator)

    u = read_output(finalstep, path)

    time = u[:,0,7]

    shock_radii = np.array(extract_shocks(finalstep, path)[0])
    shock_velocities = shock_velocity(finalstep, path)
    escape_velocities = escape_velocity(32, shock_radii)

    fig, ax = plt.subplots(constrained_layout=True)
    ax.plot(shock_radii, shock_velocities, label='Shock Velocity')
    ax.plot(shock_radii, escape_velocities, label='Escape Velocity')
    ax.set_xlabel('Shock Radius [cm]')
    ax.set_ylabel('Velocity [cm/s]')
    ax.set_yscale('log')

    def forward(x):
        return np.interp(x, shock_radii, time)

    def inverse(x):
        return np.interp(x, time, shock_radii)

    secax = ax.secondary_xaxis('top', functions=(forward, inverse))
    secax.xaxis.set_minor_locator(AutoMinorLocator())
    secax.set_xlabel('Time [s]')

    """
    plt.plot(shock_radii, shock_velocities, label='Shock Velocity')
    plt.plot(shock_radii, escape_velocities, label='Escape Velocity')
    plt.legend(loc='best')
    plt.yscale('log')
    plt.xlabel('Shock Radius [cm]')
    plt.ylabel('Velocity [cm/s]')
    """

    plt.savefig(path + 'shock_vel_rad_esc.png', dpi=150)
    plt.show()


def index_array2_larger(array1, array2):
    for i in range(len(array1)):

        a = array1[i]
        b = array2[i]

        if b > a:
            return i

def wave_luminosity(finalstep, path=''):

    import numpy as np

    shock_v = shock_velocity(finalstep, path)
    c0 = np.array(upstream_quantity(finalstep, 'c_s', path))
    rho0 = np.array(upstream_quantity(finalstep, 'density', path))

    shock_radius = np.array(extract_shocks(finalstep, path)[0])
    A = cross_sec_area(shock_radius)

    return A*rho0*c0*shock_v**2

def all_shock_eqn_plots():
    shock_eqn_plot(80000, False, 'Mach 0.5/')
    shock_eqn_plot(94000, True, 'Mach 0.75/')
    shock_eqn_plot(85000, True, 'Mach 1.0/')
    shock_eqn_plot(85000, True, 'Mach 1.25/')
    shock_eqn_plot(83000, True, 'Mach 1.5/')
    shock_eqn_plot(91000, True, 'Mach 1.75/')
    shock_eqn_plot(182000, True, 'Mach 2.0/')

def phi_grav_out_array(finalstep, path=''):
    # INCOMPLETE

    import numpy as np

    u = read_output(finalstep, path)

    G = 6.674E-8 # Gravitational constant in cgs

    radius = u[:,:,0]          # [cm]
    density = u[:,:,1]        # [g/cm^3]
    velocity = u[:,:,2]       # [cm/s]
    time = u[:,0,7]           # [s]

    shock_radius_indicies = np.array(extract_shocks(finalstep, path)[1])

    #for one time step

    radii_step = radius[0,:]
    density_step = density[0,:]
    shock_index = shock_radius_indicies[0]

    #for i in range(len(radii_step[shock_index:])):








########
