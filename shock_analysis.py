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

    return v_jumps

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

     rho0 = np.array(upstream_quantity(finalstep, 'density', path))

     c0 = np.array(upstream_quantity(finalstep, 'c_s', path))

     dv = np.array(velocity_jumps(finalstep, path))

     return shock_eqn(gamma, A, rho0, c0, dv)

def eval_shock_eqn_LHS(finalstep, path=''):

    import numpy as np

    u = read_output(finalstep, path)
    e_pos = u[:,0,12]
    e_kinetic = u[:,0,13]

    #wave_energy = e_pos
    wave_energy = e_kinetic
    shock_radius = np.array(extract_shocks(finalstep, path)[0])

    return -1*(np.gradient(e_pos, shock_radius))

def shock_eqn_plot(finalstep, path=''):

    import numpy as np
    import matplotlib.pyplot as plt

    u = read_output(finalstep, path)
    time = u[:,0,7]

    LHS = -1*eval_shock_eqn_LHS(finalstep, path)
    RHS = -1*eval_shock_eqn_RHS(finalstep, path)

    time2 = time[LHS > 0]
    LHS2 = LHS[LHS > 0]

    plt.close('all')
    plt.plot(time2, LHS2, color='red', label=r'$\frac{dE_w}{dr}$')
    plt.plot(time, RHS, color='blue', label=r'$\frac{-(\gamma + 1) A \rho_0 (\Delta v)^3}{12 c_0}$')
    plt.legend(loc='best', prop={'size': 15})
    plt.yscale('log')
    plt.xlabel('Time [s]')
    #plt.ylim(1, 0.5E40)
    plt.savefig(path + 'shock_eqn' + '.png', dpi=200)
    plt.show()




















########
