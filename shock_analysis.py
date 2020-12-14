# Functions for finding shocks in output of hydro code
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

def downstream_density(finalstep, path=''):

    import numpy as np

    u = read_output(finalstep, path)

    time = u[:,0,7]             # [s]
    density = u[:,:,1]        # [g/cm^3]

    shock_indices = extract_shocks(finalstep, path)[1]

    downstream_densities = []
    for i in range(len(time)):
        j = shock_indices[i]
        downstream_density_ = density[i,j-1]
        downstream_densities.append(downstream_density_)

    return downstream_densities

def downstream_soundspeed(finalstep, path=''):

    import numpy as np

    u = read_output(finalstep, path)

    time = u[:,0,7]             # [s]
    c_s = u[:,:,9]            # [cm/s]

    shock_indices = extract_shocks(finalstep, path)[1]

    downstream_cs = []
    for i in range(len(time)):
        j = shock_indices[i]
        downstream_cs_i = c_s[i,j-1]
        downstream_cs.append(downstream_cs_i)

    return downstream_cs

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
        t = plt.vlines(shock_positions[i], 0, -1*v_jumps[i], color='r')
        plt.xlabel(r'Radius [cm]')
        plt.ylabel(r'Velcity [cm/s]')
        plt.xscale('log')
        #plt.legend(t, [f'time = {time[i]} s'])
        camera.snap()

    animation = camera.animate()
    animation.save(path + 'shocks_animation.mp4')
