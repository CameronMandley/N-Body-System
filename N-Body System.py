#!/usr/bin/env python
# coding: utf-8

# In[137]:


def N_body_system(initial_states, masses, G, tf):
    """
    Parameters: - initial_states: A vector containing 2 vectors, a vector of N initial positions 
                  and a vector of N initial momenta both of dim 3
                - masses: masses of the respective particles with initial states initial_states
                - G: Gravitational constant
                - tf: Finishing time of the simulation
                
    Returns: The positions and momenta of N particles at time tf
    """
    # Extracting the N initial position and momentum vectors
    q = np.array([initial_states[0]])[0]
    p = np.array([initial_states[1]])[0]
    N = len(p)
    # Initializing the forces for each particle and setting them to 0
    F_q = np.array([[0.0, 0.0, 0.0]]*N)
    for ti in range(tf):
        for i in range(N):
            for j in range(i+1, N-1):
                # Updating the forces for particle i and j
                delta = G*masses[i]*masses[j]
                v_ij = q[i] - q[j]
                q_ij = np.sqrt(v_ij@v_ij)
                F_ix = q[i][0]*delta*q_ij**(-3)
                F_iy = q[i][1]*delta*q_ij**(-3)
                F_iz = q[i][2]*delta*q_ij**(-3)
                F_jx = -F_ix
                F_jy = -F_iy
                F_jz = -F_iz
                F_q[i] = np.array([F_q[i][0] + F_ix, F_q[i][1] + F_iy, F_q[i][2] + F_iz])
                F_q[j] = np.array([F_q[j][0] + F_jx, F_q[j][1] + F_jy, F_q[j][2] + F_jz])
            # Updating the positions and momenta for particles i and j
            p[i] += F_q[i]*ti
            q[i] += p[i]*ti/masses
            p[j] += F_q[j]*ti
            q[j] += p[j]*ti/masses
        
    return [q, p]


# In[138]:


get_ipython().run_line_magic('matplotlib', 'inline')
import numpy as np
import matplotlib.pyplot as plt
fig1 = plt.figure()
ax1 = plt.axes(projection='3d')
fig2 = plt.figure()
ax2 = plt.axes(projection='3d')

positions = []
momenta = []
q_0 = np.array([[random.uniform(0, 10), random.uniform(0, 10), random.uniform(0, 10)], [random.uniform(0, 10), random.uniform(0, 10), random.uniform(0, 10)], [random.uniform(0, 10), random.uniform(0, 10), random.uniform(0, 10)]])
p_0 = np.array([[random.uniform(0, 10), random.uniform(0, 10), random.uniform(0, 10)], [random.uniform(0, 10), random.uniform(0, 10), random.uniform(0, 10)], [random.uniform(0, 10), random.uniform(0, 10), random.uniform(0, 10)]])
masses = [random.uniform(0, 10), random.uniform(0, 10), random.uniform(0, 10)]
for i in range(20):
    state = N_body_system([q_0, p_0], masses, 6.67408*10**(-11), i)
    positions.append(state[0])
    momenta.append(state[1])
    print("Time: " + str(i))
    print("Positions = " + str(state[0]))
    print("Momenta = " + str(state[1]))
    print()

zline_pos = [positions[i][j][2] for i in range(len(positions)) for j in range(2)]
xline_pos = [positions[i][j][0] for i in range(len(positions)) for j in range(2)]
yline_pos = [positions[i][j][1] for i in range(len(positions)) for j in range(2)]
ax1.plot3D(xline_pos, yline_pos, zline_pos, 'gray')
ax1.set_xlabel('position: x')
ax1.set_ylabel('position: y')
ax1.set_zlabel('position: z');

zline_mom = [momenta[i][j][2] for i in range(len(momenta)) for j in range(2)]
xline_mom = [momenta[i][j][0] for i in range(len(momenta)) for j in range(2)]
yline_mom = [momenta[i][j][1] for i in range(len(momenta)) for j in range(2)]
ax2.plot3D(xline_mom, yline_mom, zline_mom, 'gray')
ax2.set_xlabel('momentum: x')
ax2.set_ylabel('momentum: y')
ax2.set_zlabel('momentum: z');

