"""RocketSimulation (c) by Roi Dvir
RocketSimulation is licensed under a Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
You should have received a copy of the license along with this work. If not, see <http://creativecommons.org/licenses/by-nc-nd/4.0/>."""


# coding: utf-8
"""
Graphs for the work
"""
from rocket_project import *  # import everything from rocket_project.py
from matplotlib import pyplot as plt
import numpy as np
# %% Rho(z) - fig 3
plt.figure()
z = np.arange(0, 20e3)
plt.plot(z/1e3, list(map(Rho, z)))
plt.grid(True)
plt.xlabel("Altitude [km]")
plt.ylabel("Air Density [kg/m^3]")
# %% G(z) - not in the work
plt.figure()
z = np.arange(0, 1e7, 1e4)
plt.plot(z/1e3, list(map(Gravity, z)))
plt.grid(True)
plt.xlabel("Altitude [km]")
plt.ylabel("Gravity [N]")
# %% Comparing Angles - Figure 4
plt.figure()
theta0_range = [22.5, 45, 45+22.5]
for theta0 in theta0_range:
    parameters = set_parameters(theta_rocket_degree=0, gas_velocity=0, m0=1,
                                Cd=0, Cl=0, Sd=0, Sl=0, dmdt=0)
    values = set_initial_values(x=0, z=0, velocity=10,
                                theta_velocity_degree=theta0,
                                fuel_mass=0, parameters=parameters)
    results = Euler_Cromer(values, parameters, t_init=0, t_final=3,
                           dt=1e-4)
    plt.plot(results.x, results.z, 'o', label=str(theta0))
    t = results.index.values
    plt.plot(values['vx'] * t, values['vz'] * t - 0.5 * 9.80665 * (t**2),'k-',
             linewidth=2,
             label=str(theta0) + " theory")
plt.legend()
plt.grid(True)
plt.xlabel("X [m]")
plt.ylabel("Z [m]")
# %% Comparing dt - not in the work
plt.figure()
dt_range = np.logspace(-1, -3, 3)
for dt in dt_range:
    parameters = set_parameters(theta_rocket_degree=0, gas_velocity=10, m0=1,
                                Cd=0, Cl=0, Sd=0, Sl=0, dmdt=0)
    values = set_initial_values(x=0, z=0, velocity=10, theta_velocity_degree=45,
                                fuel_mass=0, parameters=parameters)
    results = Euler_Cromer(values, parameters, t_init=0, t_final=3000,
                           dt=dt)
    plt.figure(num='rocket simulation')
    plt.subplot(2, 2, 1)
    results.x.plot()
    plt.subplot(2, 2, 2)
    results.z.plot()
    plt.subplot(2, 1, 2)
    plt.plot(results.x, results.z)
plt.legend(dt_range)
# %% Tsiolkovsky rocket equation - Figure 5
plt.figure()
res = []
gas_velocity = 10
fuel_mass = 1
body_mass = 1
dmdt_range = np.arange(2, 100)
for dmdt in dmdt_range:
    parameters = set_parameters(theta_rocket_degree=90,
                                gas_velocity=gas_velocity, m0=1,
                                Cd=0, Cl=0, Sd=0, Sl=0, dmdt=dmdt)
    values = set_initial_values(x=0, z=10, velocity=0, theta_velocity_degree=0,
                                fuel_mass=1, parameters=parameters)
    results = Euler_Cromer(values, parameters, t_init=0,
                           t_final=2, dt=1e-4)
    res.append(results.vz.max())
time_till_end_of_fuel = fuel_mass/dmdt_range
maximal_velocity = gas_velocity * np.log(2/1) - 9.8 * time_till_end_of_fuel
plt.plot(dmdt_range, res, 'o', label="simulation")
plt.plot(dmdt_range, maximal_velocity, label="theory")
plt.xlabel("dm/dt [kg/sec]")
plt.ylabel("maximal velocity [m/s]")
plt.grid()
plt.legend()
# %% escape velocity - Figure 6
plt.figure()
res = []
for v in [10e3, 50e3]:
    parameters = set_parameters(theta_rocket_degree=0, gas_velocity=0, m0=1,
                                Cd=0, Cl=0, Sd=0, Sl=0, dmdt=0)
    values = set_initial_values(x=0, z=1, velocity=v, theta_velocity_degree=90,
                                fuel_mass=0, parameters=parameters)
    results = Euler_Cromer(values, parameters, t_init=0,
                           t_final=20000, dt=1e-2)
    res.append(results)
plt.subplot(1, 2, 1)
plt.semilogy(res[0].z)
plt.semilogy(res[1].z)
plt.grid(True)
plt.legend(["v=10,000 m/s", "v=50,000 m/s"])
plt.subplot(1, 2, 2)
plt.plot(res[0].vz)
plt.plot(res[1].vz)
plt.grid(True)
plt.legend(["v=10,000 m/s", "v=50,000 m/s"])
# %% ballistic throw with drag - Figure 7
plt.figure()
res = []
cd_range = np.linspace(0, 1e-1, 3)
for cd in cd_range:
    parameters = set_parameters(theta_rocket_degree=0, gas_velocity=0, m0=1,
                                Cd=cd, Cl=0, Sd=1, Sl=0, dmdt=0)
    values = set_initial_values(x=0, z=0, velocity=10,
                                theta_velocity_degree=45,
                                fuel_mass=0, parameters=parameters)
    results = Euler_Cromer(values, parameters, t_init=0,
                           t_final=5, dt=1e-4)
    res.append(results)
for cd, r in zip(cd_range, res):
    plt.plot(r.x, r.z, label=f"Cd={cd}")
plt.grid(True)
plt.legend()
plt.xlabel("X [m]")
plt.ylabel("Z [m]")
plt.figure()
for cd, r in zip(cd_range, res):
    plt.plot(r.index, r.z, label=f"Cd={cd}")
plt.grid(True)
plt.legend()
plt.xlabel("time [sec]")
plt.ylabel("Z [m]")
# %% ballistic throw with lift range vs Cd - Figure 8
plt.figure()
res = []
cd_range = np.arange(0.25, 2, 0.25)
for cd in cd_range:
    parameters = set_parameters(theta_rocket_degree=0, gas_velocity=0, m0=1,
                                Cd=cd, Cl=0, Sd=1, Sl=0, dmdt=0)
    values = set_initial_values(x=0, z=0, velocity=136,
                                theta_velocity_degree=60,
                                fuel_mass=0, parameters=parameters)
    results = Euler_Cromer(values, parameters, t_init=0,
                           t_final=5, dt=1e-3)
    res.append(results)
for cd, r in zip(cd_range, res):
    plt.plot(r.x, r.z, label=f"Cd={cd}")
plt.grid(True)
plt.legend()
plt.xlabel("X [m]")
plt.ylabel("Z [m]")
# %% ballistic throw with lift - Figure 9
plt.figure()
res = []
cl_range = [0, 0.05, 0.1]
for cl in cl_range:
    parameters = set_parameters(theta_rocket_degree=0, gas_velocity=0, m0=1,
                                Cd=0, Cl=cl, Sd=0, Sl=1, dmdt=0)
    values = set_initial_values(x=0, z=0, velocity=10,
                                theta_velocity_degree=45,
                                fuel_mass=0, parameters=parameters)
    results = Euler_Cromer(values, parameters, t_init=0,
                           t_final=5, dt=1e-3)
    res.append(results)
for cl, r in zip(cl_range, res):
    plt.plot(r.x, r.z, label=f"Cl={cl}")
plt.grid(True)
plt.legend()
plt.xlabel("X [m]")
plt.ylabel("Z [m]")
plt.figure()
for cl, r in zip(cd_range, res):
    plt.plot(r.index, r.z, label=f"Cl={cl}")
plt.grid(True)
plt.legend()
plt.xlabel("time [sec]")
plt.ylabel("Z [m]")
# %% ballistic throw with drag range vs Cd - Figure 10
plt.figure()
res = []
cd_range = np.arange(0, 1, 0.01)
plt.close('all')
for cd in cd_range:
    parameters = set_parameters(theta_rocket_degree=0, gas_velocity=0, m0=1,
                                Cd=cd, Cl=0, Sd=1, Sl=0, dmdt=0)
    values = set_initial_values(x=0, z=0, velocity=10,
                                theta_velocity_degree=45,
                                fuel_mass=0, parameters=parameters)
    results = Euler_Cromer(values, parameters, t_init=0,
                           t_final=5, dt=1e-3)
    res.append(results.x.tolist()[-1])
plt.plot(cd_range, res)
plt.xlabel("Cd ")
plt.ylabel("X [m]")
plt.grid(True)
# %% ballistic throw with lift range vs Cl - Figure 11
plt.figure()
res = []
cl_range = np.linspace(0, 1, 100)
for cl in cl_range:
    parameters = set_parameters(theta_rocket_degree=0, gas_velocity=0, m0=1,
                                Cd=0, Cl=cl, Sd=0, Sl=1, dmdt=0)
    values = set_initial_values(x=0, z=0, velocity=10,
                                theta_velocity_degree=45,
                                fuel_mass=0, parameters=parameters)
    results = Euler_Cromer(values, parameters, t_init=0,
                           t_final=5, dt=1e-4)
    res.append(results.tail(1))
plt.plot(cl_range, [r.x.tolist()[0] for r in res])
plt.xlabel("Cl ")
plt.ylabel("X [m]")
plt.grid(True)
plt.figure()
plt.plot(results.index, results.x, label="X")
plt.plot(results.index, results.z, label="X")
plt.xlabel("Cl ")
plt.ylabel("meter")
plt.grid(True)
# %% aircraft - Figure 12
plt.figure()
res = []
v_range = np.linspace(200, 300, 3) * 1000 / 3600
for v in v_range:
    parameters = set_parameters(theta_rocket_degree=0, gas_velocity=0, m0=78e3,
                                Cd=0, Cl=2, Sd=0, Sl=122.6, dmdt=0)
    values = set_initial_values(x=0, z=5, velocity=v,
                                theta_velocity_degree=0,
                                fuel_mass=0, parameters=parameters)
    results = Euler_Cromer(values, parameters, t_init=0, t_final=2, dt=1e-3)
    res.append(results)
for v, r in zip(v_range, res):
    plt.plot(r.index, r.z, label=f"v = {v*3600/1000}  km/h")
plt.xlim([0, 2])
plt.ylim([0, 30])
plt.grid(True)
plt.legend()
plt.xlabel("X [m]")
plt.ylabel("Z [m]")
# %% aircraft - Figure 14 - takeoff velocity vs. lift coefficient
"""
0.5 Cl Sl rho(0) v^2 = m g  ==> v = sqrt( 2g/rho(0) ) * sqrt( m / (Cl Sl) )
"""
plt.clf()
res = []
m0 = 1
sl = 1
z0 = 5
cl_range = np.arange(1, 10, .1)
for cl in cl_range:
    v_theory = np.sqrt( 2*Gravity(0)/Rho(0) ) * np.sqrt( m0 / (cl*sl) )
    v_range = np.arange(v_theory*0.98, 1.02*v_theory, v_theory/100)
    for v in v_range:
        parameters = set_parameters(theta_rocket_degree=0, gas_velocity=0, m0=m0,
                                    Cd=0, Cl=cl, Sd=0, Sl=sl, dmdt=0)
        values = set_initial_values(x=0, z=z0, velocity=v,
                                    theta_velocity_degree=0,
                                    fuel_mass=0, parameters=parameters)
        results = Euler_Cromer(values, parameters, t_init=0, t_final=1e-3, dt=1e-4)
        if results.z.tolist()[-1] > z0:  # takeoff
            res.append(v)
            break # no need to check larger v values
plt.plot(cl_range, res, 'o', label="simulation")
plt.plot(cl_range, np.sqrt( 2*Gravity(0)/Rho(0) ) * np.sqrt( m0/(cl_range*sl)),
         label="theory")
plt.grid(True)
plt.xlabel("Cl * Sl")
plt.ylabel("takeoff velocity [m/s]")
plt.legend()
# %% aircraft - Figure 15 - takeoff velocity vs. lift mass
"""
0.5 Cl Sl rho(0) v^2 = m g  ==> v = sqrt( 2g/rho(0) ) * sqrt( m / (Cl Sl) )
"""
plt.clf()
res = []
m0 = 1
cl = 1
sl = 1
z0 = 5
m_range = np.arange(1, 10, .1)
for m0 in m_range:
    v_theory = np.sqrt( 2*Gravity(0)/Rho(0) ) * np.sqrt( m0 / (cl*sl) )
    v_range = np.arange(v_theory*0.98, 1.1*v_theory, v_theory/100)
    for v in v_range:
        parameters = set_parameters(theta_rocket_degree=0, gas_velocity=0, m0=m0,
                                    Cd=0, Cl=cl, Sd=0, Sl=sl, dmdt=0)
        values = set_initial_values(x=0, z=z0, velocity=v,
                                    theta_velocity_degree=0,
                                    fuel_mass=0, parameters=parameters)
        results = Euler_Cromer(values, parameters, t_init=0, t_final=1e-3, dt=1e-4)
        if results.z.tolist()[-1] > z0:  # takeoff
            res.append(v)
            break # no need to check larger v values
plt.plot(m_range, res, 'o', label="simulation")
plt.plot(m_range, np.sqrt( 2*Gravity(0)/Rho(0) ) * np.sqrt( m_range/(cl*sl)),
         label="theory")
plt.grid(True)
plt.xlabel("mass [kg]")
plt.ylabel("takeoff velocity [m/s]")
plt.legend()
# %% Comparing dmdt
plt.figure()
dmdt_range = np.logspace(0, 1, 2)
for dmdt in dmdt_range:
    parameters = set_parameters(theta_rocket_degree=90, gas_velocity=10, m0=1,
                                Cd=0, Cl=0, Sd=0, Sl=0, dmdt=dmdt)
    values = set_initial_values(x=0, z=10, velocity=0, theta_velocity_degree=0,
                                fuel_mass=1, parameters=parameters)
    results = Euler_Cromer(values, parameters, t_init=0, t_final=50,
                           dt=1e-4)
    print(results.tail(1))
    plt.figure(num='rocket simulation')
    plt.subplot(2, 2, 1)
    results.x.plot()
    plt.subplot(2, 2, 2)
    results.z.plot()
    plt.subplot(2, 1, 2)
    plt.plot(results.x, results.z)
plt.legend(dmdt_range)
plt.figure()
v = np.sqrt(results.vx**2 + results.vz**2)
v.plot()
results.vx.plot()
results.vz.plot()