"""RocketSimulation (c) by Roi Dvir
RocketSimulation is licensed under a Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
You should have received a copy of the license along with this work. If not, see <http://creativecommons.org/licenses/by-nc-nd/4.0/>."""

# coding: utf-8
import numpy as np
#

def set_initial_values(x, z, velocity, theta_velocity_degree, fuel_mass,
                       parameters):
    """
    function for initalizing all the variables of the simulations
    x, z: the initial position of the rocket in meters
    velocity: the magnitude of the rocket velocity at t=0 in m/s
    theta_velocity_degree: the degree between the rocket velocity
        and the x axis
    fuel_mass: fuel mass at t=0 in kg
    parameters: dict returned from set_parameters function

    returns: a dictionary of the rocket values at t=0
    """
    values = {
        "x": x,
        "z": z,
        "vx": velocity*np.cos(np.radians(theta_velocity_degree)),
        "vz": velocity*np.sin(np.radians(theta_velocity_degree)),
        "mass": fuel_mass+parameters["m0"]}
    print(values)
    return values


def set_parameters(theta_rocket_degree, gas_velocity, m0, Cd, Cl, Sd, Sl,
                   dmdt):
    """
    function for initalizing all the constant parameters of the simulations
    theta_rocket_degree: the degree between the rocket and the x axis.
        This value impact the direction of thrust
    gas_velocity: the magnitude of the rocket gas velocity in m/s
    m0: self mass of the rocket in kg. this is the mass after all fuel is gone
    Cd, Cl: drag and lift coefficients, respectivly
    Sd, Sl: drag and lift surface size in m^2
    dmdt : mass loss coefficient (referd to as alpha in the work)

    returns: a dictionary of the rocket constant parameters
    """
    parameters = {
        "theta_rocket": np.radians(theta_rocket_degree),
        "gas_velocity": gas_velocity,
        "m0": m0,
        "Cd": Cd,
        "Cl": Cl,
        "Sd": Sd,
        "Sl": Sl,
        "dmdt": dmdt}
    print(parameters)
    return parameters


def Rho(z):
    """
    air density at altitude of z meters based on
        https://en.wikipedia.org/wiki/Density_of_air#Altitude
    z: altitude in meters

    returns: air density in kg/m^3
    """
    p0 = 101.325e3  # sea level standard atmospheric pressure, Pa
    T0 = 288.15     # sea level standard temperature, K
    g = 9.80665     # earth-surface gravitational acceleration, m/s²
    L = 0.0065      # temperature lapse rate, K/m
    R = 8.31447     # ideal (universal) gas constant, J/(mol·K)
    M = 0.0289644   # molar mass of dry air, kg/mol
    T = T0-L*z
    if T < 0:
        return(0)
    p = p0 * (T/T0)**(g*M/(R*L))
    rho = p*M/(R*T)
    return(rho)


def Gravity(z, me_kg=5.972e24, re_m=6371e3):
    """
    Gravity at altitude z meters above sea level
    z: altitude in meters
    returns: gravity aceleration in m/sec^2
    """
    from scipy.constants import G
    return(G*me_kg/(re_m+z)**2)


def Drag(values, parameters):
    """
    calculate the drag force = 0.5*Cd*Sd*(v^2)
    values and parameters - the dictionaries of the current rocket state
    returns drag force in N
    """
    return(0.5 * parameters["Cd"] * Rho(values["z"]) * parameters["Sd"] *
           (values["vx"]**2 + values["vz"]**2))


def Lift(values, parameters):
    """
    calculate the lift force = 0.5*Cl*Sl*rho*(v^2)
    values and parameters - the dictionaries of the current rocket state
    returns lift force in N
    """
    return(0.5 * parameters["Cl"] * Rho(values["z"]) * parameters["Sl"] *
           (values["vx"]**2 + values["vz"]**2))


def Thrust(t, values, parameters):
    """
    calculate the thrus force = alpha * (v_gas)
    values and parameters - the dictionaries of the current rocket state
    returns thrust force in N
    """
    return dm_dt(t, values) * parameters["gas_velocity"]


def theta_velocity(vx, vz):
    # returns the velocity degree using its x and z elements
    return(np.arctan2(vz, vx))


# % EQUATIONS
def dm_dt(t, values, parameters):
    # equation for the mass change
    # dm/dt = alpha (if there is fuel to burn)
    # dm/dt = 0 (if there is no more fuel to burn)
    m = values["mass"]
    if m > parameters["m0"]:  # rocket mass contains fuel
        return(-parameters["dmdt"])
    else:  # no change in mass if we are at the rocket self mass
        return(0)


def dvz_dt(t, values, parameters):
    # acellaration in z direction - see the related work for more details
    # returns the rhs of the equation dvz/dt = az
    z = values["z"]
    vx = values["vx"]
    vz = values["vz"]
    m = values["mass"]
    theta_rocket = parameters["theta_rocket"]
    gas_velocity = parameters["gas_velocity"]
    # the froces in Z direction
    F_thrust_z = \
        -dm_dt(t, values, parameters) * gas_velocity * np.sin(theta_rocket)
    F_drag_z = Drag(values, parameters)*np.sin(theta_velocity(-vx, -vz))
    # F_lift is perpendicular to v
    F_lift_z = -Lift(values, parameters)*np.cos(theta_velocity(-vx, -vz))
    return ((F_thrust_z + F_drag_z + F_lift_z)/m-Gravity(z))


def dvx_dt(t, values, parameters):
    # acellaration in x direction - see the related work for more details
    # returns the rhs of the equation dvx/dt = ax
    vx = values["vx"]
    vz = values["vz"]
    m = values["mass"]
    theta_rocket = parameters["theta_rocket"]
    gas_velocity = parameters["gas_velocity"]
    # the froces in X direction
    F_thrust_x = \
        -dm_dt(t, values, parameters) * gas_velocity * np.cos(theta_rocket)
    # drag is opposite to velocity
    F_drag_x = Drag(values, parameters) * np.cos(theta_velocity(-vx, -vz))
    # F lift is perpendicular to v
    F_lift_x = Lift(values, parameters) * np.sin(theta_velocity(-vx, -vz))
    return ((F_thrust_x + F_drag_x + F_lift_x)/m)


def dz_dt(t, values, parameters):
    # an euqation for the position using the relation dz/dt = vz
    # returns the rhs of the equation dz/dt = vz
    vz = values["vz"]
    return(vz)


def dx_dt(t, values, parameters):
    # an euqation for the position using the relation dx/dt = vx
    # returns the rhs of the equation dx/dt = vx
    vx = values["vx"]
    return(vx)


# % Euler Cromer Solver
# equation to solve. see Euler_Cromer_Step function for more details
rhs_dict = {"mass": dm_dt, "x": dx_dt, "z": dz_dt, "vx": dvx_dt, "vz": dvz_dt}


def Euler_Cromer_Step(rhs_dict, values, parameters, t0, dt=1e-3):
    new_values = values.copy()
    """
     rhs_dict holds the equations to solve in the form {variable: eqaution}
     the update step according to euler is
     new_value = old_value + dt * (equation evaluated at this time)
     since we update each field in a seperate step (sequential update)
     we get the euler-cromer method. for example:
     key1 is update
     key2 is updated using the values of the updated key1
     which i exactly euler-cromer method.
     (in regular euler key2 uses old key1)
    """
    for key in rhs_dict.keys():
        new_values[key] += rhs_dict[key](t0, values, parameters) * dt
    return new_values


def Euler_Cromer(values, parameters, t_init, t_final, dt=1e-3,
                 rhs_dict=rhs_dict):
    import pandas as pd     # for returning a dataframe object
    from tqdm import tqdm   # for displaying a progress bar during simulation
    t_values = np.arange(t_init, t_final, dt)  # time vector to simulate
    # we collect the results at each time step into the list results
    results = []
    results.append(values.copy())
    for t in tqdm(t_values[1:]):
        # run the euler-cromer step onec to go from t to t+dt
        results.append(Euler_Cromer_Step(rhs_dict, values=results[-1],
                                         parameters=parameters, t0=t, dt=dt))
        if results[-1]["z"] < 0:  # exit condition if rocket under sea level
            results.pop()
            break
    # return a DataFrame object for better reading of the results
    df = pd.DataFrame(results, index=t_values[:len(results)])
    df.index.name = "time"
    return df
