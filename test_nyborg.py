import tol005 as tol

# delta: surface tension, N/m- Available in PIPESIM
surface_tension_g = 0.1

#density of gas
density_g = 18
#density of liquid
density_l = 700
# Velocity of gas
velocity_gas = 5
# Thermal conductivity of water, W/m/K, available in PIPESIM
thermal_conduc_w = 0.6
# Thermal conductivity of wall, W/m/K, available in PIPESIM
thermal_conduc_wall = 65
# Thickness of the wall, m, available in PIPESIM
thickness_wall = 0.025
# Gas viscosity, cp, available in PIPESIM
viscosity_g = 0.013
# Gas heat capacity, J/kg.C, available in PIPESIM
heat_capacity_g = 2.3
# Pipe diameter, available in PIPESIM
pipe_id = 0.7
# Thermal conductivity of the gas, available in PIPESIM
thermal_conduc_g = 0.07
# Total pressure of the gas, barg, available in PIPESIM
pressure_total = 3
# Diffusion number of water vapor in gas, assume the gas is pure CH4.
diffusivity_vH2O_g = 3.56e-5  # check!!
# thickness of insulation, m
thickness_insulation = 0.030
# thermal conductivity of insulation, W/m/K
thermal_conduc_insulation = 0.32
# Environment temperature
temperature_o = 25
# Fluid temperature
temperature = 100
# Molecular weight of gas
molecular_weight_g = 23
# critical pressure of gas, barg
pressure_critical_g = 46
# critical pressure of water
pressure_critical_w = 220
# critical temperature of gas
temperature_critical_g = -34
# critical temperaure of water
temperature_critical_w = 374
CO2fraction = 0.2
pco2 = CO2fraction * pressure_total

cd = tol.condensation(temperature, temperature_o, surface_tension_g, density_g, density_l, velocity_gas, thermal_conduc_w, thermal_conduc_wall, viscosity_g, heat_capacity_g, pipe_id, thermal_conduc_g, pressure_total, diffusivity_vH2O_g, thickness_wall, thickness_insulation, thermal_conduc_insulation, molecular_weight_g, pressure_critical_g, pressure_critical_w, temperature_critical_g, temperature_critical_w)
cd.run()
nyborg = tol.nyborg(cd.condensation_rate, temperature, pco2)
print("Corrosion rate={:0.3f} mm/y".format(nyborg))
print("Interface temperature={:0.3f}".format(cd.m_temperature_b_i-273))
