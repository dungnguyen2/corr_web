import math
from scipy.integrate import quad
import constant

class condensation():
    # Saturation temperature, K, assuming =T_gas_bulk, available in PIPESIM
    # Ts=100+273
    def __init__(self, temperature, temperature_o, surface_tension_g, density_g, density_l, velocity_gas,
                 thermal_conduc_w, thermal_conduc_wall, viscosity_g, heat_capacity_g, pipe_id, thermal_conduc_g,
                 pressure_total, diffusivity_vH2O_g, thickness_wall, thickness_insulation, thermal_conduc_insulation,
                 molecular_weight_g, pressure_critical_g, pressure_critical_w, temperature_critical_g,
                 temperature_critical_w):
        # Cross sectional area (micro square m2), available in PIPESIM
        # A=0.1
        # delta: surface tension, N/m- Available in PIPESIM
        self.m_surface_tension = surface_tension_g
        # heat of vaporization, J/kg- NOT available in PIPESIM- calculated using Riazi Daubert, Predict heat of vaporization of crudes and pure components
        # Revised II
        self.m_vaporization_heat = 1.632e+6
        # Density of gas, kg/m3, available in PIPESIM
        self.m_density_g = density_g

        # Density of liquid, kg/m3, available in PIPESIM
        self.m_density_l = density_l

        # Gas bulk temperature, K, available in PIPESIM
        self.m_temperature_b_g = temperature + 273

        # Gas interface temperature, K, assumed NOT available in PIPESIM
        self.m_temperature_w_i = temperature + 2 + 273
        # Dragging coefficient
        self.m_dragging_coeff = 0.44
        # Gas velocity, m/s, available in PIPESIM
        self.m_velocity_gas = velocity_gas
        # Thermal conductivity of water, W/m/K, available in PIPESIM
        self.m_k_w = thermal_conduc_w
        # Heat transfer coefficient, W/m2/K, calculated
        # self.m_ht_coeff_i = 0.15
        # Thermal conductivity of wall, W/m/K, available in PIPESIM
        self.m_k_wall = thermal_conduc_wall
        # Thickness of the wall, m, available in PIPESIM
        self.m_thickness_wall = thickness_wall
        # Gas viscosity, cp, available in PIPESIM
        self.m_viscosity_g = viscosity_g
        # Gas heat capacity, J/kg.C, available in PIPESIM
        self.m_heat_capacity_g = heat_capacity_g
        # Pipe diameter, available in PIPESIM
        self.m_pipe_d = pipe_id
        # Thermal conductivity of the gas, available in PIPESIM
        self.m_k_g = thermal_conduc_g
        # Total pressure of the gas, barg, available in PIPESIM
        self.m_pressure_total = pressure_total
        # Diffusion number of water vapor in gas, assume the gas is pure CH4.
        self.m_diffusivity_vH2O_g = diffusivity_vH2O_g  # check!!
        # thickness of insulation, m
        self.m_thickness_insulation = thickness_insulation
        # thermal conductivity of insulation, W/m/K
        self.m_k_insulation = thermal_conduc_insulation
        # Environment temperature
        self.m_temperature_o = temperature_o + 273
        # condensation factor
        self.m_condensation_factor = 1
        # Molecular weight of gas
        self.m_molecular_weight_g = molecular_weight_g
        # Specific volume of the vapor
        self.m_volume_spec = 1 / self.m_density_g
        # critical pressure of gas, barg
        self.m_pressure_critical_g = pressure_critical_g
        # critical pressure of water
        self.m_pressure_critical_w = pressure_critical_w
        # critical temperature of gas
        self.m_temperature_critical_g = temperature_critical_g + 273
        # critical temperaure of water
        self.m_temperature_critical_w = temperature_critical_w + 273
        # MW of water
        self.m_molecular_weight_w = 18

        # self.delta_T = self.m_temperature_b_g - self.m_temperature_o
        self.delta_T = 2
        self.Ts = self.m_temperature_b_g
        self.Ti_g = (self.m_temperature_b_g + self.m_temperature_w_i) / 2

    def r_min(self, Ro_l, Ro_g, hfg, delta_T, Ts, delta):
        r_min = (Ro_l - Ro_g) / (Ro_g * Ro_l) / (hfg * delta_T) * 2 * Ts * delta
        return r_min

    def r_max(self, Ro_g, Ro_l, Cd, Ug, delta, kp):
        # calculate for the horizontal direction
        a = 2 / 3 * 3.14 * (Ro_g - Ro_l * 9.81)
        b = -1 / 8 * Cd * Ro_g * Ug ** 2
        c = 3.14 * 2 * delta
        delta_equ = (b ** 2) - (4 * a * c)
        x1, x2, x = 0, 0, 0
        if delta_equ < 0:
            print('The equation has no real solution')
        elif delta_equ == 0:
            x = (-b) / (2 * a)
            print('This equation has one solution ', x)
        else:
            x1 = (-b + (delta_equ) ** 0.5) / (2 * a)
            x2 = (-b - (delta_equ) ** 0.5) / (2 * a)
        tempo = max(x1, x2, x)
        # calculate for the vertical direction
        tempo1 = 2 * 2 * kp * delta / (Cd * Ro_g * 3.14 * Ug ** 2)
        return max(tempo, tempo1)

    # Integral calculation for Total heat
    def integrand(self, r, r_max, kH2O, hi, kW, delta, Ro_g, Ro_l, hfg, dw, dl, kl, Ti_g, Tl_o):
        A = Ti_g * (1 - 2 * delta * (Ro_l - Ro_g) / (hfg * r * Ro_g * Ro_l)) - Tl_o
        B = r / (4 * 3.14 * r ** 2 * kH2O) + 1 / (2 * 3.14 * r ** 2 * hi) + dw / (2 * 3.14 * r ** 2 * kW) + dl / (
                2 * 3.14 * r ** 2 * kl)
        C = -9 / 5 * r_max ** 0.33 / 3.14
        return A / B * C

    # Potential problem
    def Total_Heat(self, r_min, r_max, kH2O, hi, kW, delta, Ro_g, Ro_l, hfg, dw, dl, kl, Ti_g, Tl_o):
        I = quad(self.integrand, r_min, r_max,
                 args=(r_max, kH2O, hi, kW, delta, Ro_g, Ro_l, hfg, dw, dl, kl, Ti_g, Tl_o))
        return I[0]

    # calculation of heat transfer coefficient
    def Hi(self, Ro, M, hfg, Ts, v):
        return 2 * Ro / (2 - Ro) * (M / (2 * 3.14 * 8.314 * Ts)) ** 0.5 * hfg / (Ts * v)

    # calculate the Reynold number

    def Reynold(self, Ug, d, Ro_g, Nuy_g):
        return (Ug * d * Ro_g / Nuy_g)

    def Prantl(self, Cp, Nuy_g, kg):
        return (Cp * Nuy_g / kg)

    def Hg(self, Ug, d, Ro_g, Nuy_g, Cp, kg):
        return 0.023 * kg / d * self.Reynold(Ug, d, Ro_g, Nuy_g) ** 0.82 * self.Prantl(Cp, Nuy_g, kg) ** 0.4

    def Qg(self, hg, Tg_b, Ti_g):
        tempo = hg * (Tg_b - Ti_g)
        return tempo

    def psat(self, A, B, C, T):
        T_temp = T - 273
        return 10 ** (A - B / (C + T_temp)) / 760

    def gas_density(self, mw_gas, press, T, Z):
        p = press * 102e5  # convert from barg to Pa
        tempo = p * mw_gas / (Z * 8.31 * T)
        return tempo

    # return math.exp(A-B/(C+T))
    def Le(self, kg, Ro_g, Cp, Dv):
        return kg / (Ro_g * Cp * Dv)

    def Beta_Ro_gas(self, h_g, Le, hg, Cp):
        return h_g / Cp * Le ** (2 / 3)

    def diff_ideal(self, MA, MB, pcA, pcB, TcA, TcB, T):
        # MW of  Gas/22.4 density of gas at standard condition
        tempo = 2.96 * 10e-6 * (1 / MA + 1 / MB) ** 0.5 * (pcA * pcB) ** (1 / 3) / (TcA * TcB) ** (1 / 12) / (MA / 22.4)
        tempo = tempo * (T / (273 + 15)) ** (3 / 2)
        return tempo

    def mass_cal(self, temp):

        self.m_diffusivity_vH2O_g = self.diff_ideal(self.m_molecular_weight_g, self.m_molecular_weight_w,
                                                    self.m_pressure_critical_g, self.m_pressure_critical_w,
                                                    self.m_temperature_critical_g, self.m_temperature_critical_w, temp)
        self.rmin = self.r_min(self.m_density_l, self.m_density_g, self.m_vaporization_heat, self.delta_T, self.Ts,
                               self.m_surface_tension)
        # print("r_min=",r_min)
        self.rmax = self.r_max(self.m_density_g, self.m_density_l, self.m_dragging_coeff, self.m_velocity_gas,
                               self.m_surface_tension,
                               kp=1.5)
        # print("r_max=",r_max)
        self.m_ht_coeff_i = self.Hi(self.m_condensation_factor, self.m_molecular_weight_g, self.m_vaporization_heat,
                                    self.Ts,
                                    self.m_volume_spec)
        self.Total_H = self.Total_Heat(self.rmin, self.rmax, self.m_k_w, self.m_ht_coeff_i, self.m_k_wall,
                                       self.m_surface_tension,
                                       self.m_density_g, self.m_density_l, self.m_vaporization_heat,
                                       self.m_thickness_wall,
                                       self.m_thickness_insulation, self.m_k_insulation, temp, self.m_temperature_o)
        # print("Total heat=", Total_Heat)
        self.h_g = self.Hg(self.m_velocity_gas, self.m_pipe_d, self.m_density_g, self.m_viscosity_g,
                           self.m_heat_capacity_g,
                           self.m_k_g)
        self.Q_g = self.Qg(self.h_g, self.m_temperature_b_g, temp)
        # print("Qg=",Qg)
        self.xb_g = self.psat(8.07131, 1730.63, 233.426, self.m_temperature_b_g) / self.m_pressure_total  # use methane
        self.xi_g = self.psat(8.07131, 1730.63, 233.426, temp) / self.m_pressure_total
        # print("xb_gas=",xb_g)
        # print("xi_gas=", xi_g)
        self.Lewis = self.Le(self.m_k_g, self.m_density_g, self.m_heat_capacity_g, self.m_diffusivity_vH2O_g)
        mass1 = self.Beta_Ro_gas(self.h_g, self.Lewis, self.h_g, self.m_heat_capacity_g) * abs(self.xb_g - self.xi_g)
        # print("mass 1=", mass1)
        mass2 = (self.Total_H - self.Q_g) / self.m_vaporization_heat
        # print("mass 2=", mass2)
        return mass1, mass2

    def run(self):
        llimit = self.m_temperature_o
        ulimit = self.m_temperature_b_g
        epsilon = 0.0000001
        while abs(ulimit - llimit) > epsilon:
            tempo = ulimit
            tempo1, tempo2 = self.mass_cal(tempo)
            fx1 = tempo1 - tempo2
            tempo = (ulimit + llimit) / 2
            tempo1, tempo2 = self.mass_cal(tempo)
            fx2 = tempo1 - tempo2
            if (fx1 * fx2) > 0:
                ulimit = tempo
                llimit = llimit
            else:
                llimit = tempo
        self.condensation_rate = tempo1
        self.m_temperature_b_i = tempo
        return self.condensation_rate


# Calculation of Henry constant
def getHenryconst_CO2(t):
    # temperature in celsuis
    # CentiGrade As Double
    # Fahrenheit As Double

    Fahrenheit = (t - 273.15) * 9 / 5 + 32
    tempo = 14.5 / 1.00258 * 10 ** (
        -(2.27 + 0.00565 * Fahrenheit - 0.00000806 * Fahrenheit ** 2 + 0.075 * constant.ionic))
    return tempo


#
# This function returns equilibrium constant for dissociation of H2CO3
# H2CO3<=>H+ + HCO3-
#

def DisConst1st_CO2(t, pco2):
    # input temperature in Kevin
    # temperature in Fahrenheit
    # Fahrenheit As Double
    Fahrenheit = (t - 273.15) * 9 / 5 + 32

    tempo = 387.6 * 10 ** (-(6.41 - 1.594 * 10 ** (-3) * (Fahrenheit) + 8.52 * 10 ** (-6) * (Fahrenheit) ** 2 - 3.07 * 10 ** (-5) * 14.5 * pco2 - 0.4772 * constant.ionic ** 0.5 + 0.118 * constant.ionic))
    return tempo


#
# This function returns equilibrium constant for dissociation of HCO3-
# * HCO3- <=>H+ +CO32-
#

#
# This function returns equilibrium constant for dissociation of HCO3-
# * HCO3- <=>H+ +CO32-
#

def DisConst2nd_CO2(t, pco2):
    # input temperature in Kevin
    # temperature in Fahrenheit
    # Fahrenheit As Double
    Fahrenheit = (t - 273.15) * 9 / 5 + 32

    tempo = 10 ** (-(10.61 - 4.97 * 0.001 * (Fahrenheit) + 1.331 * 10 ** (-5) * (Fahrenheit) ** 2 - 2.624 * 0.00001 * pco2 * 14.5 - 1.166 * constant.ionic ** (0.5) + 0.3466 * constant.ionic))
    return tempo


#
# This function returns Henry constant for CO2
# CO2(g)<=>CO2(w)

def Fe_Cal(temperature, pco2, concFeppm, concHAcppm, pH):
    # Input temperature in Celcius
    tem_K = temperature + 273
    # calculate concentration of H+ based on pH, mol/L
    #concHplus = math.exp(pH / (-0.4343))
    concHplus =  10**(-pH)
    # concentration of CO2 in water, mol/L
    concCO2 = pco2 * getHenryconst_CO2(tem_K)
    # concentration of H2CO3 in water,mol/L
    concH2CO3 = constant.Khyd * concCO2
    # concentration of HAc in water,mol/L
    concHAc = concHAcppm / (constant.MHAc) / 1000000 * 1000
    # equilibrium constant of H2CO3 dissociation: H2CO3<=>H+ + HCO3-
    K_CO2_ca = DisConst1st_CO2(tem_K, pco2)
    # concentration of HCO3- in water, mol/L
    concHCO3 = K_CO2_ca * concH2CO3 / concHplus
    # concentratio of Fe++ in water,mol/L
    concFe = concFeppm / 55.84 / 1000000 * 1000
    # equilibrium constant of HCO3- dissociation: HCO3- <=> H+ + CO3=
    K_CO2_bi = DisConst2nd_CO2(tem_K, pco2)
    # concentration of CO3= in water, mol/L
    concCO3 = K_CO2_bi * concHCO3 / concHplus  # mol/l

    K_wa = 10**-(29.3868-0.0737549*tem_K+7.47881*1e-5*tem_K**2)
    # calculate concentration of OH- based on H+, mol/L
    concOHminus = K_wa/concHplus
    K_HAc = 10**-(6.66104 - 0.0134916*tem_K + 2.37856*1e-5*tem_K**2)
    concAc = concHAc * K_HAc / concHplus
    concFe_cal = ((concHCO3 + 2*concCO3 + concOHminus) - (concHplus))/2
    concFe_cal_Hac =  concFe_cal + concAc/2
    concFe_ppm_cal = concFe_cal * 55.84 * 1000000/ 1000
    concFe_ppm_cal_Hac = concFe_cal_Hac * 55.84 * 1000000/ 1000
    #print("Calculated Fe++ (ppm)", concFe_ppm_cal)
    #print("Fe++ as input", concFeppm)
    #print("Fe++ with consideration of Hac", concFe_ppm_cal_Hac)
    return concFe_ppm_cal_Hac

def nyborg(condensation_rate, temperature, pco2, concFeppm = 10, concHAcppm =0, pH = 5.6):
    concFe_ppm_cal_Hac = Fe_Cal(temperature, pco2, concFeppm, concHAcppm, pH)
    tempo = 0.004 * condensation_rate * concFe_ppm_cal_Hac * (12.5 - 0.09 * temperature)
    return (tempo)
