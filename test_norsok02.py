import norsokm506_01 as ns
import sys
# This is written to calculate the corrosion rate at different level of CO2
# argument: CO2 fraction, temperature in C, pressure in barg.

if __name__ == "__main__":
    if len(sys.argv) >= 2:
        CO2fraction = float(sys.argv[1])
    else:
        CO2fraction = 0.2
    if len(sys.argv) == 3:
        temp = float(sys.argv[2]) # in degree C
    else:
        temp = 120
    if len(sys.argv) == 4:
        press = float(sys.argv[3]) # in barg
    else:
        press = 37.6

mass_g = 100/35.4*1762*24 #kg/hr mass flow of 100 mmscfd of gas of MW at 24
mass_l = 0.1
vol_g = 100/35.2/24*1e6/press #m3/hr #100 mmscfd at KL
vol_l = 1 #m3/hr liquid rate is very insignificant
vis_l = 1.22 #viscosity, centi point
vis_g = 0.012 #viscosity, centi point
roughness = 0.00005
dia = 16*0.025 #diameter of 16 inches pipe
density_g = mass_g / vol_g
density_l = 800


v_sg = 100/35.2/24/(0.76*dia**2)*1e6/press #m/s #velocity of gas at 100 mmscfd at 37 barg
v_sl = 0.1 #m/s # no liquid, assume a fake velocity

bicarbonate = 2003
ionstrength = 39.66

holdup = 1 #percentage of liquid occupied in the pipe.
print("                                     ")
print("------------------------------------")
#print("Temperature=", temp)
kt=ns.Kt(temp)
#print(f"KT = {kt:.2f}")
FugCO2=ns.FugacityofCO2(CO2fraction, press, temp)
print(f"Fugacity of CO2 = {FugCO2:.2f}")
shearstress=ns.Shearstress(v_sg, v_sl , mass_g, mass_l,vol_g,vol_l,holdup,vis_g,vis_l,roughness,dia, density_g, density_l)
print(f"Shear stress = {shearstress:.2f}")
ph=ns.pHCalculator(temp, press, CO2fraction*press, bicarbonate, ionstrength, 2)
print(f"PH = {ph:.2f}")
fph=ns.fpH_Cal(temp,float(ph))
#print(f"f_pH = {fph}")
corr=ns.Cal_Norsok(CO2fraction, press, temp, v_sg, v_sl , mass_g, mass_l,vol_g,vol_l,holdup,vis_g,vis_l,roughness,dia, fph,bicarbonate, ionstrength, 2, density_g, density_l)
print("density of gas = ", density_g)
print(f"Corrosion rate {corr:.2f} mm/year")
print("------------------------------------")
print("                                     ")