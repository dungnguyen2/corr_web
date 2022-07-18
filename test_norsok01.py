import norsokm506_01 as ns


v_sg = 9
v_sl = 1
mass_g = 234.5
mass_l = 542.3
vol_g = 14.8
vol_l = 637
vis_l = 1.4
vis_g = 0.03
roughness = 0.00005
dia = 0.475
temp = 64
press = 37.6
bicarbonate = 2003
ionstrength = 39.66
CO2fraction = 0.2
holdup=10
print("                                     ")
print("------------------------------------")
print("Temperature=", temp)
kt=ns.Kt(temp)
print("KT=", kt)
FugCO2=ns.FugacityofCO2(CO2fraction, press, temp)
print(f"Fugacity of CO2 = {FugCO2:.2f}")
shearstress=ns.Shearstress(v_sg, v_sl , mass_g, mass_l,vol_g,vol_l,holdup,vis_g,vis_l,roughness,dia)
print(f"Shear stress = {shearstress:.2f}")
ph=ns.pHCalculator(temp, press, CO2fraction*press, bicarbonate, ionstrength, 2)
print(f"PH = {ph:.2f}")
fph=ns.fpH_Cal(temp,float(ph))
print("f_pH=",fph)
corr=ns.Cal_Norsok(CO2fraction, press, temp, v_sg, v_sl , mass_g, mass_l,vol_g,vol_l,holdup,vis_g,vis_l,roughness,dia, fph,bicarbonate, ionstrength, 2)
print(f"Corrosion rate {corr:.2f} mm/year")
print("------------------------------------")
print("                                     ")