#!python3
import norsokm506_01 as ns
from flask import Flask, render_template, request

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def index():
    ns_rate = ''
    if request.method == 'POST' and 'temp' in request.form:
        temp = float(request.form.get('temp'))
        press = float(request.form.get('press'))
        press = 35
        CO2fraction = float(request.form.get('CO2fraction'))        
        holdup = float(request.form.get('holdup'))
        gasrate = float(request.form.get('gasrate'))
        mw  = float(request.form.get('mw'))
        mass_g = gasrate/press*1762*mw
        vol_l = float(request.form.get('vol_l'))
        density_l = float(request.form.get('density_l'))
        mass_l = vol_l * density_l
        vis_l = float(request.form.get('vis_l'))
        vis_g = float(request.form.get('vis_g'))
        dia = float(request.form.get('dia'))
        bicarbonate = float(request.form.get('bicarbonate'))
        ionstrength = float(request.form.get('ionstrength'))
        roughness = float(request.form.get('roughness'))
        ns_rate = ns_corr(temp, press, CO2fraction, holdup, mass_g, mw, vol_l, density_l, gasrate, mass_l, vis_l, vis_g, dia, bicarbonate, ionstrength, roughness)
    return render_template("bmi_calc.html",
	                        ns_rate=ns_rate)

def ns_corr(temp, press, CO2fraction, holdup, mass_g, mw, vol_l, density_l, gasrate, mass_l, vis_l, vis_g, dia, bicarbonate, ionstrength, roughness):
    
    #mass_g = 100/35.4*1762*24 #kg/hr mass flow of 100 mmscfd of gas of MW at 24
    #mass_l = 0.1
    vol_g = gasrate/35.2/mw*1e6/press #m3/hr #100 mmscfd at KL
    #vol_l = 1 #m3/hr liquid rate is very insignificant
    #vis_l = 1.22 #viscosity, centi point
    #vis_g = 0.012 #viscosity, centi point
    #roughness = 0.00005
    #dia = 16*0.025 #diameter of 16 inches pipe
    density_g = mass_g / vol_g
    #density_l = 800
    
    
    v_sg = gasrate/35.2/mw/(0.76*dia**2)*1e6/press #m/s #velocity of gas at 100 mmscfd at 37 barg
    #v_sl = 0.1 #m/s # no liquid, assume a fake velocity
    v_sl = vol_l/(0.76*dia**2)*3600
    
    #bicarbonate = 2003
    #ionstrength = 39.66
    
    #holdup = 1 #percentage of liquid occupied in the pipe.
    
    #CO2fraction = 0.2
    ph=ns.pHCalculator(temp, press, CO2fraction*press, bicarbonate, ionstrength, 2)
    fph=ns.fpH_Cal(temp,float(ph))
    tempo = ns.Cal_Norsok(CO2fraction, press, temp, v_sg, v_sl , mass_g, mass_l,vol_g,vol_l,holdup,vis_g,vis_l,roughness,dia, fph,bicarbonate, ionstrength, 2, density_g, density_l)
    
    return tempo

if __name__ == '__main__':
    app.debug = True
    app.run()
    app.run(debug=True)