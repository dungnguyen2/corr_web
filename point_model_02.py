# ===========================================================================
#    pycorp - Corrosion Prediction Software
#    Copyright (C) 2008 Srdjan Nesic, Dusan Sormaz, Hui Li, Jing Huang.
#	 Original name: FreeCorp.
#	 Developed in Python by Dung Nguyen
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#    for information on how to use the program visit
#    <http://corrosioncenter.ohiou.edu/freecorp>


import configparser
import constant
import anodic, cathodic
import math

config = configparser.ConfigParser()
config.read('model.ini')


# This class module contains all caculations associated with
# the electrochemical model for carbon dioxide corrosion in the
# presence of iron carbonate film and chemical reaction
# model for hydrogen sulfide corrosion

class pointmodel():
    # *******************************************************************************************
    # Definition of member variables
    # m_t: temperature in Kelvin
    # m_p: total pressure in bar
    # m_diam:pipe or cylinder diameter in m
    # m_vel: velocity of water in m/s
    # m_concHAcppm:concentration of HAc in water in ppm
    # m_concHCO3ppm: concentration of HCO3- in ppm
    # m_concFeppm: concentration of Fe++ in ppm
    # m_concH2Sppm: concentration of H2S in gas phase in ppm
    # m_concO2ppb: concentration of O2 in water in ppb
    # m_pco2: partial pressure of CO2 in bar
    # m_timemax: time required for transient H2S model,unit second
    #
    # assocated with H2S concentration profile
    # m_ConcProX():array used to store distance away from steel surface
    # m_ConcProY():array used to store concentration of H2S at different distance
    # m_isH2S:flag used to determine if H2S exists
    #
    # shared by Reaction Model
    # m_concH2CO3:concentration of H2CO3 in mol/L
    # m_concHAc: concentration of HAc in mol/L
    # m_pH: pH of water
    #
    # member variables
    # m_waterdensity: density of water in kg/m3
    # m_waterviscosity: viscosity of water in kg/m/s
    # m_diffTcorrection: modification coefficient for diffusion coefficient in terms of temperature and viscosity
    # m_kfhydration:forward rate constant of CO2 hydration reaction CO2+H2O-->H2CO3
    # m_concCO2:concentration of CO2 in water in mol/L
    # m_concHplus:concentration of H+ in mol/L
    # m_concFe:concentration of Fe++ in mol/L
    # m_concCO3:concentration of CO3= in mol/L
    # m_concHCO3:concentration of HCO3- in mol/L
    # m_ScaleFactor:scale factor related to iron carbonate film
    # m_FeSfactor: scale factor related to FeS film
    # m_SSFeCO3: supersaturation with respect to iron carbonate
    # m_K_CO2_ca: equilibrium constant for dissociation of H2CO3, H2CO3<=>H+ + HCO3-
    # m_K_CO2_bi:equilibrium constant for dissociation of HCO3-, HCO3-<=>H+ + CO3=

    # m_quant_os: amount of outer sulfide scale,mol/m2'm_Epsilon_os:porosity of outer sulfide scale,=0.9
    # m_Psi_os:tortuosity of outer sulfide scale,=0.002
    # m_delta_os:thickness of outer sulfide scale,m

    # m_CRH2S_ave:time averaged H2S corrosion rate,mm/y
    # m_CRCO2_ave:CO2 corrosion rate,mm/y
    # m_CRH2S_new:H2S corrosion rate after first time step,mm/y
    # m_Corpot:corrosion potential, V vs SHE
    # m_icorr:current density,A/m2
    # m_CR_new:H2S corrosion rate at each time step, mm/y
    #
    # collection of reaction curves
    # m_CathodicReactionArray: collection of cathodic reactions
    # m_AnodicReactionArray:collection of anodic reactions
    #
    # m_rpmTag:tag for rotating cylinder
    #
    # tags for selected species in H2S Model, used for expert system
    # m_HplusTag:tag for existence of H+ in H2S model
    # m_CO2Tag:tag for existence of CO2 in H2S model
    # m_HACTag:tag for existence of HAc in H2S model

    # details for H2S Model
    # m_CRpercentH2S_H2S:contribution of H2S in final corrosion rate
    # m_CRpercentHplus_H2S: contribution of H+ in final corrosion rate
    # m_CRpercentH2CO3_H2S: contribution of H2CO3 in final corrosion rate
    # m_CRpercentHAc_H2S:contribution of HAc in final corrosion rate
    # m_CRpercentH2O_H2S:contribution of H2O in final corrosion rate
    # ******************************************************************************************

    # In connection with user input

    #
    # This subroutine calculates some common variables shared by various subroutines in this module,
    # such as water density, water viscosity, some reaction rate constants,
    # some reaction equilibrium constants and supersaturation of iron carbonate,etc.
    def __init__(self, temperature, pressure, velocity, diameter, pco2, concH2Sppm, concFeppm, concHAcppm, pH, time=24):
        # declare a variable for temperature in celsius
        self.m_t = temperature
        self.m_pH = pH
        self.m_pco2 = pco2
        self.m_CO2Tag = True if self.m_pco2 > 0 else False
        self.m_p = pressure
        self.m_concHAcppm = concHAcppm
        self.m_HplusTag = True
        self.m_concFeppm = concFeppm
        self.m_H2STag=True if concH2Sppm > 0 else False
        self.m_concH2Sppm = concH2Sppm
        self.m_HAcTag = True if concHAcppm > 0 else False
        # velocity of the water
        self.m_vel = velocity
        # diamter of tge pipeline
        self.m_diam = diameter
        # assume reaction time takes 24 hours to complete
        self.m_timemax = time
        CentiGrade = self.m_t - 273.15
        # water density
        self.m_waterdensity = 1152.3 - 0.5116 * (self.m_t)
        # water viscosity
        self.m_waterviscosity = 1.002 * 10 ** (
                (1.3272 * (20 - CentiGrade) - 0.001053 * (CentiGrade - 20) ** 2) / (CentiGrade + 105)) / 1000
        # modification coefficient for diffusion coefficient in terms of temperature and viscosity
        self.m_diffTcorrection = (self.m_t) / 298.16 * 0.0009 / self.m_waterviscosity
        # forward reaction rate constant of CO2 hydration: CO2+H2O->H2CO3
        self.m_kfhydration = 10 ** (169.2 - 53 * math.log(self.m_t) / math.log(10) - 11715 / self.m_t)
        # calculate concentration of H+ based on pH, mol/L
        self.m_concHplus = math.exp(self.m_pH / (-0.4343))
        # concentration of CO2 in water, mol/L
        self.m_concCO2 = self.m_pco2 * self.getHenryconst_CO2()
        # concentration of H2CO3 in water,mol/L
        self.m_concH2CO3 = constant.Khyd * self.m_concCO2
        # concentration of HAc in water,mol/L
        self.m_concHAc = self.m_concHAcppm / constant.MHAc / 1000000 * 1000
        # equilibrium constant of H2CO3 dissociation: H2CO3<=>H+ + HCO3-
        self.m_K_CO2_ca = self.DisConst1st_CO2()
        # concentration of HCO3- in water, mol/L
        self.m_concHCO3 = self.m_K_CO2_ca * self.m_concH2CO3 / self.m_concHplus
        # concentration of HCO3 in water,ppm
        self.m_concHCO3ppm = self.m_concHCO3 * constant.MHCO3 * 1000000 / 1000
        # concentratio of Fe++ in water,mol/L
        self.m_concFe = self.m_concFeppm / 55.84 / 1000000 * 1000
        # equilibrium constant of HCO3- dissociation: HCO3- <=> H+ + CO3=
        self.m_K_CO2_bi = self.DisConst2nd_CO2()
        # concentration of CO3= in water, mol/L
        self.m_concCO3 = self.m_K_CO2_bi * self.m_concHCO3 / self.m_concHplus  # mol/l
        # supersaturation of iron carbonate
        self.m_SSFeCO3 = self.m_concFe * self.m_concCO3 / (10 ** (
                -59.3498 - 0.041377 * self.m_t - 2.1963 / self.m_t + 24.5724 * math.log(
            self.m_t) / 2.303 + 2.518 * constant.ionic ** 0.5 - 0.657 * constant.ionic))
        # Assuming the tag is for pipeline
        self.m_rpmTag = False
        # Assume no Oxygen soluble in the water
        self.m_concO2ppb = 0
        # Initial scalefactor induced by FeCO3 and FeS equal to 1
        self.m_ScaleFactor = 1
        self.m_FeSfactor = 1

        # calculate scalefactor induced by FeS if H2S exists
        if self.m_concH2Sppm != 0:
            self.CalFeSfactor()
        # calculate scalefactor induced by FeCO3
        self.CalScaleFactor()
        self.InitializeCurves()

    #
    # In this sub routine, two new colletions are instantiated to
    # store cathodic and anodic reactions,respectively.
    #

    def InitializeCurves(self):
        #
        #	 'collection for cathodic reactions
        self.m_CathodicReactionArray = []
        #	 'collection for anodic reactions
        self.m_AnodicReactionArray = []

    #
    # This function returns final corrosion rate.
    # In this subroutine, two kinds of corrosion rate are calculated
    # based on electrochemical and chemical reaction model, respectively
    # given H2S exists in the system. These corrosion rates are then compared
    # with each other,the higher one is reported as final corrosion rate.
    # Accordingly, the mechanism that gives higher corrosion rate is taken
    # to be the dominant corrosion mechanism.

    def getCR(self):
    #counter for total number of reactions
    #get common variables
    #Call Initialize
    #check if there are anodic reactions with inappropraite parameters
        self.m_AnodicReactionArray[:] = [reaction for reaction in self.m_AnodicReactionArray if reaction.m_status]
        # check if there are cathodic reactions with inappropraite parameters
        self.m_CathodicReactionArray[:] = [reaction for reaction in self.m_CathodicReactionArray if reaction.m_status]
        #calculate corrosion potential in electrochemical model
        self.m_Corpot=self.CalCorpot()
        #calculate corrosion rate generated by electrochemical model
        self.CalCO2CorrosionRate()
        #if H2S exists, calcualte corrosion rate generated by H2S model
        if self.m_concH2Sppm > 0:
           self.transient_model()

        #compare H2S corrosion rate and CO2 corrosion rate if H2S exists and
        #take higher one as the final corrosion rate
        if self.m_concH2Sppm > 0:
        #m_CRH2S_new is the H2S corrosion rate generated in the first time step,
        #this is used for comparison to ensure a higher H2S corrosion rate at
        #high H2S concentration,which is independent of the time span.
            if self.m_CRCO2_ave < self.m_CRH2S_new:
                    return self.m_CR_new
                    self.m_isH2S = True
            else:
                return self.m_CRCO2_ave
                self.m_isH2S = False

        else:
            #if H2S does not exist, CO2 corrosion rate would be the final corrosion rate
            return self.m_CRCO2_ave
            self.m_isH2S = False

    #
    # This subroutine calculates CO2 corrosion rate based on electrochemical theoery.
    #

    def CalCO2CorrosionRate(self):

        # counter for number of reactions
        # Reaindex As Integer
        # corrrosion potential
        # Ecorr As Double
        # corrosion current density
        # icorr As Double

        # Corrosion potential
        Ecorr = self.m_Corpot
        # corrosion current density equals the sum of current density of all cathodic reactions
        icorr = 0
        for Reaindex in self.m_CathodicReactionArray:
            icorr = icorr + Reaindex.get_currentdensity(Ecorr, self.m_ScaleFactor, self.m_FeSfactor)

        # Calculate the contributions of corrosive species involved in cathodic reactfgetzions to overall corrosion rate
        for Reaindex in self.m_CathodicReactionArray:
            if Reaindex.m_status:
                Reaindex.m_contribution = Reaindex.get_currentdensity(Ecorr, self.m_ScaleFactor,
                                                                      self.m_FeSfactor) / icorr
            else:
                Reaindex.m_contribution = 0
        self.m_icorr = icorr
        # conversion from current density in A/m2 to corrsoion rate of iron in mm/y
        self.m_CRCO2_ave = 1.155 * icorr

    #
    # This subroutine calculates corrosion potential based on
    # the fact that total anodic current density equals total cathodic current density at equilibrium.
    # A bi-section algorithm is applied for solving the equation: itotal(anodic)=itotal(cathodic).
    #

    def CalCorpot(self):
        # two variables used for storing values of function for upper and lower limit.
        # fx1 As Double, fx2 As Double
        # used to store the potential in the middle of upper and lower potentials.
        # root As Double
        # set up upper and lower limit of potential in which real potential resides
        # llimit As Double, ulimit As Double
        # tolerence used in this numerical methode
        epsilon = 0.0000001

        llimit = -1  # lower limit of potential is -1V
        ulimit = 0  # upper limit of potential is 0V

        # calculate scale factor related to FeCO3
        self.m_ScaleFactor = 1
        self.m_FeSfactor = 1
        self.CalScaleFactor()

        # solve corrosion potential using bi-section algorithm
        while (ulimit - llimit) > epsilon:
            root = ulimit
            fx1 = self.getNetCurrentDensity(root)
            root = (ulimit + llimit) / 2
            fx2 = self.getNetCurrentDensity(root)
            if (fx1 * fx2) > 0:
                ulimit = root
                llimit = llimit
            else:
                llimit = root

        return root

    #
    # This function returns a net current density for a given
    # corrosion potential.
    #
    def getNetCurrentDensity(self, Ecorr):

        # counter for number of reactions
        # Reaindex As Integer
        # Func As Double
        tempo = 0
        for Reaindex in self.m_CathodicReactionArray:
            tempo = tempo - Reaindex.get_currentdensity(Ecorr, self.m_ScaleFactor, self.m_FeSfactor)
        for Reaindex in self.m_AnodicReactionArray:
            tempo = tempo + Reaindex.get_currentdensity(Ecorr, self.m_ScaleFactor, self.m_FeSfactor)
        return tempo

    #
    # $This function returns a pH value based on HCO3- concentration
    #

    def CalpH(self):
        # concentration of CO2 in water,mol/L
        self.m_concCO2 = self.m_pco2 * self.getHenryconst_CO2()
        # concentration of H2CO3 in water,mol/L
        self.m_concH2CO3 = constant.Khyd * self.m_concCO2
        # equilibrium constant for H2CO3 dissociation: H2CO3 <=> H+ + HCO3-
        self.m_K_CO2_ca = self.DisConst1st_CO2()
        # concentration of HCO3- in water,mol/L
        self.m_concHCO3 = self.m_concHCO3ppm / constant.MHCO3 / 1000000 * 1000
        # concentration of H+,mol/L
        self.m_concHplus = self.m_K_CO2_ca * self.m_concH2CO3 / self.m_concHCO3

        self.m_pH = -0.4343 * math.log(self.m_concHplus)

    #
    # This subroutine calculates Scale factor of FeCO3.
    #

    def CalScaleFactor(self):

        # solubility limit of FeCO3
        # KspFeCO3 As Double
        # precipitation rate of FeCO3
        # PRFeCO3 As Double
        # temperature in celsius
        # CentiGrade As Double
        # CentiGrade = m_t - 273.15

        # Calculates scale formation

        # calculate solubility limit for iron carbonate
        KspFeCO3 = 10 ** (-59.3498 - 0.041377 * self.m_t - 2.1963 / self.m_t + 24.5724 * math.log(
            self.m_t) / 2.303 + 2.518 * constant.ionic ** 0.5 - 0.657 * constant.ionic)

        # calculate supersaturation of iron carbonate
        self.m_SSFeCO3 = self.m_concFe * self.m_concCO3 / KspFeCO3

        # precipitation rate of iron carbonate,mol/m^2/s
        if self.m_SSFeCO3 > 1:
            PRFeCO3 = math.exp(21.3 - 64851.4 / constant.r / self.m_t) * KspFeCO3 * (self.m_SSFeCO3 - 1)

        else:
            PRFeCO3 = 0  # no precipitation when supersaturation less than 1

        # empirically determined function for scale factor based on supersaturation of FeCO3
        if self.m_SSFeCO3 > 1:
            self.m_ScaleFactor = self.m_SSFeCO3 ** (-0.7)
        else:
            self.m_ScaleFactor = 1

    #
    # This subroutine calculates FeS Factor using an empirical equation.
    #

    def CalFeSfactor(self):
        ConstFeS = 0.05
        self.m_FeSfactor = 1 / (1 + ConstFeS * self.m_concH2Sppm)

    #
    # This function returns equilibrium constant for dissociation of H2CO3
    # H2CO3<=>H+ + HCO3-
    #

    def DisConst1st_CO2(self):
        # temperature in Fahrenheit
        # Fahrenheit As Double
        Fahrenheit = (self.m_t - 273.15) * 9 / 5 + 32

        tempo = 387.6 * 10 ** (-(6.41 - 1.594 * 10 ** (-3) * (Fahrenheit) + 8.52 * 10 ** (-6) \
                                 * (Fahrenheit) ** 2 - 3.07 * 10 ** (-5) * 14.5 * self.m_pco2 - 0.4772 \
                                 * constant.ionic ** 0.5 + 0.118 * constant.ionic))
        return tempo

    #
    # This function returns equilibrium constant for dissociation of HCO3-
    # * HCO3- <=>H+ +CO32-
    #

    def DisConst2nd_CO2(self):
        # temperature in Fahrenheit
        # Fahrenheit As Double
        Fahrenheit = (self.m_t - 273.15) * 9 / 5 + 32

        tempo = 10 ** (-(10.61 - 4.97 * 0.001 * (Fahrenheit) + 1.331 * 10 ** (-5) * (Fahrenheit) ** 2 \
                         - 2.624 * 0.00001 * self.m_pco2 * 14.5 - 1.166 * constant.ionic ** (
                             0.5) + 0.3466 * constant.ionic))

        return tempo

    #
    # This function returns Henry constant for CO2
    # CO2(g)<=>CO2(w)

    def getHenryconst_CO2(self):
        # temperature in celsuis
        # CentiGrade As Double
        # Fahrenheit As Double

        Fahrenheit = (self.m_t - 273.15) * 9 / 5 + 32
        CentiGrade = self.m_t - 273.15

        tempo = 14.5 / 1.00258 * 10 ** (
            -(2.27 + 0.00565 * Fahrenheit - 0.00000806 * Fahrenheit ** 2 + 0.075 * constant.ionic))
        return tempo

    #
    # This function returns mass transfer coefficient
    #

    def getKm(self, Diff):
        if self.m_rpmTag == True:
            # mass transfer coefficient for rotating cylinder
            tempo = 0.0791 * (self.m_vel * self.m_diam * self.m_waterdensity / self.m_waterviscosity) ** 0.7 \
                    * (self.m_waterviscosity / self.m_waterdensity / Diff) ** 0.356 * Diff / self.m_diam

        else:
            # mass transfer coefficient for pipelines
            tempo = 0.0165 * (self.m_vel * self.m_diam * self.m_waterdensity / self.m_waterviscosity) ** 0.86 \
                    * (self.m_waterviscosity / self.m_waterdensity / Diff) ** 0.33 * Diff / self.m_diam
        return tempo

    #
    # This function returns fluxes of species(except CO2) through mass transfer layers in H2S system
    #

    def getFlux(self, flux1, flux2, conc_b, Diff, km, conc_s, expterm):
        root = flux1
        fx1 = conc_b - root * (
                self.m_delta_os / (Diff * self.m_Epsilon_os * self.m_Psi_os) + 1 / km) - conc_s * math.exp(
            root / expterm)
        while abs(fx1) / conc_s > 0.0001:
            root = flux1
            fx1 = conc_b - root * (
                    self.m_delta_os / (Diff * self.m_Epsilon_os * self.m_Psi_os) + 1 / km) - conc_s * math.exp(
                root / expterm)
            root = (flux1 + flux2) / 2
            fx2 = conc_b \
                  - root * (self.m_delta_os / (Diff * self.m_Epsilon_os * self.m_Psi_os) + 1 / km) \
                  - conc_s * math.exp(root / expterm)
            if (fx1 * fx2) > 0:
                flux1 = root
                flux2 = flux2
            else:
                flux2 = root
        return root

    #
    # This function returns fluxes of CO2 through mass transfer layers in H2S system
    # equation 102 in Freecorp 1.0

    def getFluxCO2(self, flux1, flux2, conc_b, Diff, km, conc_s, para, expterm):

        root = flux1
        conc_s = root / para
        fx1 = conc_b - root * (self.m_delta_os / (Diff * self.m_Epsilon_os * self.m_Psi_os) + 1 / km) \
              - conc_s * math.exp(root / expterm)
        while abs(fx1) / conc_s > 0.0001:
            root = flux1
            conc_s = root / para
            fx1 = conc_b - root * (self.m_delta_os / (Diff * self.m_Epsilon_os * self.m_Psi_os) + 1 / km) \
                  - conc_s * math.exp(root / expterm)

            root = (flux1 + flux2) / 2
            conc_s = root / para

            fx2 = conc_b \
                  - root * (self.m_delta_os / (Diff * self.m_Epsilon_os * self.m_Psi_os) + 1 / km) \
                  - conc_s * math.exp(root / expterm)

            if (fx1 * fx2) > 0:
                flux1 = root
                flux2 = flux2
            else:
                flux2 = root

        return root

    #
    # This function checks if species exist.
    # Returns false if the species don't exist.
    # Return ture if the species exist.
    #

    def Isexist(self, speciesconc):

        try:
            if (speciesconc > 0):
                return True
            else:
                return False
        except:
            return False

    #
    # This function returns a mass transfer limiting current density for H+ reduction
    #
    def getlimitcur_Hplus(self):

        # the constant diffusion coefficient of H+ under reference conditions
        # DiffrefHplus As Double
        # real diffusion coefficient in the system
        # DiffHplus As Double
        # mass transfer coeffcient of H+
        # kmHplus As Double

        DiffrefHplus = 0.00000000931
        DiffHplus = DiffrefHplus * self.m_diffTcorrection
        kmHplus = self.getKm(DiffHplus)
        # calculate mass transfer limiting current density of H+ reduction,A/m2
        return kmHplus * self.m_ScaleFactor * self.m_FeSfactor * constant.F * 10 ** (-self.m_pH) * 1000

    #
    # This function returns the chemical reaction limiting current density for H2CO3 reduction
    # It is limited by the slow reaction CO2+H2O->H2CO3

    def getlimitcur_H2CO3(self):

        # the constant diffusion coefficient of H2CO3 under reference conditions
        # DiffrefH2CO3 As Double
        # real diffusion coefficient of H2CO3 in the system
        # DiffH2CO3 As Double
        # mass transfer coeffcient of H2CO3
        # kmH2CO3 As Double
        # mass transfer layer thickness
        # Deltamasstransfer As Double
        # reaction layer thickness
        # Deltareaction As Double
        # flow factor related to relative effect of mass transfer and chemical reaction
        # flowfactor As Double
        # limiting current density of H2CO3 reduction
        # ilimreactionH2CO3 As Double

        DiffrefH2CO3 = 0.0000000013
        DiffH2CO3 = DiffrefH2CO3 * self.m_diffTcorrection
        ilimreactionH2CO3 = constant.F * self.m_concCO2 * 1000 * (
                DiffH2CO3 * self.m_ScaleFactor * self.m_FeSfactor * constant.Khyd * self.m_kfhydration) ** 0.5

        kmH2CO3 = self.getKm(DiffH2CO3)

        # when FeCO3 film exists, flow factor becomes insignificant
        if self.m_ScaleFactor != 1:
            Deltamasstransfer = DiffH2CO3 / kmH2CO3
            Deltareaction = (DiffH2CO3 * constant.Khyd / self.m_kfhydration) ** 2
            try:
                flowfactor = (1 + math.exp(-2 * Deltamasstransfer / Deltareaction))/(1 - math.exp(-2 * Deltamasstransfer / Deltareaction))
            except:
                flowfactor = 1
        else:
            flowfactor = 1

        return ilimreactionH2CO3 * flowfactor  # A/m2

    #
    # This function returns a mass transfer limiting current density for HAc reduction
    #

    def getlimitcur_HAc(self):

        # the constant diffusion coefficient of HAc under reference conditions
        # DiffrefHAc As Double
        # real diffusion coefficient of HAc in the system
        # DiffHAC As Double
        # mass transfer coeffcient of HAc
        # kmHAc As Double

        DiffrefHAc = 0.0000000005
        DiffHAC = DiffrefHAc * self.m_diffTcorrection
        kmHAc = self.getKm(DiffHAC)
        tempo = kmHAc * self.m_ScaleFactor * self.m_FeSfactor * constant.F * self.m_concHAc * 1000  # A/m2
        return tempo

    #
    # This function returns the mass tranfer limiting current density for O2 reduction
    #

    def getlimitcur_O2(self):

        # the constant diffusion coefficient of O2 under reference conditions
        # DiffrefO2 As Double
        # real diffusion coefficient of O2 in the system
        # DiffO2 As Double
        # mass transfer coeffcient of O2
        # kmO2 As Double
        DiffrefO2 = 0.000000002
        DiffO2 = DiffrefO2 * self.m_diffTcorrection
        kmO2 = self.getKm(DiffO2)
        return 4 * kmO2 * self.m_ScaleFactor * self.m_FeSfactor * constant.F * self.m_concO2ppb / 1000 / 32  # A/m2

    #
    # This function returns a HAc concentration value based on given
    # Ac- concentration and current pH value.
    #

    def getHAc(self, concAcppm):
        # concentration of Ac-,H+, HAc,respectively
        # concAc, concHplus, concHAc As Double
        # equilibrium constant for HAc dissociation: HAc<=> H+ + Ac-
        # K_HAc_ha As Double

        concHplus = 10 ** (-self.m_pH)  # mol/L
        concAc = concAcppm / 59 / 1000000 * 1000  # mol/L
        K_HAc_ha = 10 ** (-6.66104 + 0.0134916 * (self.m_t) - 2.37856 * 10 ** (-5) * (self.m_t) ** 2)
        concHAc = concAc * concHplus / K_HAc_ha  # mol/L
        return concHAc * MHAc * 1000000 / 1000  # ppm

    #
    # This function returns a Ac- concentration value by given
    # HAc concentration and current pH value.
    #

    def getAc(self, concHAcppm):

        # concentration of Ac-,H+, HAc,respectively
        # concAc, concHplus, concHAc As Double
        # equilibrium constant for HAc dissociation: HAc<=> H+ + Ac-
        # K_HAc_ha As Double

        concHplus = 10 ** (-self.m_pH)  # mol/L
        K_HAc_ha = 10 ** (-6.66104 + 0.0134916 * (self.m_t) - 2.37856 * 10 ** (-5) * (self.m_t) ** 2)
        concAc = concHAcppm * 1000 * K_HAc_ha / concHplus / MHAc / 1000000  # mol/L
        return concAc * 59 * 1000000 / 1000  # ppm

    #
    # This function returns a pH value based on given
    # HCO3 concentration.
    # HCO3-<=>H+ + CO3=
    #

    def getPH(self, concHCO3ppm):
        # concentration of HCO3-,CO2,H2CO3,H+,respectively
        # concHCO3, concCO2, concH2CO3, concHplus As Double
        # equilibrium constant for dissociation of HCO3-
        # K_CO2_ca As Double
        # temperature in celsius and fahrenheit, respectively
        # CentiGrade, Fahrenheit As Double

        CentiGrade = self.m_t - 273.15
        Fahrenheit = (self.m_t - 273.15) * 9 / 5 + 32

        concHCO3 = concHCO3ppm / constant.MHCO3 / 1000000 * 1000  # mol/L

        concCO2 = self.m_pco2 * 14.5 / 1.00258 * 10 ** (-(2.27 + 0.00565 * constant.Fahrenheit \
                                                          - 0.00000806 * constant.Fahrenheit ** 2 + 0.075 * ionic))  # mol/L

        concH2CO3 = constant.Khyd * concCO2  # mol/L
        K_CO2_ca = 387.6 * 10 ** (-(6.41 - 1.594 * 10 ** (-3) * (constant.Fahrenheit) \
                                    + 8.52 * 10 ** (-6) * (Fahrenheit) ** 2 - 3.07 * 10 ** (-5) * 14.5 \
                                    - 0.4772 * ionic ** 0.5 + 0.118 * ionic))
        concHplus = K_CO2_ca * concH2CO3 / concHCO3  # mol/L
        return -0.4343 * math.log(concHplus)

    #
    # This subroutine add one cathodic reaction into the
    # cathodic reaction colletion.
    #

    def AddCathodicReactionCurve(self, newCurve):
        #  i As Integer
        limitcur_Hplus = self.getlimitcur_Hplus()
        limitcur_H2CO3 = self.getlimitcur_H2CO3()
        limitcur_HAc = self.getlimitcur_HAc()
        limitcur_O2 = self.getlimitcur_O2()
        newCathodicReaction = cathodic.cathodic(newCurve, self.m_t, self.m_pH, limitcur_Hplus, self.m_concH2CO3,
                                                limitcur_H2CO3 \
                                                , self.m_concHAc, limitcur_HAc, limitcur_O2)
        self.m_CathodicReactionArray.append(newCathodicReaction)

    #
    # This subroutine add one anodic reaction to the
    # anodic reaction  colletion.
    #

    def AddAnodicReactionCurve(self, newCurve):

        newAnodicReaction = anodic.anodic(newCurve, self.m_t)
        self.m_AnodicReactionArray.append(newAnodicReaction)

    #
    # This subroutine add common anodic and cathodic reactions curves
    # and delete the existing curves
    def AddCommonReactionCurve(self):
        self.InitializeCurves()
        self.AddAnodicReactionCurve("Iron_dissolution")
        self.AddCathodicReactionCurve("H2CO3_reduction")
        self.AddCathodicReactionCurve("H2O_reduction")
        self.AddCathodicReactionCurve("HAc_reduction")
        #		self.AddCathodicReactionCurve("O2_reduction")
        self.AddCathodicReactionCurve("H_+reduction")

    #
    # This subroutine calculates the H2S corrosion rate in a transient model
    # Three layers are assumed to exist on steel surface: inner Mackinawite film,
    # outer Mackinawite scale and liquid mass transfer boundary layer.
    # A fast chemical reaction between H2S and iron is the cause of corrosion which
    # is limited by mass transfer of species through three layers on steel surface.

    def transient_model(self):

        CentiGrade = self.m_t - 273.15
        # -----------------------------------------------------------------------------------------------
        # VARIABLE INTIALIZATION
        time = 0  # s
        # duration of simulation
        # time step
        deltat = 1  # s
        # time between writing to sheet
        Timetowrite = 60  # s
        Timetowrite_counter = 100
        # amount of outer sulfide scale
        self.m_quant_os = 0  # mol / m2
        # outer scale porosity
        self.m_Epsilon_os = 0.9
        # outer scale tortuosity
        self.m_Psi_os = 0.002
        # outer scale initial thickness
        self.m_delta_os = 0  # m

        # -----------------------------------------------------------------------------------------------
        # PHYSICAL PROPERTIES
        # H2S gas concentration
        pH2S = self.m_concH2Sppm / 1000000 * self.m_p * constant.MCO2 / constant.MH2S  # ppm

        # -----------------------------------------------------------------------------------------------
        #  H+ related constants
        # kinetic constant for H+ solid state diffusion
        A_Hplus = 0.0000004  # kmol/m2/s
        # activation energy for H+ solid state diffusion
        E_Hplus = 15500  # J/mol
        # surface concetration of H+ (reference value)
        conc_s_Hplus = 0.0000001  # mol/l
        # standard diffusion coefficient for H+ in water
        DiffHplus = 0.00000000931  # m^2/s
        DiffHplus = DiffHplus * self.m_diffTcorrection  # the mass transfer coefficient for H+
        kmHplus = self.getKm(DiffHplus)
        # Calculation of the H+ concetration in the bulk
        conc_b_Hplus = self.m_concHplus
        # Maximum flux of H+ (without any film or scale)
        flux_Hplus_max = kmHplus * conc_b_Hplus
        # exponential term in the Hplus flux equation
        #    expterm_Hplus = A_Hplus * Exp(-E_Hplus / (R * 1000) / (m_t))
        expterm_Hplus = A_Hplus

        # -----------------------------------------------------------------------------------------------
        # H2S related constants
        # kinetic constant for H2S solid state diffusion
        A_H2S = 0.00000002  # kmol/m2/s
        # activation energy for H2S solid state diffusion
        E_H2S = 15500  # J/mol
        # surface concetration of dissolved H2S (reference value)
        conc_s_H2S = 0.0000001  # mol/l
        # Solubility of H2S
        K_H2S_sol = 10 ** (-(
                634.27 + 0.2709 * self.m_t - 0.00011132 * self.m_t ** 2 - 16719 / self.m_t - 261.9 * math.log(
            self.m_t) / math.log(10)))

        # dissoved H2S (first) dissociation constant
        K_H2S = 10 ** (-(15.345 - 0.045676 * (self.m_t) + 5.9666 * 0.00001 * (self.m_t) * (self.m_t)))
        # bisulfide ion HS- (second) dissociation contant
        # sulfide ion S= concentration

        # diffusion coefficient for H2S in water
        Diff_H2S = 0.0000000013  # m^2/s
        # Calculation of the mass transfer coefficient for H2S
        Diff_H2S = Diff_H2S * self.m_diffTcorrection
        km_H2S = self.getKm(Diff_H2S)
        # Calculation of the dissolved H2S concetration in the bulk
        # bulk concetration of dissolved H2S
        conc_b_H2S = K_H2S_sol * pH2S  # mol/l

        # solubility product constant for iron sulfide for reaction: H2S + Fe2+ -> FeS + 2H+
        KspFeS = 10 ** (2848.779 / self.m_t - 6.347)
        # supersaturation of iron sulfide
        SSFeS = self.m_concFe * conc_b_H2S / conc_b_Hplus ** 2 / KspFeS

        # Maximum flux of H2S (without any film or scale)
        flux_H2S_max = km_H2S * conc_b_H2S
        # exponential term in the H2S flux equation
        # expterm_H2S = A_H2S * Exp(-E_H2S / (R * 1000) / (m_t))
        expterm_H2S = A_H2S

        # -----------------------------------------------------------------------------------------------
        #  CO2 related constants
        if self.Isexist(self.m_pco2):
            # kinetic constant for CO2 solid state diffusion
            A_CO2 = 0.000000002  # kmol/m2/s
            # activation energy for CO2 solid state diffusion
            E_CO2 = 15500  # J/mol
            # standard diffusion coefficient for CO2 in water
            DiffCO2 = 0.00000000196  # m^2/s
            # Calculation of the mass transfer coefficient for CO2
            DiffCO2 = DiffCO2 * self.m_diffTcorrection
            kmCO2 = self.getKm(DiffCO2)
            # standard diffusion coefficient for H2CO3 in water
            DiffrefH2CO3 = 0.0000000013  # m^2/s
            DiffH2CO3 = DiffrefH2CO3 * self.m_diffTcorrection
            # Calculation of the dissolved CO2 concetration in the bulk
            conc_b_CO2 = self.m_pco2 * self.getHenryconst_CO2()
            # Maximum flux of CO2 (without any film or scale)
            flux_CO2_max = conc_b_CO2 * (DiffH2CO3 * constant.Khyd * self.m_kfhydration) ** 0.5
            # exponential term in the CO2 flux equation
            # expterm_CO2 = A_CO2 * Exp(-E_CO2 / (R * 1000) / (m_t))
            expterm_CO2 = A_CO2  # HOLD

        # -----------------------------------------------------------------------------------------------
        #  HAc related constants
        if self.Isexist(self.m_concHAcppm):
            # kinetic constant for HAc solid state diffusion
            A_HAc = 0.000000002  # kmol/m2/s
            # activation energy for HAc solid state diffusion
            E_HAc = 15500  # J/mol
            # surface concetration of HAc (reference value)
            conc_s_HAc = 0.0000001  # mol/l
            # standard diffusion coefficient for HAc in water
            DiffHAC = 0.0000000005  # m^2/s
            # Calculation of the mass transfer coefficient for HAc
            DiffHAC = DiffHAC * self.m_diffTcorrection
            kmHAc = self.getKm(DiffHAC)
            # Calculation of the  HAc concetration in the bulk
            conc_b_HAc = self.m_concHAc  # mol/l
            # Maximum flux of HAc (without any film or scale)
            flux_HAc_max = kmHAc * conc_b_HAc
            # exponential term in the CO2 flux equation
            # expterm_HAc = A_HAc * math.exp(-E_HAc / (constant.R * 1000) / (m_t))
            expterm_HAc = A_HAc
        CR_t, CR_H2S, CR_Hplus, CR_CO2, CR_HAc, delta_cum = 0, 0, 0, 0, 0, 0
        # -----------------------------------------------------------------------------------------------
        # Calculate corrosion rate of H2S
        while time < self.m_timemax * 3600:

            # calculate corrosion rate due to H2S
            if self.m_H2STag:
                flux1 = 0.00001
                flux2 = 0
                flux_H2S = flux1
                flux_H2S = self.getFlux(flux1, flux2, conc_b_H2S, Diff_H2S, km_H2S, conc_s_H2S, expterm_H2S)  # mol/(m2.s)

                CR_H2S = flux_H2S * 1000 * constant.MFe / constant.densityFe * 3600 * 24 * 365  # mm/y
                CR_H2S_max = flux_H2S_max * 1000 * constant.MFe / constant.densityFe * 3600 * 24 * 365  # mm/y
                # outer SULFIDE scale thickness growth by solid state reaction and spaling
                self.m_delta_os = self.m_delta_os + flux_H2S * deltat * constant.MFeS / constant.densityFeS / (
                        1 - self.m_Epsilon_os) * 0.5  # m
                self.m_quant_os = self.m_quant_os + flux_H2S * deltat * 1000  # mol/m2

            # calcualte corrosion rate due to H+
            if self.m_HplusTag:
                flux1 = 0.0001
                flux2 = 0
                flux_Hplus = flux1
                flux_Hplus = self.getFlux(flux1, flux2, conc_b_Hplus, DiffHplus, kmHplus, conc_s_Hplus,
                                          expterm_Hplus)  # mol/(m2.s)
                CR_Hplus = flux_Hplus * 1000 * constant.MFe / constant.densityFe * 3600 * 24 * 365  # mm/y
                CR_Hplus_max = flux_Hplus_max * 1000 * constant.MFe / constant.densityFe * 3600 * 24 * 365  # mm/y
            # calcualte H2S corrosion rate due to CO2
            if self.m_CO2Tag:
                if self.Isexist(self.m_pco2):
                    flux1 = 0.000001
                    flux2 = 0
                    flux_CO2 = flux1
                    Parameter_CO2 = (
                                            DiffH2CO3 * self.m_Epsilon_os * self.m_Psi_os * constant.Khyd * self.m_kfhydration) ** 0.5
                    # surface concentration of CO2,mol/m3
                    conc_s_CO2 = flux_CO2 / Parameter_CO2
                    flux_CO2 = self.getFluxCO2(flux1, flux2, conc_b_CO2, DiffCO2, kmCO2, conc_s_CO2, Parameter_CO2,
                                               expterm_CO2)  # mol/(m2.s)
                    CR_CO2 = flux_CO2 * 1000 * constant.MFe / constant.densityFe * 3600 * 24 * 365  # mm/y
                    CR_CO2_max = flux_CO2_max * 1000 * constant.MFe / constant.densityFe * 3600 * 24 * 365  # mm/y
            # calcualte corrosion rate due to HAC
            if self.m_HAcTag:
                if self.Isexist(self.m_concHAcppm):
                    flux1 = 0.000001
                    flux2 = 0
                    flux_HAc = flux1
                    flux_HAc = self.getFlux(flux1, flux2, conc_b_HAc, DiffHAC, kmHAc, conc_s_HAc,
                                            expterm_HAc)  # mol/(m2.s)

                CR_HAc = flux_HAc * 1000 * constant.MFe / constant.densityFe * 3600 * 24 * 365  # mm/y
                CR_HAc_max = flux_HAc_max * 1000 * constant.MFe / constant.densityFe * 3600 * 24 * 365  # mm/y

            # time step,s
            deltat = 1.1 * deltat

            # total time,s
            time = time + deltat
            #    Timetowrite_counter = Timetowrite_counter + deltat

            # -----------------------------------------------------------------------------------------------
            # total corrosion rate
            CR_t = CR_H2S + CR_Hplus + CR_CO2 + CR_HAc  # mm/y
            self.m_CR_new = CR_t

            if time < 2:
                self.m_CRH2S_new = self.m_CR_new  # corrosion rate after the first time step
            # cummulative corrosion damage
            delta_cum = delta_cum + (CR_H2S + CR_Hplus + CR_CO2 + CR_HAc) * deltat / 3600 / 24 / 365  # mm
        # ****************************************************
        # Calculate contribution of species to overall corrosion
        # rate based on respective fluxes.
        self.m_CRpercentH2S_H2S = CR_H2S / CR_t * 100
        self.m_CRpercentHplus_H2S = CR_Hplus / CR_t * 100
        self.m_CRpercentH2CO3_H2S = CR_CO2 / CR_t * 100
        self.m_CRpercentHAc_H2S = CR_HAc / CR_t * 100
        self.m_CRpercentH2O_H2S = 0
        self.m_CRpercentO2_H2S = 0
        # ****************************************************
        # time averaged corrosion rate
        self.m_CRH2S_ave = delta_cum / time * 3600 * 24 * 365

        # generate data for plotting concentration profile of H2S
        # at steel surface
        self.m_ConcProX = [0, 0, 0, 0, 0]
        self.m_ConcProY = [0, 0, 0, 0, 0]
        self.m_ConcProX[0] = 0.001
        self.m_ConcProY[0] = 0.0000001 * 1000

        # at inner Makinawite film
        self.m_ConcProX[1] = 0.01
        self.m_ConcProY[1] = math.exp(flux_H2S / expterm_H2S) * 0.0000001 * 1000

        # at outer Makinawite scale
        self.m_ConcProX[2] = 0.01 + self.m_delta_os * 1000000
        self.m_ConcProY[2] = (conc_b_H2S - flux_H2S / km_H2S) * 1000

        # at liquid boundary layer
        self.m_ConcProX[3] = 0.01 + self.m_delta_os * 1000000 + Diff_H2S / km_H2S * 1000000
        self.m_ConcProY[3] = conc_b_H2S * 1000

        # bulk solution
        if (0.01 + self.m_delta_os * 1000000 + Diff_H2S / km_H2S * 1000000) < 1000:
            self.m_ConcProX[4] = 1000
            self.m_ConcProY[4] = conc_b_H2S * 1000
        else:
            self.m_ConcProX[4] = 0.01 + self.m_delta_os * 1000000 + Diff_H2S / km_H2S * 1000000 + 1000
            self.m_ConcProY[4] = conc_b_H2S * 1000

    #
    # This function returns data points used in the cencentration profile.
    # return an array of data points
    #
    def getConcentrationProfile(XorY):
        if XorY == False:
            return self.m_ConcProY
        else:
            return self.m_ConcProX
