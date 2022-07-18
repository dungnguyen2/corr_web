#===========================================================================
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


#===========================================================================
#    Revision History
#
#    For revision history see code header in PointModelDialog Form
#===========================================================================
import configparser
import constant
import math
config = configparser.ConfigParser()
config.read('model.ini')

#
# This class module deals with cathodic reactions.
#
class cathodic:

#******************************************************************************
#description of member variables
#m_curvename:name of the reaction
#m_curvetype: type of the reaction, anodic or cathodic
#m_Erev: reversible potential,V
#m_b:Tafel slope,V/decade
#m_io:exchange current density,A/m2
#m_ilim:limiting current density,A/m2
#m_Epas:passivation potential,V
#m_ipas:passive current density,A/m2
#m_Epit:pitting potential,V
#m_PointModel:variable for PointModel
#m_contribution:contribution of species to overall corrosion rate
#Constant:
#******************************************************************************


#
#This subroutine initializes some common parameters used in
#further calculations, such as name and type of the reaction,
#reversible potential,Tafel slope,exchange current density
#and limiting current density,etc.
#Reference exchange current density and activation energy
#are read from the file,other parameters are calculated.
#

	def __init__(self, curvename, temperature, pH, getlimitcur_Hplus, concH2CO3, getlimitcur_H2CO3, concHAc, getlimitcur_HAc, getlimitcur_O2):

	#reversible potential,exchange current density,Tafel slope,limiting current density,respectively
	#ErevInit, ioInit, bInit, ilimdiffusionInit
	#reversible potential and Tafel slope for H+ reduction,respectively
	#ErevHplus, bHplus

	#reference exchange current density and activation energy,respectively
	#ioref As Double, ActivEnergy
	#reference concentration
	#concref	
	
		self.m_curvename=curvename
		self.m_status=True
		try:
			ErevHplus = -2.303 * constant.r * (temperature) / constant.F * pH  # V
			bHplus = 2.303 * constant.r * (temperature) / constant.F / 0.5  # V/decade
			# parameters for H+ reduction
			if	curvename== "H_+reduction":
				ioref = float(config["H_+reduction"]["ioref"])
				ActivEnergy = float(config["H_+reduction"]["ActivEnergy"])
				ioInit = ioref* 10 ** (-0.5 * (pH - 4))* math.exp(-ActivEnergy / constant.r * (1 / (temperature) - 1 / 298.16)) #A/m^2
				ErevInit = ErevHplus
				bInit = bHplus
				ilimdiffusionInit = getlimitcur_Hplus

			# parameters for H2CO3 reduction
			elif curvename == "H2CO3_reduction":
				# kinetics of H2CO3 reduction
				ioref = float(config["H2CO3_reduction"]["ioref"])
				ActivEnergy = float(config["H2CO3_reduction"]["ActivEnergy"])
				concref = float(config["H2CO3_reduction"]["concref"])
				ioInit = ioref* concH2CO3 / concref* 10 ** (0.5 * (pH - 5))* math.exp(-ActivEnergy / constant.r * (1 / (temperature) - 1 / 298.16))  #A/m^2
				ErevInit = ErevHplus
				bInit = bHplus
				ilimdiffusionInit = getlimitcur_H2CO3 #A/m2

			# parameters for HAc reduction
			elif curvename== "HAc_reduction":
				ioref = float(config["HAc_reduction"]["ioref"])
				ActivEnergy = float(config["HAc_reduction"]["ActivEnergy"])
				concref = float(config["HAc_reduction"]["concref"])
				ioInit = ioref* (concHAc / concref) ** 0.5* math.exp(-ActivEnergy / constant.r * (1 / (temperature) - 1 / 293.16)) #A/m^2
				ErevInit = ErevHplus
				bInit = bHplus
				ilimdiffusionInit = getlimitcur_HAc #A/m2

			#parameters for H2O reduction
			elif curvename== "H2O_reduction":
				ioref = float(config["H2O_reduction"]["ioref"])
				ActivEnergy = float(config["H2O_reduction"]["ActivEnergy"])
				ioInit = ioref* math.exp(-ActivEnergy / constant.r * (1 / (temperature) - 1 / 298.16))  #A/m^2
				ErevInit = ErevHplus
				bInit = bHplus
				ilimdiffusionInit = 1E+20  #A/m2

			#parameters for O2 reduction
			elif curvename== "O2_reduction":
				ioref = float(config["O2_reduction"]["ioref"])
				ActivEnergy = float(config["O2_reduction"]["ActivEnergy"])
				ioInit = ioref* math.exp(-ActivEnergy / constant.r * (1 / (temperature) - 1 / 298.16))  #A/m^2
				ErevInit = 0.5 # V
				bInit = 2.303 * constant.r * (temperature) / constant.F / 1 #V/decade
				ilimdiffusionInit = getlimitcur_O2 #A/m2

			else:
			   #for user-defined reactions,read all the parameters from the file
				if float(config["User_Cathodic_Reaction"]["Temperature(K)"]) != temperature:
					print("Reaction ",curvename," is not compatible!")
					return
				ioInit = float(config["User_Cathodic_Reaction"]["ExchangeCurrentDensity(A/m2)"])
				ErevInit = float(config["User_Cathodic_Reaction"]["ReversablePotential(V)"])
				bInit = float(config["User_Cathodic_Reaction"]["TafelSlope(b)"])
				ilimdiffusionInit = float(config["User_Cathodic_Reaction"]["LimitingCurrentDensity(A/m2)"])

			self.m_Erev = ErevInit
			self.m_io = ioInit
			self.m_b = bInit
			self.m_ilim = ilimdiffusionInit
			self.m_contribution=0
			#assume reaction is succesful
			self.m_status=True
			#print("Reaction ", curvename," loaded succesfully!")
		except:
			#print("Loading reaction ",self.m_curvename," failed!")
			self.m_status=False

#
#This function returns the cathodic current density at a given potential
#

	def get_currentdensity(self, potential, ScaleFactor, FeSfactor):
		#The Tafel portion of cathodic current density
		try:
			# Tafel current density,A/m2
			itaf = self.m_io * 10 ** ((self.m_Erev - potential) / self.m_b)
			#Tafel current density modified by FeCO3 and FeS film,A/m2
			itaf = itaf * ScaleFactor * FeSfactor
			#total cathodic current density,A/m2
			return 1 / (1 / itaf + 1 / self.m_ilim)
			#print("Reaction ", self.m_curvename, " Calculation succesfully completed!")
		except:
			#print("Reaction ",self.m_curvename," Calculation Fails!")
			self.m_status=False
			return 0

#
#This function returns an array of cathodic current density for any given potential
#

	def get_currentdensity_array(self, minPotential, maxPotential, Npoints, ScaleFactor, FeSfactor):

	#counter for number of data points
		IncOvp = abs(maxPotential - minPotential) / Npoints
	#Calculates the cathodic current density at every given potential,A/m2
		retrunArray=[]
		for i  in range(Npoints):
			retrunArray.append(self.get_currentdensity(maxPotential - (i - 1) * IncOvp, ScaleFactor, FeSfactor))
		return retrunArray