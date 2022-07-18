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

#
#This class module deals with anodic reactions.
#
import configparser
import constant
import math
config = configparser.ConfigParser()
config.read('model.ini')

class anodic:
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
#_Epit:pitting potential,V
#m_PointModel:variable for PointModel
#m_contribution:contribution of species to overall corrosion rate
#Constant:
#******************************************************************************

	#This subroutine initializes some common parameters used in
	#further calculations, such as name and type of the reaction,
	#reversible potential,Tafel slope,exchange current density,
	#limiting current density,passivation potential,
	#passive current density and pitting potential,etc.
	#Reference exchange current density and activation energy
	#are read from the file,other parameters are calculated.
   
	def __init__(self, curvename, temperature):
		#reversible potential,exchange current density,Tafel slope,limiting current density,respectively
		#ErevInit, ioInit, bInit, ilimdiffusionInit
		#passivation potential,passive current density,pitting potential,respectively
		#EpasInit, ipasInit, EpitInit
		#reference exchange current density and activation energy,respectively
		#ioref, ActivEnergy
		#data read from the file as a string
		self.m_curvename=curvename
		self.m_status=True
		try:
			if curvename== "Iron_dissolution":
				ioref = float(config["Iron_dissolution"]["ioref"])
				ActivEnergy = float(config["Iron_dissolution"]["ActivEnergy"])
				ioInit = ioref* 10 ** (4 - 4)* math.exp(-ActivEnergy / constant.r * (1 / (temperature) - 1 / 298.16)) #A/m^2		
				ErevInit = -0.488
				bInit = 2.303 * constant.r * (temperature) / constant.F / 1.5
				EpasInit = 0
				ipasInit = 0
				EpitInit = 0
				#if EpasInit == 0:
				#	EpasInit = 1000
				#if EpitInit == 0:
				#	EpitInit = 1000
			else:
			 #for user-defined reactions,read all the data from the file
				if float(config["User_Anodic_Reaction"]["Temperature(K)"]) != temperature:
					print("Reaction ", curvename, " is not compatible!")
					return
				ioInit = float(config["User_Anodic_Reaction"]["ExchangeCurrentDensity(A/m2)"])
				ErevInit = float(config["User_Anodic_Reaction"]["ReversablePotential(V)"])
				bInit = float(config["User_Anodic_Reaction"]["TafelSlope(b)"])
				ilimdiffusionInit = float(config["User_Anodic_Reaction"]["LimitingCurrentDensity(A/m2)"])
				EpasInit = float(config["User_Anodic_Reaction"]["PassivationPotential(V)"])
				ipasInit = float(config["User_Anodic_Reaction"]["PassivationCurrentDensity(A/m2)"])
				EpitInit = float(config["User_Anodic_Reaction"]["PittingPotential(V)"])		
			self.m_Erev = ErevInit
			self.m_io = ioInit
			self.m_b = bInit
			#self.m_ilim = ilimdiffusionInit
			self.m_Epas = EpasInit
			self.m_ipas = ipasInit
			self.m_Epit = EpitInit
			self.contribution=0
		except:
			#print("Loading Reaction ",curvename, " failed!")
			self.m_status=False
			return
		  
	#This function returns the anodic current density for a given potential
	   
	def get_currentdensity(self, potential, ScaleFactor, FeSfactor):
		try:
			#Sets up constant to construct a smooth anodic polarization curve
			#with active,passive and tran-passive zone.
			Apas = 50000
			pas = 2.5
			pit = 2.5
			#declare some temporary variables for construction of anodic polarization curve
			# Calculates anodic current density, A/m2
			itaf = self.m_io * 10 ** ((potential - self.m_Erev) / self.m_b) #current density in the active zone
			if potential > self.m_Epas:
				Spas = math.exp(-Apas * (potential - self.m_Epas) ** pas) #be careful with exp fucntion
			else:
				Spas = 1
			if potential > self.m_Epit:
				ipit = self.m_ipas * 10 ** ((potential - self.m_Epit) / self.m_b) #current density in the tran-passive zone
				Apit = 1 / self.m_b
				Xpit = math.exp(Apit * (potential - self.m_Epit) ** pit)
				Spit = Xpit / (1 + Xpit)
			else:
				Spit = 1
			ipit=0
			#universal equation for calculating current density in the whole range
			tempo0 = Spit * (Spas * itaf + (1 - Spas) * self.m_ipas) + (1 - Spit) * ipit
			#current density is modified by the presence of FeCO3 and Fes film
			tempo = tempo0 * ScaleFactor * FeSfactor
			return tempo
		except:
			self.m_status=False
			#print("Reaction ", self.m_curvename," Calculation Fails!")
			return 0
	   
	#
	#This function generates an array of anodic current density based on given potentials
	#Used for polarization curve plotting
	#

	def get_currentdensity_array(self, minPotential, maxPotential, Npoints, ScaleFactor, FeSfactor):

		IncOvp = abs(maxPotential - minPotential) / Npoints

		#Calculates the anodic current density at every given potential,A/m2
		retrunArray=[]
		for i in range(1,Npoints):
			retrunArray.append(self.getCurrentDensity(minPotential + (i - 1) * IncOvp, ScaleFactor, FeSfactor))
			
		return retrunArray
