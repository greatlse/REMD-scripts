#! /bin/env python
# -*- coding: UTF-8 -*-

'''
#==========================================================
#     Filename: remd_temp_calculator.py
#         Desc: generate mdp file for REMD with GMX
#       Author: senliu
#      Created: 2015-Oct-17-14:51
#   LastChange: 
#      History: 
#          Ref:http://githup.com/mchelem/remd_temperature
#==========================================================
'''

import logging
import requests
import BeautifulSoup
from sys import argv

'''
#========================USAGE===============================
#    Usage: python remd_temp_calculator.py parafile num
# parafile: the parameters for retrieving data from webserver
#      num: the number temp item that you want
#=========================END==============================
'''

class mdpfilegen(object):
	__esol_mdpfile='''
define =
; Run paramaters
integrator = md ; leap-frog integrator
dt = 0.002 ; 2fs
nsteps = 50000 ; 0.002 * 50000 = 100ps
; Output control
nstxout=0                                
nstvout=0
nstenergy=100
nstlog=100
nstxtcout=500
xtc-grps=Protein
energygrps=Protein Water_and_ions
; neighborsearch param
nstlist = 10
ns_type=grid
pbc = xyz
rlist=1.0
rlistlong=-1
; options for elec and vdw
coulombtype = PME
rcoulomb = 1.0
epsilon_r=1
epsilon_rf=1
vdw-type=Cut-off
rvdw=1.2
fourierspacing=0.12
pme_order =4
ewald_rtol = 1e-05
optimize_fft = yes
; Options for weak coupling algorithms
tcoupl = v-rescale
ld_seed=2015
tc_grps= Protein Water_and_ions
tau_t =0.1 0.1
ref_t= 300 300
Pcoupl=No
Pcoupltype=Isotropic
nstpcouple = -1
tau-p=1
ref-p=1
compressibility= 4.5e-5
;simulated anneling
;annealing = single single
;annealing_npoints= 3 3
;anealing_time = 0 50 100 0 50 100
;annealing_temp = 100 200 300 100 200 300
; generate velocities for startup run
gen_vel = yes
gen_temp = 300
gen_seed= 399
; options for bond
constraints = hbonds
constraint_algorithm=lincs
continuation = no
lincs_order=4
lincs_iter = 1
'''
	def __init__(self, tlist):
		self.tlist=tlist
		esol=self.__esol_mdpfile.split("\n")
		self.esol=esol[1:-1]
		self.time=0.1
	
	def __listwrite(self, flist, fname):
		f=file(fname, 'w')
		f.writelines(flist)
		f.close()

	def treplace(self, pname, tm):
		# pname means that a item in mdp file for replace, such as
		# "ref_t" or "gen_temp"
		repart=['ref_t','gen_temp']
		if pname in repart:
			for index, line in enumerate(self.esol):
				if pname in line:
					tindex=index
			tm=round(tm,2)
			self.esol[tindex]=pname + " = " + str(tm) + " " + str(tm)
				#fname="heat_" + str(i) + ".mdp"
				#tmp="\n".join(self.esol)
				#self.__listwrite(tmp, fname)
			return self.esol
		else:
			err='''
			Nothing can be replaced !
			Please confirm that the item is either "ref_t" or "gen_temp"
			'''
			raise TypeError(err)

	def equilmdpout(self):
		for i, tm in enumerate(self.tlist):
			self.treplace("ref_t", tm)
			self.treplace("gen_temp", tm)
			fname="heat_" + str(i) + ".mdp"
			tmp="\n".join(self.esol)
			self.__listwrite(tmp, fname)

	def preplace(self, tm):
		self.__esol_mdpfile=self.__esol_mdpfile.replace("Pcoupl=No", "Pcoupl= berendsen")
		self.__esol_mdpfile=self.__esol_mdpfile.replace("continuation = no", "continuation = yes")
		self.__esol_mdpfile=self.__esol_mdpfile.replace("gen_vel = yes", "gen_vel = no")
		esol=self.__esol_mdpfile.split("\n")
		self.esol=esol[1:-1]
		self.treplace("ref_t", tm)
		self.treplace("gen_temp", tm)
		self.step_replace(self.time)
	
	def mdmdpout(self):
		for i, tm in enumerate(self.tlist):
			self.preplace(tm)
			fname="prod_" + str(i) +".mdp"
			tmp="\n".join(self.esol)
			self.__listwrite(tmp, fname)
	
	def step_replace(self, time):
		# time : ns
		self.time=time
		dtitem=[ d for d in self.esol if d.startswith("dt")]
		dt=float(dtitem[0].split(";")[0].split("=")[1])
		step=int(float(time)*1000/dt)
		stem_item=[ item for item in self.esol if item.startswith("nsteps")]
		sindex=self.esol.index(stem_item[0])
		self.esol[sindex]="nsteps = " + str(step)



class temp_requst(object):
	default_params=dict(
		Pdes=0.25,
		Tol=1e-4,
		Tlow=300,
		Thign=500,
		Nw=0,
		WC=3,
		Np=200,
		PC=1,
		Hff=0,
		Vs=0,
		Alg=0,
		)
	constraints_water = {
			'fully flexible':0,
			'flexible angle':2,
			'rigid':3,
			}
	constraints_protein = {
			'fully flexible': 0,
			'bonds to hydrogens only': 1,
			'all bonds': 2,
			}
	params_mapping = {
			'exchange probability': 'Pdes',
			'lower temperature limit': 'Tlow',
			'number of water molecules': 'Nw',
			'number of protein atoms': 'Np',
			'hydrogens in protein': ('Hff', {'all h':0, 'polar h':1}),
			'tolerance': 'Tol',
			'upper temperature limit': 'Thign',
			'constraints in water': ('WC', constraints_water),
			'constraints in the protein': ('PC', constraints_protein),
			'virtual sites in protein': ('Vs', {'none':0, 'virtual hydrogen':1}),
			}

	def __init__(self, input_param):
		''' Load parameters from a dictionary
		The parameters can be specified both as the original from parameters
		on the server or as the labels for these parameters
		Original form parameters:
			param['PC'] = 2 # Constraints in the protein: All bonds
		Human readable patameters:
			param['constraints in the protein'] = 'all bonds'
		'''
		self.logger = logging.getLogger(__name__)
		self.TGENERATOR_URL = 'http://folding.bmc.uu.se/remd/tgenerator.php'

		keys = [key.lower() for key in input_param.keys()]
		if not ('np' in keys or 'number of protein atoms' in keys):
			raise AttributeError("You must specify the number of protein atoms")

		self.param = self.default_params.copy()
		if 'Np' in input_param:
			self.param.update(input_param)
		else:
			for key, value in input_param.iteritems():
				try:
					#Remove spaces and make everything lower case
					key = ' '.join(key.lower().split())
					mapping = params_mapping[key]
					if isinstance(value, str):
						value=' '.join(value.lower().split())
					if isinstance(mapping, tuple):
						self.param[mapping[0]]= mapping[1][value]
					else:
						self.param[mapping] = value
				except KeyError as e:
					logger.warn('Invalid variable: ' + str(e))
	def get_temeratures(self):
		'''
		Retrieve the temperatures (K) as a list
		Sample Output:[300.0, 332.18, 366.98, 404.54, 445.26, 489.27, 536.92]
		'''
		temperatures = []
		res = requests.post(self.TGENERATOR_URL, data=self.param)
		soup = BeautifulSoup.BeautifulSoup(res.text)
		heading = soup.findAll(
				lambda x: 'we also give the temperatures below' in x.text)
		temperatures = heading[0].nextSibling.nextSibling.text
		
		return [float(x) for x in temperatures.split(',')]

	def __read_table(self,table):
		'''
		Read an HTML table filled with numbers to a list
		'''
		result = []
		for row in table.findAll('tr'):
			cols = row.findAll('td')
			if cols and len(cols) > 1:
				result.append([])
				for col in cols[1:]:
					test = ''.join(col.findAll(test=True))
					test = 'nan' if not text else text
					result[-1].append(float(text))
		return result

	def get_temeratures_energies(self):
		'''
		Get a table containing the following parameters
		The parameters, from left to right are the following:
		Temperature (K)
		μ(kJ/mol)
		σ(kJ/mol)
		μ12(kJ/mol)
		σ12(kJ/mol)
		P12
		'''
		res = requests.post(self.TGENERATOR_URL, data=self.param)
		soup = BeautifulSoup.BeautifulSoup(res.text)
		table = soup.findAll('table')[1]

		return self.__read_table(table)

def templist(pafile):
	'''
	Get the temperatures list from web server based on the parameters file
	'''
	para={}
	fp=open(pafile,'r')
	palist=fp.readlines()
	fp.close()

	for line in palist:
		pline=line.strip('\n').split(':')
		if not ( '' in pline):
			para[pline[0].strip()]=pline[1].strip()
	tr=temp_requst(para)
	tlist=tr.get_temeratures()
	return tlist

def tlist_control(num, tlist):
	if isinstance(num, str):
		num=int(num)
	tlist=tlist[:-1:2]
	tlist=tlist[:num]

	return tlist

def main():
	ff=argv[1]
	tlist=templist(ff)
	num=argv[2]
	tlist=tlist_control(num, tlist)
	first=mdpfilegen(tlist)
	#first.step_replace(100)
	first.equilmdpout()
	first.mdmdpout()

if __name__ == "__main__":
	main()

