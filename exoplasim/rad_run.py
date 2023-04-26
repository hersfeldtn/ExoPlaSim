#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Mo Cohen
"""

import exoplasim as exo 
import os

run_length = 1
top_dir = '/home/s1144983/Repos/exodev/exoplasim/'
density = 1262
radius = 500e-09
planet = 'wolf'

wolfdict = {'radius' : 1.66, 'solcon' : 1777, 'startemp' : 3408, 'N2' : 0.988,
            'CO2' : 0.01, 'CH4': 0.002, 'H2' : 0.0, 'rotperiod' : 17.9, 'gravity' : 12.1,
            'starrad' : 0.32, 'starspec' : top_dir + 'stellarspectra/wolf.dat',
            'aerofile' : top_dir + 'hazeconstants/wolf_constants_500.dat',
            'name' : 'wolf', 'eccentricity': 0.0, 'bulk' : 1}
            
current_dir = os.getcwd()
print('Current directory is ' + str(current_dir))
if current_dir != top_dir:
os.chdir(top_dir)
new_dir = os.getcwd()

print('Current directory is ' + str(new_dir))
print('Simulating ' + planet)
print('Particle density: ' + str(density))
print('Particle radius: ' + str(radius))

if planet == wolfdict['name']:
config = wolfdict
else:
print('No config dictionary available for ' + str(planet))


model = exo.Model(workdir = planet +'_' + str(density) + '_' + str(radius),
                  modelname = planet +'_' + str(density) + '_' + str(radius),
                  resolution='T21',
                  layers=10,
                  ncpus=16,
                  precision=8,
                  outputtype='.npz')

model.configure(flux = config['solcon'], startemp = config['startemp'],
               starspec = config['starspec'],
               pH2 = config['H2'], pHe=0.0, pN2 = config['N2'],
               pO2 = 0.0, pCO2 = config['CO2'], pAr = 0.0, pNe = 0.0,
               pKr = 0.0, pH2O = 0.0, pCH4=config['CH4'], rotationperiod = config['rotperiod'],
               synchronous=True, substellarlon=180.0, gravity = config['gravity'],
               radius = config['radius'], eccentricity = config['eccentricity'], 
               obliquity = 0.0, aquaplanet = True, stratosphere = True,
               timestep = 15.0, snapshots = None,
               otherargs={'NQSPEC@plasim_namelist':'1', 'NLOWIO@plasim_namelist':'1'},
               aerosol=True, aerobulk = config['bulk'], asource=1,
               fcoeff = 1e-07, apart = radius, rhop = density, aerofile=config['aerofile'])

model.exportcfg()
print('Starting run')

model.run(years = run_length)
print('Finalising run')

model.finalize(planet +'_' + str(density) + '_' + str(radius) + '_output')

