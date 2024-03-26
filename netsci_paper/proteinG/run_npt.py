"""
Run a NPT (constant number of atoms, pressure, and temperature)
simulation of protein G identical to the simulations run in:

Lange OF, Grubmüller H. Generalized correlation for biomolecular 
dynamics. Proteins. 2006 Mar 1;62(4):1053-61. 
doi: 10.1002/prot.20784. PMID: 16355416.
"""

import simtk.openmm as mm
from sys import stdout
from simtk.openmm import app
from simtk import unit
import os

gro_input = '1pgb_solv_ions.gro'
gro_top = 'topol.top'
# Load the system from GROMACS files
gro = app.GromacsGroFile(gro_input)
top = app.GromacsTopFile(gro_top, periodicBoxVectors=gro.getPeriodicBoxVectors())

# Set up the system
system = top.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometers, constraints=app.AllBonds)

# Set up the integrator
temperature = 300*unit.kelvin
pressure = 1*unit.atmospheres
integrator = mm.LangevinIntegrator(temperature, 0.1/unit.picoseconds, 2.0*unit.femtoseconds)

# Set up the barostat
barostat = mm.MonteCarloBarostat(pressure, temperature, 50)
system.addForce(barostat)

# Create the simulation context
platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed'}
simulation = app.Simulation(top.topology, system, integrator, platform, properties)

# Set the initial positions and velocities
simulation.context.setPositions(gro.positions)
# Equilibration
simulation.minimizeEnergy()
simulation.reporters.append(app.StateDataReporter(stdout, reportInterval=1000, step=True, potentialEnergy=True, temperature=True))
simulation.step(2500000)
print("equilibration done")

# Production run
#simulation.reporters.append(app.StateDataReporter(stdout, reportInterval=1000, step=True, potentialEnergy=True, temperature=True, volume=True))
simulation.reporters.append(app.DCDReporter('1pgb_out.dcd', 5000))
simulation.step(97500000)

