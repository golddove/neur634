# Interpreter Python 2.7
# Author: Risheek Rajolu
# Date: September 4, 2018

import moose
import numpy

#moose.__version__

#neuron = moose.Neutral('neuron')
#moose.Compartment('neuron/soma')
#soma = moose.element('neuron/soma')
#dend = moose.Compartment('neuron/dend')

#soma.length = 25e-6
#soma.diameter = 20e-6
#soma.Em = -65e-3
#soma.initVm = -65e-3
#sA = numpy.pi*soma.diameter**2 #Spherical surface area
#resistivity = 200e-3 # Converted from 2000 ohm-cm^2
#soma.Rm = resistivity * 1.0 / sA #Resistance
#capacitivity = 10e-3 # Converted from 1 uF-cm^2
#soma.Cm = capacitivity * sA #Capacitance

#print(neuron)
#print(soma)
#print("l=" + str(soma.length))
#print("d=" + str(soma.diameter))
#print("Rm=" + str(soma.Rm))
#print("Vm=" + str(soma.Vm))
#print("Cm=" + str(soma.Cm))
#print("Em=" + str(soma.Em))
#print("initVm=" + str(soma.initVm))
#print(dend)

#moose.reinit()
#moose.start(0.3)
#print(soma.Vm)
#soma.inject=1e-9
#moose.reinit()
#moose.start(0.3)
#print(soma.Vm)
#soma.inject = 0

def create_compartment(name,length,diameter,surfaceArea,xArea,membranceResistivity,axialResistivity,capacitivity): # Em, initVm?
	newC = moose.Compartment(name)
	newC.length = length
	newC.diameter = diameter
	newC.Em = Em
	newC.initVm = initVm
	newC.Rm = resistivity * 1.0 / surfaceArea
	newC.Cm = capacitivity * surfaceArea
	newC.Ra = axialResistivity * length * 1.0 / xArea

def create_spherical_compartment(name,length,diameter,membraneResistivity,axialResistivity,capacitivity):
	surfaceArea = numpy.pi*diameter**2
	xArea = surfaceArea * 1.0 /4 
	create_compartment(name,length,diameter,surfaceArea,xArea,membraneResistivity,axialResistivity,capacitivity)

def create_pulse(pulsename,pulsedelay,pulsewidth,pulselevel,pulsecomp):
	pulse = moose.PulseGen(pulsename)
	pulse.delay[0]=pulsedelay
	pulse.width[0]=pulsewidth
	pulse.level[0]=pulselevel
	pulse.level[1]=1e-9
	pulse.delay[1]=1e9
	moose.connect(pulse,'output',pulsecomp,'injectMsg')

#output = moose.Neutral("/output")
#neuron = moose.Neutral("/neuron")
# RA = 4 ohm-meters
#soma = create_spherical_compartment("neuron/soma",20e-6,20e-6,2,4,10e-3) # -65e-3, -65e-3
#dend = create_spherical_compartment("neuron/dend",40e-6,8e-6,2,4,10e-3) # -65e-3, -65e-3
#pulse = create_pulse('pulse',0,0,0,soma)

#vmtab = moose.Table("/output/somaVm")
#moose.showfield(vmtab)

#moose.connect(vmtab, "requestOut", soma, "getVm")

#moose.setClock(soma.tick, 2e-5)
#moose.showmsg(soma)
#moose.reinit()
#moose.start(900e-3)

#import pylab
#t = pylab.linespace(0,900e-3,len(vmtab.vector))
#pylab.plot(t,vmtab.vector)
#pylab.show()

def read_genesis_file(filename, name):
	moose.loadModel(filename, name)

def read_swc_file(filename, parentName):
	moose.loadModel(filename, parentName)
	for component in moose.element(parentName).children():
		pass

import sys

pfile = "layer2.p"
swcfile = "538ser3sl5-cell1-2-a.CNG.swc"
container1="cell1"
container2="cell2"
cell1 = moose.loadModel(pfile,container1)
cell2 = moose.loadModel(swcfile,container2)
print("P file:")
moose.le("/cell1")
print("SWC file:")
moose.le("/cell2")
moose.showfield("/cell2/soma")
moose.showfield("/cell2/soma_1_0")
moose.showfield("/cell2/soma_2_0")

moose.showmsg("cell1/soma")
moose.showmsg("cell2/soma")

#import neuron
from neuron import h, gui

print("\n\n\n DOING NEURON NOW! \n\n\n")

soma = h.Section(name="soma")
dir(soma)
soma.insert("pas")
soma.nseg=3
h.psection(sec=soma)

dend = h.Section(name="dend")
dend.connect(soma,1)
dend.nseg = 10

print(h.topology())

soma.L = soma.diam = 12.6157
dend.L = 200
dend.diam = 1
h.area(0.1, sec = dend)

pas = soma(0.5).pas
pas.g=0.002
print("pas: " + str(soma(0.5).g_pas))
h.psection(sec=dend)

def printAll():
    print("\nALL SECTIONS:")
    for sec in h.allsec():
        print(sec.name())
        for seg in sec:
            print("\t" + str(seg.x))
            for mech in seg:
                print("\t\t" + mech.name())
    print("\n")

printAll()

dend.insert("pas")

stim = h.IClamp(soma(0.5))
stim.delay = 5
stim.dur = 1
stim.amp = 0.1
h.psection(sec=soma)

v_vec = h.Vector()
t_vec = h.Vector()
v_vec.record(soma(0.5)._ref_v)
t_vec.record(h._ref_t)

h.tstop = 40.0
h.run()

from matplotlib import pyplot
pyplot.figure(figsize=(8,4))
pyplot.plot(t_vec,v_vec)
pyplot.xlabel("time (ms)")
pyplot.ylabel("mV")
pyplot.show()

