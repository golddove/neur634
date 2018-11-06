# Interpreter Python 2.7
# Author: Risheek Rajolu
# Date: September 4, 2018

import moose
import numpy
import sys

moose.__version__


neuron = moose.Neutral('neuron')
moose.Compartment('neuron/soma')
soma = moose.element('neuron/soma')
dend = moose.Compartment('neuron/dend')

soma.length = 25e-6
soma.diameter = 20e-6
soma.Em = -65e-3
soma.initVm = -65e-3
sA = numpy.pi*soma.diameter**2 #Spherical surface area
resistivity = 200e-3 # Converted from 2000 ohm-cm^2
soma.Rm = resistivity * 1.0 / sA #Resistance
capacitivity = 10e-3 # Converted from 1 uF-cm^2
soma.Cm = capacitivity * sA #Capacitance

print(neuron)
print(soma)
print("l=" + str(soma.length))
print("d=" + str(soma.diameter))
print("Rm=" + str(soma.Rm))
print("Vm=" + str(soma.Vm))
print("Cm=" + str(soma.Cm))
print("Em=" + str(soma.Em))
print("initVm=" + str(soma.initVm))
print(dend)

moose.reinit()
moose.start(0.3)
print(soma.Vm)
soma.inject=1e-9
moose.reinit()
moose.start(0.3)
print(soma.Vm)
soma.inject = 0

def create_compartment(name):
    pass

def create_pulse(pulsename,pulsedelay,pulsewidth,pulselevel,pulsecomp):
	pulse = moose.PulseGen(pulsename)
	pulse.delay[0]=pulsedelay
	pulse.width[0]=pulsewidth
	pulse.level[0]=pulselevel
	pulse.level[1]=1e-9
	pulse.delay[1]=1e9
	moose.connect(pulse,'output',pulsecomp,'injectMsg')
pulse = create_pulse('pulse',0,0,0,soma)
print(moose.showmsg(soma))

def main(*args):
    pass

main(*sys.argv)
