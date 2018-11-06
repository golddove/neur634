import moose
import numpy
import sys
import matplotlib.pyplot as pyplot

def create_compartment(name,length,diameter,surfaceArea,xArea,membraneResistivity,axialResistivity,capacitivity): # Em, initVm?
    newC = moose.Compartment(name)
    newC.length = length
    newC.diameter = diameter
    #newC.Em = Em
    #newC.initVm = initVm
    newC.Rm = float(membraneResistivity) / surfaceArea
    newC.Cm = capacitivity * surfaceArea
    newC.Ra = float(axialResistivity * length) / xArea
    return newC

def create_spherical_compartment(name,length,diameter,membraneResistivity,axialResistivity,capacitivity):
    surfaceArea = numpy.pi*diameter**2
    xArea = float(surfaceArea) / 4 
    return create_compartment(name,length,diameter,surfaceArea,xArea,membraneResistivity,axialResistivity,capacitivity)

def create_pulse(pulsename,pulsedelay,pulsewidth,pulselevel,pulsecomp):
    pulse = moose.PulseGen(pulsename)
    pulse.delay[0]=pulsedelay
    pulse.width[0]=pulsewidth
    pulse.level[0]=pulselevel
    pulse.delay[1]=1e9
    moose.connect(pulse,"output",pulsecomp,"injectMsg")
    return pulse

def create_table(tablename, tablecomp, compproperty):
    tab = moose.Table(tablename)
    moose.connect(tab,"requestOut",tablecomp,compproperty)
    return tab

def plot_table(table):
    x = numpy.fromiter((i*table.dt for i in range(table.size)), float, table.size)
    y = table.vector
    pyplot.plot(x, y)

def main():
    neuron = moose.Neutral("/neuron")
    inputs = moose.Neutral("/inputs")
    outputs = moose.Neutral("/outputs")

    soma = create_spherical_compartment("/neuron/soma",25e-6,20e-6,2.8,4.0,0.03)
    pulse = create_pulse('/inputs/somaPulse',50e-3,100e-3,1e-9,soma)
    vmtable = create_table('outputs/somaVmTable',soma,"getVm")
    print(soma.tick)
    print(soma.dt)
    moose.reinit()
    moose.start(300e-3)
    moose.showfields(vmtable)
    moose.showmsg(soma)
    moose.showmsg(pulse)
    moose.showmsg(vmtable)

    plot_table(vmtable)
    pyplot.show()

    #t = pylab.linspace(0, 300e-3, len(vmtable.vector))
    #pylab.plot(t,)

main()
