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
    plot_tables([table])

def plot_tables(tables):
    #pyplot.clf() # clear figure
    for table in tables:
        x = numpy.fromiter((i*table.dt for i in range(table.size)), float, table.size)
        y = table.vector
        pyplot.plot(x, y)

def run_sim():
    moose.reinit()
    moose.start(900e-3)

def main():
    neuron = moose.Neutral("/neuron")
    inputs = moose.Neutral("/inputs")
    outputs = moose.Neutral("/outputs")
    
    soma = create_spherical_compartment("/neuron/soma",25e-6,20e-6,2.8,4.0,0.01)
    pulse = create_pulse("/inputs/somaPulse",50e-3,100e-3,1e-9,soma)
    vmtable = create_table('outputs/somaVmTable',soma,"getVm")

    for di, simDt in enumerate([5e-5, 1e-3, 5e-3]):
        for ci, dendCount in enumerate([1, 5]):
            totalLength = 100e-6
            dendVmtable = None
            dends = []
            for i in range(dendCount):
                dend = create_spherical_compartment("/neuron/dendrite" + str(i),totalLength/dendCount,2e-6,2.8,4.0,0.01)
                moose.setClock(dend.tick, simDt)
                prev = soma if i==0 else dends[i-1]
                moose.connect(prev,"axialOut",dend,"handleAxial")
                dends.append(dend)
                if i == dendCount//2:
                    dendVmtable = create_table('outputs/dendVmTable'+str(di*2+ci),dend,"getVm")
            moose.setClock(soma.tick, simDt)
            run_sim()
            print(vmtable)
            print(dendVmtable)
            plot_tables([vmtable, dendVmtable])
    
    pyplot.show()

main()
