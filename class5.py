import moose
import matplotlib.pyplot as pyplot
import numpy

def load_neuron_file(fileName, cellPath, Rm, Ra, Cm, initVm):
    cell = moose.loadModel(fileName, cellPath)
    for child in cell[0].children:
        for comp in child:
            comp.Rm = Rm
            comp.Ra = Ra
            comp.Cm = Cm
            comp.initVm = initVm

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
    pyplot.clf() # clear figure
    for table in tables:
        x = numpy.fromiter((i*table.dt for i in range(table.size)), float, table.size)
        y = table.vector
        pyplot.plot(x, y)
def run_sim():
    moose.reinit()
    moose.start(900e-3)

neuron = moose.Neutral("/neuron")
inputs = moose.Neutral("/inputs")
outputs = moose.Neutral("/outputs")

swcNeuron = load_neuron_file("538ser3sl5-cell1-2-a.CNG.swc", "/swcNeuron", 0.8, 0.9, 0.01, 0.072)
moose.le("/swcNeuron")

pulse = create_pulse("/inputs/somaPulse",10e-3,10e-3,0.1e-9,moose.element("/swcNeuron/soma"))
table = create_table('outputs/somaVmTable',moose.element("/swcNeuron/soma"),"getVm")
dtable = create_table('outputs/dendVmTable',moose.element("/swcNeuron/dend_3_0"),"getVm")

run_sim()

plot_tables([table, dtable])
pyplot.show()

ch = moose.HHChannel("ch")
moose.showfield(ch)
ch.Xpower = 3
ch.Ypower = 1
ch.Ek = 0.05
ch.Gbar = 1

