import moose


simtime = 0.1
simdt = 0.25e-5
plotdt = 0.25e-3
for i in range(10):
    moose.setClock(i, simdt)
moose.setClock(8, plotdt)

spike = moose.RandSpike("/spike0")
spike.rate = 1300
spike.refractT = 1e-3

tab = moose.Table("/pre0")
moose.connect(spike,"spikeOut",tab,"spike")

moose.reinit()
moose.start(simtime)

print(tab.vector)
