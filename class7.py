# ionchannel.py ---
#
# Filename: ionchannel.py
# Description:
# Author: Subhasis Ray
# Maintainer:
# Created: Wed Sep 17 10:33:20 2014 (+0530)
# Version:
# Last-Updated:
#           By:
#     Update #: 0
# URL:
# Keywords:
# Compatibility:
#
#

# Commentary:
#
#
#
#

# Change log:
#
#
#
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street, Fifth
# Floor, Boston, MA 02110-1301, USA.
#
#

# Code:
import numpy as np
import matplotlib.pyplot as plt
import moose
import collections

EREST_ACT = -70e-3 #: Resting membrane potential


AbParams = collections.namedtuple("AbParams",[
    "aRate",
    "aB",
    "aC",
    "aVHalf",
    "aVSlope",
    "bRate",
    "bB",
    "bC",
    "bVHalf",
    "bVSlope"
    ])

CaDepParams = collections.namedtuple("CaDepParams", [
    "Kd",
    "power",
    "tau"
    ])

ChannelSettings = collections.namedtuple("ChannelSettings",[
    "xPower",
    "yPower",
    "zPower",
    "eRev",
    "name",
    "xParam",
    "yParam",
    "zParam",
    "chan_type"
    ])

SK_Z_params=CaDepParams(0.57e-3, 5.2, 4.9e-3)

#: The parameters for defining m as a function of Vm
Na_m_params = AbParams(
        1e5 * (25e-3 + EREST_ACT),  # 'A_A':
        -1e5,                       # 'A_B':
        -1.0,                       # 'A_C':
        -25e-3 - EREST_ACT,         # 'A_D':
        -10e-3,                     # 'A_F':
        4e3,                        # 'B_A':
        0.0,                        # 'B_B':
        0.0,                        # 'B_C':
        0.0 - EREST_ACT,            # 'B_D':
        18e-3                       # 'B_F':
        )

#: Parameters for defining h gate of Na+ channel
Na_h_params = AbParams(
        70.0,                        # 'A_A':
        0.0,                       # 'A_B':
        0.0,                       # 'A_C':
        0.0 - EREST_ACT,           # 'A_D':
        0.02,                     # 'A_F':
        1000.0,                       # 'B_A':
        0.0,                       # 'B_B':
        1.0,                       # 'B_C':
        -30e-3 - EREST_ACT,        # 'B_D':
        -0.01                    # 'B_F':
        )



#: K+ channel in Hodgkin-Huxley model has only one gate, n and these
#are the parameters for the same
K_n_params = AbParams(
        1e4 * (10e-3 + EREST_ACT),   #  'A_A':
        -1e4,                      #  'A_B':
        -1.0,                       #  'A_C':
        -10e-3 - EREST_ACT,         #  'A_D':
        -10e-3,                     #  'A_F':
        0.125e3,                   #  'B_A':
        0.0,                        #  'B_B':
        0.0,                        #  'B_C':
        0.0 - EREST_ACT,            #  'B_D':
        80e-3                       #  'B_F':
        )

KDr_X_params = AbParams(
        28.2,
        0,
        0.0,
        0,
        -12.5e-3,
        6.78,
        0.0,
        0.0,
        0.0,
        33.5e-3)

CaL_X_params = AbParams(
-880,
-220e3,
-1.0,
4.0003e-3,
-7.5e-3,
-284,
71e3,
-1.0,
-4.0003e-3,
5e-3)

PoolSettings = collections.namedtuple("CaParamTuple", [
"CaBasal",
"CaThick",
"CaTau",
"BufCapacity",
"name"
])

Ca_pool_settings = PoolSettings(50e-6, 1e-6, 20e-3, 20, "CaPool")

Na_settings = ChannelSettings(3,1,0,0.06,"na",Na_m_params,Na_h_params,None,"")

K_settings = ChannelSettings(4,0,0,-0.1,"k",K_n_params,None,None,"")

KDr_settings = ChannelSettings(2,0,0,90e-3,"kdr",KDr_X_params,None,None,"")

SK_settings = ChannelSettings(0,0,1,-87e-3,"skca",None,None,SK_Z_params,"ca_dependent")

CaL_settings = ChannelSettings(1,0,0,130e-3,CaL_X_params,None,None,"CaL","ca_permeable")

#: We define the rate parameters, which are functions of Vm as
#: interpolation tables looked up by membrane potential.
#: Minimum x-value for the interpolation table
VMIN = -30e-3 + EREST_ACT
#: Maximum x-value for the interpolation table
VMAX = 120e-3 + EREST_ACT
#: Number of divisions in the interpolation table
VDIVS = 3000

def create_channel_proto(channelSettings, vDivs = 3000, vMin = -30e-3 + EREST_ACT, vMax = 120e-3 + EREST_ACT, caDivs=10000, caMin=0, caMax=1):
    if(not(moose.exists("/library"))):
        moose.Neutral('/library')
    channel = moose.HHChannel('/library/' + channelSettings.name)
    channel.tick = -1
    if(channelSettings.xPower > 0):
        channel.Xpower = channelSettings.xPower
        xGate = moose.HHGate(channel.path + "/gateX")
        xGate.setupAlpha(channelSettings.xParam + (vDivs, vMin, vMax))
    if(channelSettings.yPower > 0):
        channel.Ypower = channelSettings.yPower
        yGate = moose.HHGate(channel.path + "/gateY")
        yGate.setupAlpha(channelSettings.yParam + (vDivs, vMin, vMax))
    if(channelSettings.zPower > 0):
        channel.Zpower = channelSettings.zPower
        zGate = moose.HHGate(channel.path + "/gateZ")
        zGate.min = caMin
        zGate.max = caMax
        caTerm = (np.linspace(caMin, caMax, caDivs)/channelSettings.zParam.Kd)**channelSettings.zParam.power
        inf_z=caTerm/(1+caTerm) # open probability at steady state for a given Ca concentration
        tau_z=channelSettings.zParam.tau*np.ones(caDivs) # constant tau for any Ca concentration
        zGate.tableA = inf_z/tau_z #alpha
        zGate.tableB = 1/tau_z #alpha + beta?
        channel.useConcentration = True
    return channel

def create_ca_pool_proto(poolSettings):
    if(not(moose.exists("/library"))):
        moose.Neutral('/library')
    pool = moose.CaConc('/library/' + poolSettings.name)
    pool.CaBasal=poolSettings.CaBasal
    pool.ceiling=1
    pool.floor = 0
    pool.thick = poolSettings.CaThick
    pool.tau = poolSettings.CaTau
    return pool

def add_calcium(cellname, poolSettings, chan_name, chan_type, calname):
    FARADAY = 96485.33289
    pool_proto = create_ca_pool_proto(poolSettings)
    for comp in moose.wildcardFind("%s/#[TYPE=Compartment]"%(cellname)):
        pool = moose.copy(pool_proto, comp, poolSettings.name,1)
        print(pool)
        pool.length = comp.length
        pool.diameter = comp.diameter
        SA = np.pi*pool.length*pool.diameter
        vol = SA*pool.thick
        pool.B = 1/(FARADAY*vol*2)/poolSettings.BufCapacity
        chan = moose.element(comp.path+"/"+chan_name)
        if chan_type=="ca_permeable":
            m = moose.connect(chan,"IKOut", )

# channelSpecs must be a dictionary with keys:
#    settings - ChannelSettings instance
#    gbar - conductance of channel for this component
#    Ek -
def addChannelToComps(channelSpecs, comps, number=1):
    container = comps[0].parent.path
    protoChannel = create_channel_proto(channelSpecs["settings"])
    print(container)
    print(protoChannel)
    channel = moose.copy(protoChannel, container, channelSpecs["settings"].name+'_{}'.format(comps.name), number)
    channel.Gbar = channelSpecs["gbar"] * np.pi*comps.length*comps.diameter # is the list unnecessary?
    channel.Ek = channelSpecs["Ek"]
    print(moose.showfield(moose.element(channel)))
    moose.connect(channel, "channel", comps, "channel", "OneToOne")

def create_1comp_neuron(path, number=1):
    """Create single-compartmental neuron with Na+ and K+ channels.

    Parameters
    ----------
    path : str
        path of the compartment to be created

    number : int
        number of compartments to be created. If n is greater than 1,
        we create a vec with that size, each having the same property.

    Returns
    -------
    comp : moose.Compartment
        a compartment vec with `number` elements.

    """
    comps = moose.vec(path=path, n=number, dtype='Compartment')
    diameter = 30e-6
    length = 50e-6
    comps.diameter = diameter
    comps.length = length
    sarea = np.pi * diameter * length
    xarea = np.pi * diameter * diameter / 4.0
    Em = EREST_ACT + 10.613e-3
    comps.Em = Em
    comps.initVm = EREST_ACT
    #: CM = 1 uF/cm^2
    comps.Cm = 1e-6 * sarea * 1e4
    #: RM = 0.3 mS/cm^2
    comps.Rm = 1 / (0.3e-3 * sarea * 1e4)
    channels = [{
            "settings" : Na_settings,
            "gbar" : 1200,
            "Ek" : 115e-3 + EREST_ACT
        }, {
            "settings" : K_settings,
            "gbar" : 360,
            "Ek" : -12e-3 + EREST_ACT
        #}, {
        #    "settings" : KDr_settings,
        #    "gbar" : 0,
        #    "Ek" : -12e-3 + EREST_ACT #??
        }, {
            "settings" : SK_settings,
            "gbar" : 2,
            "Ek" : -12e-3 + EREST_ACT #??
        }, {
            "settings" : CaL_settings,
            "gbar" : 15,
            "Ek" : -12e-3 + EREST_ACT #??
        }]
    for channel in channels:
        addChannelToComps(channel, comps)
    return comps

def create_compartment(name,length,diameter,surfaceArea,xArea,membraneResistivity,axialResistivity,capacitivity,Em,initVm): # Em, initVm?
	newC = moose.Compartment(name)
	newC.length = length
	newC.diameter = diameter
	newC.Em = Em
	newC.initVm = initVm
	newC.Rm = membraneResistivity * 1.0 / surfaceArea
	newC.Cm = capacitivity * surfaceArea
	newC.Ra = axialResistivity * length * 1.0 / xArea
	return newC

def create_spherical_compartment(name,length,diameter,membraneResistivity,axialResistivity,capacitivity,Em,initVm):
	surfaceArea = np.pi*diameter**2
	xArea = surfaceArea * 1.0 /4
	return create_compartment(name,length,diameter,surfaceArea,xArea,membraneResistivity,axialResistivity,capacitivity,Em,initVm)

def addChannels(comp):
        container = comp.parent.path
        nachan = moose.copy(create_na_proto(), container, 'na_{}'.format(comp.name), 1)
        nachan.Gbar = [120e-3 * np.pi*comp.diameter**2 * 1e4] * len(nachan)
        nachan.Ek = 115e-3 + EREST_ACT
        moose.connect(nachan, 'channel', comp, 'channel', 'OneToOne')
        kchan = moose.copy(create_k_proto(), container, 'k_{}'.format(comp.name), 1)
        kchan.Gbar = 36e-3 * np.pi*comp.diameter**2 * 1e4
        kchan.Ek = -12e-3 + EREST_ACT
        moose.connect(kchan, 'channel', comp, 'channel', 'OneToOne')

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
        x = np.fromiter((i*table.dt for i in range(table.size)), float, table.size)
        y = table.vector
        plt.plot(x, y)
    plt.show()

def current_step_test(simtime, simdt, plotdt):
    """
    Create a single compartment and set it up for applying a step
    current injection.

    We use a PulseGen object to generate a 40 ms wide 1 nA current
    pulse that starts 20 ms after start of simulation.

    """
    model = moose.Neutral('/model')
    comp = create_1comp_neuron('/model/neuron')
    #comp = create_spherical_compartment("/model/neuron",50e-6,30e-6,0.3e-3,0.3e-3,1e-6,EREST_ACT + 10.613e-3,EREST_ACT)
    #addChannels(comp)
    #stim = moose.PulseGen('/model/stimulus')
    #stim.delay[0] = 20e-3
    #stim.level[0] = 1e-9
    #stim.width[0] = 40e-3
    #stim.delay[1] = 1e9
    #moose.connect(stim, 'output', comp, 'injectMsg')

    stim = create_pulse("/model/stimulus", 20e-3, 40e-3, 1e-9, comp)

    data = moose.Neutral('/data')
    current_tab = create_table("/data/current", stim, "getOutputValue")
    vm_tab = create_table("/data/Vm", comp, "getVm")
    for i in range(10):
        moose.setClock(i, simdt)
    moose.setClock(8, plotdt)
    moose.reinit()
    moose.start(simtime)
    ts = np.linspace(0, simtime, len(vm_tab.vector))
    vm_tab.vector = vm_tab.vector * 1e3
    current_tab.vector = current_tab.vector * 1e9
    return ts, current_tab, vm_tab

def main():
    """
This demo shows how to set the parameters for a Hodgkin-Huxley type ion channel.

Hodgkin-Huxley type ion channels are composed of one or more gates
that allow ions to cross the membrane. The gates transition between
open and closed states and this, taken over a large population of
ion channels over a patch of membrane has first order kinetics, where
the rate of change of fraction of open gates (n) is given by::

    dn/dt = alpha(Vm) * (1 - n) - beta(Vm) * n

where alpha and beta are rate parameters for gate opening and
closing respectively that depend on the membrane potential.
The final channel conductance is computed as::

    Gbar * m^x * h^y ...

where m, n are the fraction of open gates of different types and x,
y are the number of such gates in each channel. We can define the
channel by specifying the alpha and beta parameters as functions of
membrane potential and the exponents for each gate.
The number gates is commonly one or two.

Gate opening/closing rates have the form::

    y(x) = (A + B * x) / (C + exp((x + D) / F))

where x is membrane voltage and y is the rate parameter for gate
closing or opening.
    """
    bloop = moose.loadModel("layer2.p", "bloop")
    add_calcium(bloop.name, Ca_pool_settings)
    moose.le(bloop)
    simtime = 0.1
    simdt = 0.25e-5
    plotdt = 0.25e-3
    ts, current, vm = current_step_test(simtime, simdt, plotdt)
    plot_tables([vm, current])

if __name__ == '__main__':
    main()

#
# ionchannel.py ends here
