import numpy as np
import matplotlib.pyplot as plt
import moose
import collections

EREST_ACT = -70e-3 #: Resting membrane potential

AbParams = collections.namedtuple("AbParams",[
    "aRate",    # 'A_A'  
    "aB",       # 'A_B'
    "aC",       # 'A_C'
    "aVHalf",   # 'A_D'
    "aVSlope",  # 'A_F'
    "bRate",    # 'B_A'
    "bB",       # 'B_B'
    "bC",       # 'B_C'
    "bVHalf",   # 'B_D'
    "bVSlope"   # 'B_F'
    ])

#: The parameters for defining m as a function of Vm
Na_m_params = AbParams(
      aRate = 1e5 * (25e-3 + EREST_ACT), 
         aB = -1e5,
         aC = -1.0,
     aVHalf = -25e-3 - EREST_ACT,
    aVSlope = -10e-3,
      bRate = 4e3,
         bB = 0.0,
         bC = 0.0,
     bVHalf = 0.0 - EREST_ACT,
    bVSlope = 18e-3
    )
#: Parameters for defining h gate of Na+ channel
Na_h_params = AbParams(
      aRate = 70.0,
         aB = 0.0,
         aC = 0.0,
     aVHalf = 0.0 - EREST_ACT,
    aVSlope = 0.02,
      bRate = 1000.0,
         bB = 0.0,
         bC = 1.0,
     bVHalf = -30e-3 - EREST_ACT,
    bVSlope = -0.01
    )

#: K+ channel in Hodgkin-Huxley model has only one gate, n and these
#are the parameters for the same
K_n_params = AbParams(
      aRate = 1e4 * (10e-3 + EREST_ACT),
         aB = -1e4,
         aC = -1.0,
     aVHalf = -10e-3 - EREST_ACT,
    aVSlope = -10e-3,
      bRate = 0.125e3,
         bB = 0.0,
         bC = 0.0,
     bVHalf = 0.0 - EREST_ACT,
    bVSlope = 80e-3
    )

KDr_X_params = AbParams(
      aRate = 28.2,
         aB = 0,
         aC = 0.0,
     aVHalf = 0,
    aVSlope = -12.5e-3,
      bRate = 6.78,
         bB = 0.0,
         bC = 0.0,
     bVHalf = 0.0,
    bVSlope = 33.5e-3
    )

CaL_X_params = AbParams(
      aRate = -880,
         aB = -220e3,
         aC = -1.0,
     aVHalf = 4.0003e-3,
    aVSlope = -7.5e-3,
      bRate = -284,
         bB = 71e3,
         bC = -1.0,
     bVHalf = -4.0003e-3,
    bVSlope = 5e-3
    )

CaDepParams = collections.namedtuple("CaDepParams", [
    "Kd",
    "power",
    "tau"
    ])

SK_Z_params = CaDepParams(
       Kd = 0.57e-3,
    power = 5.2,
      tau = 4.9e-3
    )

PoolSettings = collections.namedtuple("PoolSettings", [
    "CaBasal",
    "CaThick",
    "CaTau",
    "BufCapacity",
    "name"
    ])

Ca_pool_settings = PoolSettings(
        CaBasal = 50e-6,    
        CaThick = 1e-6,    
          CaTau = 20e-3,    
    BufCapacity = 20,    
           name = "CaPool"
    )

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


Na_settings = ChannelSettings(
       xPower = 3,
       yPower = 1,
       zPower = 0,
         eRev = 0.06,
         name = "na",
       xParam = Na_m_params,
       yParam = Na_h_params,
       zParam = None,
    chan_type = ""
    )

K_settings = ChannelSettings(
       xPower = 4,
       yPower = 0,
       zPower = 0,
         eRev = -0.1,
         name = "k",
       xParam = K_n_params,
       yParam = None,
       zParam = None,
    chan_type = ""
    )

KDr_settings = ChannelSettings(
       xPower = 2,
       yPower = 0,
       zPower = 0,
         eRev = 90e-3,
         name = "kdr",
       xParam = KDr_X_params,
       yParam = None,
       zParam = None,
    chan_type = ""
    )

SK_settings = ChannelSettings(
       xPower = 0,
       yPower = 0,
       zPower = 1,
         eRev = -87e-3,
         name = "skca",
       xParam = None,
       yParam = None,
       zParam = SK_Z_params,
    chan_type = "ca_dependent"
    )

CaL_settings = ChannelSettings(
       xPower = 1,
       yPower = 0,
       zPower = 0,
         eRev = 130e-3,
         name = "cal",
       xParam = CaL_X_params,
       yParam = None,
       zParam = None,
    chan_type = "ca_permeable"
    )

# Ca_pool_settings = PoolSettings(50e-6, 1e-6, 20e-3, 20, "CaPool")

# Na_settings = ChannelSettings(3,1,0,0.06,"na",Na_m_params,Na_h_params,None,"")

# K_settings = ChannelSettings(4,0,0,-0.1,"k",K_n_params,None,None,"")

# KDr_settings = ChannelSettings(2,0,0,90e-3,"kdr",KDr_X_params,None,None,"")

# SK_settings = ChannelSettings(0,0,1,-87e-3,"skca",None,None,SK_Z_params,"ca_dependent")

# CaL_settings = ChannelSettings(1,0,0,130e-3,"cal",CaL_X_params,None,None,"ca_permeable")



#: We define the rate parameters, which are functions of Vm as
#: interpolation tables looked up by membrane potential.
#: Minimum x-value for the interpolation table
VMIN = -30e-3 + EREST_ACT
#: Maximum x-value for the interpolation table
VMAX = 120e-3 + EREST_ACT
#: Number of divisions in the interpolation table
VDIVS = 3000

def create_channel_proto(channelSettings,
        vDivs  = 3000,
        vMin   = -30e-3 + EREST_ACT,
        vMax   = 120e-3 + EREST_ACT,
        caDivs = 10000,
        caMin  = 0,
        caMax  = 1):
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
            m = moose.connect(chan,"IKOut", pool, "current")
        elif chan_type=="ca_dependent":
            m = moose.connect(pool,"concOut", pool, "concen")
        else:
            print("unknown calcium concentration type")

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
    if(channelSpecs["settings"].chan_type!=""):
        print(channelSpecs["settings"].chan_type)
        add_calcium(channel.name, Ca_pool_settings, channel.name, channelSpecs["settings"].chan_type, 'CaPool')
    return channel

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
        #    "Ek" : -12e-3 + EREST_ACT
        }, {
            "settings" : SK_settings,
            "gbar" : 0,
            "Ek" : -12e-3 + EREST_ACT
        }, {
            "settings" : CaL_settings,
            "gbar" : 0,
            "Ek" : -12e-3 + EREST_ACT
        }]
    for channel in channels:
        addChannelToComps(channel, comps)
    return comps

# NOTE: This function is currently unused.
#def create_compartment(name,length,diameter,surfaceArea,xArea,membraneResistivity,axialResistivity,capacitivity,Em,initVm): # Em, initVm?
#	newC = moose.Compartment(name)
#	newC.length = length
#	newC.diameter = diameter
#	newC.Em = Em
#	newC.initVm = initVm
#	newC.Rm = membraneResistivity * 1.0 / surfaceArea
#	newC.Cm = capacitivity * surfaceArea
#	newC.Ra = axialResistivity * length * 1.0 / xArea
#	return newC

# NOTE: This function is currently unused.
#def create_spherical_compartment(name,length,diameter,membraneResistivity,axialResistivity,capacitivity,Em,initVm):
#	surfaceArea = np.pi*diameter**2
#	xArea = surfaceArea * 1.0 /4
#	return create_compartment(name,length,diameter,surfaceArea,xArea,membraneResistivity,axialResistivity,capacitivity,Em,initVm)

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

def main():
    simtime = 0.1
    simdt = 0.25e-5
    plotdt = 0.25e-3
    model = moose.Neutral('/model')
    comp = create_1comp_neuron('/model/neuron')
    #comp = create_spherical_compartment("/model/neuron",50e-6,30e-6,0.3e-3,0.3e-3,1e-6,EREST_ACT + 10.613e-3,EREST_ACT)
    stim = create_pulse("/model/stimulus", 20e-3, 40e-3, 1e-9, comp)

    data = moose.Neutral('/data')
    current_tab = create_table("/data/current", stim, "getOutputValue")
    vm_tab = create_table("/data/Vm", comp, "getVm")
    for i in range(10):
        moose.setClock(i, simdt)
    moose.setClock(8, plotdt)
    moose.reinit()
    moose.start(simtime)
    #ts = np.linspace(0, simtime, len(vm_tab.vector))
    plot_tables([vm_tab, current_tab])

if __name__ == '__main__':
    main()
