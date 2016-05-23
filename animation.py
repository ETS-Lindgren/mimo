from math import cos as cosine
from math import sin as sine
from math import pi as PI
from math import log10

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class UMA(object):

    def __init__(self):
        Power_array = [0, -2.2, -1.7, -5.2, -9.1, -12.5]
        Power_array = [10**(each/10.0) for each in Power_array]
        N_path = len(Power_array)  ## number of clusters
        AoA_array = [65.7, 45.6, 143.2, 32.5, -91.1, -19.2]
        AoD_array = [82,   80.5, 79.6,  98.6, 102.1, 107.1]
        AS_AoD = 5.0
        AS_AoA = 35.0
        self.Frequency = 7.51e6  ## in Hz
        self.V_MA = 100000.0 / 60 / 60    ## 100km/hr, mobile speed in m/s
        self.V_direction = 120.0
        XPR = 9           ## in dB

        subpath10_1 = [-75.4274, -53.1816, -40.1824, -30.9538, -23.7899, -17.9492, -13.0045, -8.7224, -4.9447, -1.4649]
        subpath10_2 = [-1*each for each in subpath10_1[::-1]]
        subpath20 = subpath10_1 + subpath10_2
        subpath20_AoA = list(subpath20)        
        subpath20_AoD = [each / AS_AoA * AS_AoD for each in subpath20] ## normalizaing subpath, so can be extened later for either AoA or AoD

        self.AoD_all = []
        for each_AoD in AoD_array:
            for each_subpath in subpath20_AoD:
                self.AoD_all.append(each_AoD + each_subpath)

        self.AoA_all = []
        for each_AoA in AoA_array:
            for each_subpath in subpath20_AoA:
                self.AoA_all.append(each_AoA + each_subpath)

        ## BS antenna model
        V_power_BS = [ 0.5 for each in range(120)]
        H_power_BS = [ pow(cosine(each/180.0*PI), 2) for each in self.AoD_all]

        self.V_power_BS = V_power_BS
        self.H_power_BS = H_power_BS
        
        ## XPR
        XPR = 10**(XPR/10.0)
        co_pol = XPR / (XPR + 1.0)
        cr_pol = 1 / (XPR + 1.0)

        ## BS power after XPR
        power_from_BS_V = [ V_power_BS[i]*co_pol + H_power_BS[i]*cr_pol for i in range(120) ]
        power_from_BS_H = [ H_power_BS[i]*co_pol + V_power_BS[i]*cr_pol for i in range(120) ]

        self.power_from_BS_V = power_from_BS_V
        self.power_from_BS_H = power_from_BS_H

        ## calculate the V and H power on each subpath from BS to the UE
        power_array_all = []
        for each in Power_array:
            for each_subpath in range(20):
                power_array_all.append(each/20.0)
        self.power_from_BS_to_UE_V = []
        self.power_from_BS_to_UE_H = []
        
        for each in range(120):
            self.power_from_BS_to_UE_V.append( power_from_BS_V[each]*power_array_all[each] )
            self.power_from_BS_to_UE_H.append( power_from_BS_H[each]*power_array_all[each] )

        self.cross_pol = 10*log10(sum(self.power_from_BS_to_UE_V) / sum(self.power_from_BS_to_UE_H))
        
    def ani(self, cross_len = 2.0, nos = 80):
        """cross length (corss_len) is assigned as two lambdas by default.  nos means number of sections."""
        self.ani_x = {}
        self.ani_y = {}
        spatial_max = 0
        delta = float(cross_len) / nos
        
        for n in range(nos + 1):
            step = n*delta
            data = []
            for k, each_AoA in enumerate(self.AoA_all):
                step2 = step*sine(each_AoA/180.0*PI)*2*PI
                phase = cosine(step2) + 1j*sine(step2)
                data.append(self.power_from_BS_to_UE_V[k]*phase)
                self.ani_x[n] = np.array([each.real for each in data])
                self.ani_y[n] = np.array([each.imag for each in data])
            
    def update(self, line, count):
        x_data = self.ani_x[count]
        y_data = self.ani_y[count]
#        print(x_data)
#        print(y_data)
        line.set_data(x_data, y_data)
        return line,
    
                
if __name__ == '__main__':
    import sys
    if len(sys.argv) == 2:
        time_int = sys.argv[1]
    else:
        time_int = 50
    model = UMA()
    nos = 200
    model.ani(5,200)
    fig = plt.figure(figsize = (7.5,7.2))
    line, = plt.plot([],[], '.')
    
    def update_data(count):
        return model.update(line, count)        

    plt.xlim(-0.025, 0.025)
    plt.ylim(-0.025, 0.025)
    plt.xlabel('real')
    plt.ylabel('imag')    
    plt.title('Spatial Correlation - UMA')


    from time import sleep
    line_ani = animation.FuncAnimation(fig, update_data, np.arange(0,nos-1), interval=time_int, blit=True)

    plt.show()        
    
    print("\n")
    print("Modeling of UMA\n")
    print("\n")
