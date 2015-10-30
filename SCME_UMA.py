from math import cos as cosine
from math import sin as sine
from math import pi as PI
from math import log10

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
        
    def spatial_corr(self, cross_len = 2.0, nos = 80):
        """Spatial Correlatoin, cross length (corss_len) is assigned as two lambdas by default.  nos means number of sections."""
        self.spatial_correlation = []
        spatial_max = 0
        step = float(cross_len) / nos
        for each in range(nos + 1):
            temp = 0
            for k, each_AoA in enumerate(self.AoA_all):
                each_step = each*step*sine(each_AoA/180.0*PI)*2*PI
                phase = cosine(each_step) + 1j*sine(each_step)
                temp = temp + self.power_from_BS_to_UE_V[k]*phase
            self.spatial_correlation.append(abs(temp))
        self.spatial_correlation = [each / self.spatial_correlation[0] for each in self.spatial_correlation]

        step_tenth_lambda = 0.1
        step_half_lambda = 0.5
        step_full_lambda = 1.0
        self.power_from_BS_to_UE_tenth_lambda = []
        self.power_from_BS_to_UE_half_lambda = []
        self.power_from_BS_to_UE_full_lambda = []
        for k, each_AoA in enumerate(self.AoA_all):
            step = step_tenth_lambda*sine(each_AoA/180.0*PI)*2*PI
            phase = cosine(step) + 1j*sine(step)            
            self.power_from_BS_to_UE_tenth_lambda.append(self.power_from_BS_to_UE_V[k]*phase)
            
            step = step_half_lambda*sine(each_AoA/180.0*PI)*2*PI
            phase = cosine(step) + 1j*sine(step)
            self.power_from_BS_to_UE_half_lambda.append(self.power_from_BS_to_UE_V[k]*phase)
            
            step = step_full_lambda*sine(each_AoA/180.0*PI)*2*PI
            phase = cosine(step) + 1j*sine(step)
            self.power_from_BS_to_UE_full_lambda.append(self.power_from_BS_to_UE_V[k]*phase)
            
if __name__ == '__main__':
    model = UMA()
    model.spatial_corr()
    
    with open("data_AoD.txt", "w") as f, open("data_BS_V", "w") as fv, open("data_BS_H", "w") as fh:
        for k, each in enumerate(model.AoD_all):
            f.write( str(each) + "\t" + str(model.power_from_BS_V[k]) + "\n")
            fv.write( str(each) + "\t" + str(model.power_from_BS_V[k]) + "\n")
            fh.write( str(each) + "\t" + str(model.power_from_BS_H[k]) + "\n")
            
    with open("data_AoD_polar.txt", "w") as f:
        for k, each in enumerate(model.AoD_all):        
            f.write( str(each/180.0*PI) + "\t" + str(10*log10(model.power_from_BS_V[k])) + "\n")
            
    with open("data_AoA.txt", "w") as g, open("data_AoA_tenth.txt", "w") as gt, open("data_AoA_half.txt", "w") as gh, open("data_AoA_full.txt", "w") as gf:
        for k, each in enumerate(model.AoA_all):
            g.write( str(each) + "\t" + str(model.power_from_BS_to_UE_V[k]) + "\n")
            tenth_complex = model.power_from_BS_to_UE_tenth_lambda[k]
            half_complex = model.power_from_BS_to_UE_half_lambda[k]
            full_complex = model.power_from_BS_to_UE_full_lambda[k]
            gt.write( str(each) + "\t" + str(tenth_complex.real) + "\t" + str(tenth_complex.imag) + "\n")            
            gh.write( str(each) + "\t" + str(half_complex.real) + "\t" + str(half_complex.imag) + "\n")
            gf.write( str(each) + "\t" + str(full_complex.real) + "\t" + str(full_complex.imag) + "\n")


    with open("data_AoA_polar.txt", "w") as g:
        for k, each in enumerate(model.AoA_all):            
            g.write( str(each/180.0*PI) + "\t" + str(10*log10(model.power_from_BS_to_UE_V[k])) + "\n")

    
    with open("data_spatial_corr_UMA_V.txt", "w") as h:
        for k, each in enumerate(model.spatial_correlation):
            h.write( str(k*0.025) + "\t" + str(each) + "\n")

    print("\n")
    print("Modeling of UMA\n")
    print("Spatial Correlation is output for Vertical Polarization: \ndata_spatial_corr_UMA_V.txt")
    print("\n")
