import numpy as np
import parameters as para



class Weight:
    def __init__(self, para):
        # parameters.py is a file containing all constant parameters of the aircraft
        self.para = para
        # take_off weight that is assumed initially        


    #########################################################################
    "CLASS II WEIGHT ESTIMATION"
    #########################################################################

    def wing_weight(self):   #Span, sweep, ultimate load factor, thickness over chord, cord length at root, wing loading, gross weight
        k_w = 4.9*10**(-3)
        b_s = self.para.b / np.cos(self.para.lambda_mid)
        b_ref = 1.905
        t_r = self.para.t_c * self.para.cwr

        W_w = k_w * b_s**0.75 * (1 + (b_ref/self.para.b_s)**0.5) * self.para.n_ult**0.55 * ((b_s/t_r)/self.para.W_loading)**0.3 * self.para.W_G
        
        #ADD 30% IF BRACED WING USED, 10% IF STRUT USED?
        
        return W_w
    
    def tail_weight(self): #ultimate load factor, tail surface area
        W_t = 0.64 * (self.para.n_ult * self.para.s_tail**2)**0.75

        #IF TAILPLANE AREA UNKNOWN, WEIGHT ASSUMED TO BE 3.5-4% OF EMPTY WEIGHT
        return W_t

    def gear_weight(self):
        A1 = input("Give A for main gear (Table 8-6 Torenbeek)", )
        B1 = input("Give B for main gear (Table 8-6 Torenbeek)", )
        C1 = input("Give C for main gear (Table 8-6 Torenbeek)", )
        D1 = input("Give D for main gear (Table 8-6 Torenbeek)", )
        A2 = input("Give A for nose gear (Table 8-6 Torenbeek)", )
        B2 = input("Give B for nose gear (Table 8-6 Torenbeek)", )
        C2 = input("Give C for nose gear (Table 8-6 Torenbeek)", )
        D2 = input("Give D for nose gear (Table 8-6 Torenbeek)", )

        #Main gear:
        W_uc1 = 1.08 * (A1 + B1 * self.para.W_TO**0.75 + C1 * self.para.W_TO + D1 * self.para.W_TO**1.5)

        #Nose gear:
        W_uc2 = 1.08 * (A2 + B2 * self.para.W_TO**0.75 + C2 * self.para.W_TO + D2 * self.para.W_TO**1.5)

        return W_uc1 + W_uc2

    def nacelle_weight(self): #Take off power in hp
        W_n = 1.134 * self.para.P_TO**0.5
        return W_n


    def equipment_weight(self):
        return 0.08 * self.para.W_TO #MORE DETAILED ESTIMATION CAN BE MADE BUT NOT NECESSARY FOR TRADE-OFF

    def fuselage_weight(self):
        k_wf = 0.23
        W_f = self.para.k_wf*(self.para.V_D*(self.para.l_t/(self.para.b_f+self.para.h_f)))**0.5*self.para.S_G**1.2

        # Add 7% if the main landing gear is attached to the fuselage, but 4% can be subtracted from the fuselage weight if there is no attachment structure for the main landing gear
        # Add 10% for freighter aircraft
        # For booms l_t is defined as distance between local wing chord and horizontal tailplane

        return W_f

    def control_surface_weight(self):
        k_sc = 0.44*0.768
        W_sc = k_sc*self.para.W_TO**(2/3)

        # Add 20% for LE flap or slat
        # Add 15% for lift dumper controls

        return W_sc

    def propulsion_weight(self):
        k_pg = 1.16 # tractor single propeller aircraft
        W_pg = k_pg*self.para.N_e*(self.para.W_e+0.109*self.para.P_TO)

        # If number of cylinder and volume of cylinder are known use figure 4-12 Torenbeek

        return W_pg


if __name__ == "__main__":
    weight = Weight(para)







