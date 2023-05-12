import numpy as np
import parameters as para



class Weight:
    def __init__(self, para):
        # parameters.py is a file containing all constant parameters of the aircraft
        self.para = para
        self.W_TO = 750
        # take_off weight that is assumed initially


    #########################################################################
    "CLASS II WEIGHT ESTIMATION"
    #########################################################################

    def wing_weight(self):   #Span, sweep, ultimate load factor, thickness over chord, cord length at root, wing loading, gross weight
        k_w = 4.9E-3
        b_s = self.para.b / np.cos(self.para.lambda_mid)
        b_ref = 1.905
        t_r = self.para.t_c * self.para.cwr

        self.W_w = k_w * b_s**0.75 * (1 + (b_ref/b_s)**0.5) * self.para.n_ult**0.55 * ((b_s/t_r)/self.para.W_loading)**0.3 * self.para.W_G
        self.W_w = k_w * b_s**0.75 * (1 + (b_ref/b_s)**0.5) * self.para.n_ult**0.55 * ((b_s/t_r)/self.para.W_loading)**0.3 * self.para.W_G
        
        #ADD 30% IF BRACED WING USED, 10% IF STRUT USED?

    def tail_weight(self): #ultimate load factor, tail surface area
        self.W_t = 0.64 * (self.para.n_ult * self.para.s_tail**2)**0.75

        #IF TAILPLANE AREA UNKNOWN, WEIGHT ASSUMED TO BE 3.5-4% OF EMPTY WEIGHT

    def gear_weight(self):
        A1 = float(input("Give A for main gear (Table 8-6 Torenbeek)", ))
        B1 = float(input("Give B for main gear (Table 8-6 Torenbeek)", ))
        C1 = float(input("Give C for main gear (Table 8-6 Torenbeek)", ))
        D1 = float(input("Give D for main gear (Table 8-6 Torenbeek)", ))
        A2 = float(input("Give A for nose gear (Table 8-6 Torenbeek)", ))
        B2 = float(input("Give B for nose gear (Table 8-6 Torenbeek)", ))
        C2 = float(input("Give C for nose gear (Table 8-6 Torenbeek)", ))
        D2 = float(input("Give D for nose gear (Table 8-6 Torenbeek)", ))

        #Main gear:
        W_uc1 = 1.08 * (A1 + B1 * self.para.W_TO**0.75 + C1 * self.para.W_TO + D1 * self.para.W_TO**1.5)

        #Nose gear:
        W_uc2 = 1.08 * (A2 + B2 * self.para.W_TO**0.75 + C2 * self.para.W_TO + D2 * self.para.W_TO**1.5)
        self.W_uc = W_uc1 + W_uc2

    def nacelle_weight(self): #Take off power in hp
        self.W_n = 1.134 * self.para.P_TO**0.5


    def equipment_weight(self):
        self.W_eq = 0.08 * self.para.W_TO
        #MORE DETAILED ESTIMATION CAN BE MADE BUT NOT NECESSARY FOR TRADE-OFF

    def fuselage_weight(self):
        k_wf = 0.23
        self.W_f = k_wf*(self.para.V_D*(self.para.l_t/(self.para.b_f+self.para.h_f)))**0.5*self.para.S_G**1.2
        self.W_f = k_wf*(self.para.V_D*(self.para.l_t/(self.para.b_f+self.para.h_f)))**0.5*self.para.S_G**1.2

        # Add 7% if the main landing gear is attached to the fuselage, but 4% can be subtracted from the fuselage weight if there is no attachment structure for the main landing gear
        # Add 10% for freighter aircraft
        # For booms l_t is defined as distance between local wing chord and horizontal tailplane

    def control_surface_weight(self):
        k_sc = 0.44*0.768
        self.W_sc = k_sc*self.para.W_TO**(2/3)

        # Add 20% for LE flap or slat
        # Add 15% for lift dumper controls


    def propulsion_weight(self):
        k_pg = 1.16 # tractor single propeller aircraft
        self.W_pg = k_pg*self.para.N_e*(self.para.W_e+0.109*self.para.P_TO)

        # If number of cylinder and volume of cylinder are known use figure 4-12 Torenbeek

    def weight_empty(self):
        self.W_OEW = self.W_pg + self.W_sc + self.W_f + self.W_eq + self.W_n + self.W_uc + self.W_t + self.W_w

        print(f"W_OEW:{self.W_OEW}")
        return self.W_OEW


    def cg_calc(self):
        self.wing_cg = 0
        if self.para.sweep_angle == 0:
            self.wing_cg = 0.4 * self.cwr + self.l_LE #40% of root chord plus Leading Edge location
        else:
            self.wing_cg = 0 # to be done later, depends on spar locations (table 8-15 Torenbeek)
        
        self.fus_cg = 0.335 * self.para.l_f #32-35% of fuselage length
        self.tail_cg = 0.42 * self.para.cwr + self.para.l_LE #42% of root chord plus Leading Edge location

        self.engine_cg = 0 # to be done later
        # c.g. or rotax 912is is at 32.7 cm from the front of the engine

        self.landing_gear_cg = 0 # to be done later
        # can be at airplane c.g. -> iteration needed, or use location main and nose landing gear

        self.cg = (self.wing_cg * self.W_w + self.fus_cg * self.W_f + self.tail_cg * self.W_t + self.engine_cg * self.W_pg + self.landing_gear_cg * self.W_uc) / (self.W_w + self.W_f + self.W_t + self.W_pg + self.W_uc)

if __name__ == "__main__":
    weight = Weight(para)
    weight.wing_weight()
    weight.tail_weight()
    weight.gear_weight()
    weight.nacelle_weight()
    weight.equipment_weight()
    weight.fuselage_weight()
    weight.control_surface_weight()
    weight.propulsion_weight()
    weight.weight_empty()
    weight.cg_calc()