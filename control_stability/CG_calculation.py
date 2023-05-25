import numpy as np
import matplotlib.pyplot as plt



#Code that will be used for generation of loading diagram

#Parameters

reference_aircraft = {
    "Name" : "ATR72-600_reference", #Name of aircraft
    "OEW": 13178.4,                   #kg
    "MTOW": 22800.00,               #kg
    "Fuel_max": 5000,               #kg
    "max_payload": 7500,            #kg
    "pax": 72,                      #Number of passengers
    "pax_weight": 95,               #Weight of 1 passengers [kg]
    "front_volume": 4.6,            #Volume of front cargo bay [m^3]
    "aft_volume": 4.8,              #Volume of aft cargo bay [m^3]
    "MAC": 2.45,                    #Mean Aerodynamic Chord of aircraft
    "xlemac": 11.43,                #x position of the leading edge of the mean aerodynamic chord
    "xcg_OEW": 1.027711073,         #CG location OEW /MAC
    "seat_pitch": 29,               #Seat pitch [inch]
    "xcg_front_seat": 6.25,         #x position first seat from nose [m]
    "xcg_front_CB": 4.358164881,    #xcg front cargo bay [m]
    "xcg_aft_CB": 21.47884442,      #xcg aft cargo bay [m]
    "xcg_fuel": 12.655              #xcg of the fuel tanks [m]
}

modified_aircraft = {
    "Name" : "ATR72-600_modified", #Name of aircraft
    "OEW"  : 13837.13,             #kg
    "MTOW" : 22800.00,             #kg
    "Fuel_max" : 5000,             #kg #This needs to be reconsidered, less fuel
    "max_payload" : 7500 - 760,          #kg
    "pax" : 64,                    #Number of passengers
    "pax_weight" : 95,             #Weight of 1 passengers [kg]
    "front_volume" : 4.6,          #Volume of front cargo bay [m^3]
    "aft_volume" : 4.8,            #Volume of aft cargo bay [m^3]
    "MAC" : 2.45,                  #Mean Aerodynamic Chord of aircraft
    "xlemac" : 11.43,              #x position of the leading edge of the mean aerodynamic chord
    "xcg_OEW" : 1.300000062,       #CG location OEW /MAC
    "seat_pitch" : 29,             #Seat pitch [inch]
    "xcg_front_seat" : 6.25,       #x position first seat from nose [m]
    "xcg_front_CB" : 4.358164881,  #xcg front cargo bay [m]
    "xcg_aft_CB" : 21.47884442,    #xcg aft cargo bay [m]
    "xcg_fuel" : 12.655            #xcg of the fuel tanks [m]
}


def loading_diagram(aircraft_dict,plotting=False):
    #Unpack variables
    OEW = aircraft_dict["OEW"]
    ac_name = aircraft_dict["Name"]
    MTOW = aircraft_dict["MTOW"]
    Fuel_max = aircraft_dict["Fuel_max"]
    max_payload = aircraft_dict["max_payload"]
    pax = aircraft_dict["pax"]
    pax_weight = aircraft_dict["pax_weight"] * 2
    front_volume = aircraft_dict["front_volume"]
    aft_volume = aircraft_dict["aft_volume"]
    payload_weight = max_payload - aircraft_dict['pax_weight'] * pax
    front_cargo_weight = front_volume/(front_volume + aft_volume) * payload_weight
    aft_cargo_weight = aft_volume/(front_volume + aft_volume) * payload_weight
    MAC = aircraft_dict["MAC"]
    xlemac = aircraft_dict["xlemac"]
    xcg_OEW = aircraft_dict["xcg_OEW"]
    seat_pitch = aircraft_dict["seat_pitch"] * 2.54 * 1/100
    xcg_front_seat = aircraft_dict["xcg_front_seat"]
    xcg_front_CB = aircraft_dict["xcg_front_CB"]
    xcg_aft_CB = aircraft_dict["xcg_aft_CB"]
    xcg_fuel = aircraft_dict["xcg_fuel"]

    # Functions
    def convert_mac(x):
        return x - xlemac

    def cargoload(current_mass, current_cg):
        # Calculate cargoloading CG points

        front_load_cg = [current_cg]  # List of CG points for front loading
        aft_load_cg = [current_cg]  # List of CG points for aft loading

        aft_load_mass = [current_mass]  # List of masses for aft loading
        front_load_mass = [current_mass]  # List of masses for front loading

        front_load_cg.append((OEW_moment + front_moment) / (OEW + front_cargo_weight))  # CG point for front loading
        aft_load_cg.append((OEW_moment + aft_moment) / (OEW + aft_cargo_weight))  # CG point for aft loading

        front_load_mass.append(OEW + front_cargo_weight)  # Mass for front loading
        aft_load_mass.append(OEW + aft_cargo_weight)  # Mass for aft loading

        final_cg = (OEW_moment + aft_moment + front_moment) / (
                OEW + front_cargo_weight + aft_cargo_weight)  # Final CG point
        final_mass = OEW + aft_cargo_weight + front_cargo_weight  # Final mass

        front_load_cg.append((OEW_moment + aft_moment + front_moment) / (OEW + front_cargo_weight + aft_cargo_weight))
        aft_load_cg.append((OEW_moment + aft_moment + front_moment) / (OEW + front_cargo_weight + aft_cargo_weight))

        front_load_mass.append(final_mass)
        aft_load_mass.append(final_mass)

        return front_load_cg, aft_load_cg, front_load_mass, aft_load_mass, final_cg, final_mass

    def seat_loading(current_mass, current_cg):
        # Calculate seat loading CG points
        loading_mass = current_mass
        front_load_cg = [current_cg]
        aft_load_cg = [current_cg]
        front_load_mass = [current_mass]
        aft_load_mass = [current_mass]

        xcg_pax_front = xcg_front_seat
        xcg_pax_aft = xcg_front_seat + seat_pitch * (rows-1)

        # Front loading
        for i in range(rows ):
            xcg_calc = (front_load_cg[-1] * loading_mass + pax_weight * xcg_pax_front) / (loading_mass + pax_weight)
            loading_mass += pax_weight
            xcg_pax_front += seat_pitch
            front_load_cg.append(xcg_calc)
            front_load_mass.append(loading_mass)

        loading_mass = current_mass

        # Aft loading
        for i in range(rows ):
            xcg_calc = (aft_load_cg[-1] * loading_mass + pax_weight * xcg_pax_aft) / (loading_mass + pax_weight)
            loading_mass += pax_weight
            xcg_pax_aft -= seat_pitch
            aft_load_cg.append(xcg_calc)
            aft_load_mass.append(loading_mass)

        final_cg = front_load_cg[-1]
        final_mass = loading_mass

        return front_load_cg, aft_load_cg, front_load_mass, aft_load_mass, final_cg, final_mass

    def fuel_loading(current_mass, current_cg):
        loading_mass = [current_mass]
        load_cg = [current_cg]
        Fuel_ad = Fuel_max
        if current_mass + Fuel_max > MTOW:
            Fuel_ad = MTOW - current_mass

        xcg_calc = (current_mass * current_cg + Fuel_ad * xcg_fuel) / (current_mass + Fuel_ad)
        load_cg.append(xcg_calc)
        loading_mass.append(current_mass + Fuel_ad)

        final_cg = xcg_calc
        final_mass = current_mass + Fuel_ad

        return load_cg, loading_mass, final_cg, final_mass

    #Convert all required variables to MAC
    xcg_front_seat = convert_mac(xcg_front_seat)  # x position first seat from nose
    xcg_front_CB = convert_mac(xcg_front_CB)  # xcg front cargo bay
    xcg_aft_CB = convert_mac(xcg_aft_CB)  # xcg aft cargo bay
    xcg_fuel = convert_mac(xcg_fuel)  # xcg of the fuel tanks

    rows = int(pax / 4)

    OEW_moment = OEW * xcg_OEW  # Moment of OEW
    front_moment = xcg_front_CB * front_cargo_weight  # Moment of front cargo bay
    aft_moment = xcg_aft_CB * aft_cargo_weight  # Moment of aft cargo bay

    # Perform calculations
    mass = OEW
    xcg = xcg_OEW
    CB_front, CB_aft, CB_front_mass, CB_aft_mass, xcg, mass = cargoload(mass, xcg)
    PW_front, PW_aft, PW_front_mass, PW_aft_mass, xcg, mass = seat_loading(mass, xcg)
    PA_front, PA_aft, PA_front_mass, PA_aft_mass, xcg, mass = seat_loading(mass, xcg)
    Fuel_center, Fuel_mass, xcg, mass = fuel_loading(mass, xcg)

    # Make lists of most and most forward CG
    xcg_forward = CB_front + PW_front + PA_front + Fuel_center
    xcg_aft = CB_aft + PW_aft + PA_aft + Fuel_center

    # Calculate most forward and most aft CG
    min_cg = min(xcg_forward) * 0.98    #Add 2% margin
    max_cg = max(xcg_aft) * 1.02        #Add 2% margin

    max_xcg_mac = max_cg / MAC
    min_xcg_mac = min_cg / MAC

    print(f'Max cg = {max_cg} [m]')
    print(f'Min cg = {min_cg} [m]')
    print(f'Max cg = {max_xcg_mac} [%MAC]')
    print(f'Min cg = {min_xcg_mac} [%MAC]')

    # Convert CG points to MAC relative position
    CB_front = [(x) / MAC for x in CB_front]
    CB_aft = [(x) / MAC for x in CB_aft]
    PW_front = [(x) / MAC for x in PW_front]
    PW_aft = [(x) / MAC for x in PW_aft]
    PA_front = [(x) / MAC for x in PA_front]
    PA_aft = [(x) / MAC for x in PA_aft]
    Fuel_center = [(x) / MAC for x in Fuel_center]

    # Plot results
    if plotting:
        plt.figure()
        plt.plot(CB_front, CB_front_mass, 'o-', color='maroon', label='Cargo')
        plt.plot(CB_aft, CB_aft_mass, 'o-', color='crimson')
        plt.plot(PW_front, PW_front_mass, 'o-', color='teal', label='Window seats')
        plt.plot(PW_aft, PW_aft_mass, 'o-', color='darkblue')
        plt.plot(PA_front, PA_front_mass, 'o-', color='fuchsia', label='Aisle seats')
        plt.plot(PA_aft, PA_aft_mass, 'o-', color='blue')
        plt.plot(Fuel_center, Fuel_mass, 'o-', color='cyan', label='Fuel')
        plt.xlabel(r'CG location $\tilde{x_{cg}}$[mac]')
        plt.ylabel('Mass [kg]')
        plt.legend()
        name = ac_name + '_loading diagram'
        plt.savefig(name+'.pdf')
        plt.show()

    return max_xcg_mac, min_xcg_mac



#max_xcg_mac_ref, min_xcg_mac_ref = loading_diagram(reference_aircraft,plotting=True)
max_xcg_mac_mod, min_xcg_mac_mod = loading_diagram(modified_aircraft,plotting=True)








