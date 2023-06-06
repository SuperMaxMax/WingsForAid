def fuelusesortie(ac_obj, atm_obj, W_F, n_drops, n_boxes, cruiseheight, V_cruise = None, Range = None, dropregion = None, V_des = None):
    if Range == None:
        Range = ac_obj.R / 2
    else:
        Range *= 1000
    if V_des == None:
        CL_des_opt = np.sqrt(ac_obj.CD0*np.pi*ac_obj.A*ac_obj.e)
    else:
        V_des = V_des
    if V_cruise == None:
        CL_cr_opt = np.sqrt(ac_obj.CD0*np.pi*ac_obj.A*ac_obj.e)
    else:
        V_cruise = V_cruise
    # determine take-off weight
    W_TO = ac_obj.W_OE + W_F + n_boxes * ac_obj.boxweight
    W = W_TO
    print("=====================================================")
    print(f"Take-off weight: {W_TO} [kg] | OEW: {ac_obj.W_OE} [kg] | Fuelweight: {W_F} [kg] / {W_F/ac_obj.fueldensity} [L] | Payload: {n_boxes} boxes / {n_boxes*ac_obj.boxweight} [kg]")
    # determine weight right after take off using fuel fractions
    W_a_TO = W_TO * ac_obj.W1W_TO * ac_obj.W2W1 * ac_obj.W3W2
    W = W_a_TO
    W_F_used = W_TO - W_a_TO
    W_F -= W_F_used
    # W_a_TO happens at screen height h = 15 m
    h = 15.0         
    p, rho = atm_parameters(atm_obj, h)[0], atm_parameters(atm_obj, h)[2]
    print(p, rho, "at screenheight")
    V = 1.3 * np.sqrt(2*W/(rho*ac_obj.Sw*ac_obj.CL_max_TO))
    dt= 1.0
    t = 0.0
    x = 0.0
    while h < ac_obj.accelheight:
        Pa = ac_obj.power * ac_obj.prop_eff * p/(atm_obj.p0) * hp_to_watt           # 100% power setting
        CL = 2*W*atm_obj.g/(rho*ac_obj.Sw*V**2)
        CD = dragpolar(ac_obj, CL)
        Pr = 1/2 * rho * V**3 * ac_obj.Sw * CD
        ROC= (Pa - Pr)/(W*atm_obj.g)
        gamma = ROC / V
        x  += V * np.cos(gamma) * dt
        h  += ROC * dt
        W_F_used += (Pa/ac_obj.prop_eff) * ac_obj.SFC
        W  -= (Pa/ac_obj.prop_eff) * ac_obj.SFC
        p, rho = atm_parameters(atm_obj, h)[0], atm_parameters(atm_obj, h)[2]
        print(p, rho, "climb to accel height")
        t += dt
    print(rho)
    CL_opt_climb = np.sqrt(3*ac_obj.CD0*np.pi*ac_obj.A*ac_obj.e)
    V_opt = np.sqrt(2*W*atm_obj.g/(rho*ac_obj.Sw*CL_opt_climb))
    while V < V_opt:
        print(V)
        Pa = ac_obj.power * ac_obj.prop_eff * p/(atm_obj.p0) * hp_to_watt
        print(Pa)
        CL = 2*W*atm_obj.g/(rho*ac_obj.Sw*V**2)
        CD = dragpolar(ac_obj, CL)
        Pr = 1/2 * rho * V**3 * ac_obj.Sw * CD
        print(Pr)
        x  += V*dt
        dVdt = (1/(V*W))(Pa - Pr)
        print(dVdt)
        V += dVdt * dt
        print(V)
        W_F_used += (Pa/ac_obj.prop_eff) * ac_obj.SFC
        W  -= (Pa/ac_obj.prop_eff) * ac_obj.SFC
        t += dt
    print(V, V_opt, "assess whether equal")
    if dropregion == None:
        target_dist = Range / n_drops
    
        if 0.0 <= target_dist <= 10000.0:
            h_cruise = 1000 * 0.3048                                                              # m
        elif 10000.0 <= target_dist <= 50000.0:
            h_cruise = 1000 + (target_dist - 10000.0) * 0.1
        elif 50000.0 <= target_dist <= 100000.0:
            h_cruise = 5000 + (target_dist - 50000.0) * 0.05
        elif 100000.0 <= target_dist <= 200000.0:
            h_cruise = 7500 + (target_dist - 100000.0) * 0.025
        else:
            h_cruise = cruiseheight
        h_cruise = round(h_cruise/500) * 500 * 0.3048
        for i in range(n_drops):
            while h < h_cruise:
                Pa = ac_obj.power * ac_obj.prop_eff * p/(atm_obj.p0) * hp_to_watt           # 100% power setting
                CL = 2*W*atm_obj.g/(rho*ac_obj.Sw*V**2)
                CD = dragpolar(ac_obj, CL)
                Pr = 1/2 * rho * V**3 * ac_obj.Sw * CD
                ROC= (Pa - Pr)/(W*atm_obj.g)
                print(ROC, "Rate of Climb")
                gamma = ROC / V
                x  += V * np.cos(gamma) * dt
                h  += ROC * dt
                W_F_used += (Pa/ac_obj.prop_eff) * ac_obj.SFC
                W  -= (Pa/ac_obj.prop_eff) * ac_obj.SFC
                p, rho = atm_parameters(atm_obj, h)[0], atm_parameters(atm_obj, h)[2]
                print(p, rho, "optimum climb, dropregion = None")
                t += dt
            V_cruise = 2*W*atm_obj.g/(rho*ac_obj.Sw*CL_cr_opt)
            while V < V_cruise:
                Pa = ac_obj.power * ac_obj.prop_eff * p/(atm_obj.p0) * hp_to_watt
                CL = 2*W*atm_obj.g/(rho*ac_obj.Sw*V**2)
                CD = dragpolar(ac_obj, CL)
                Pr = 1/2 * rho * V**3 * ac_obj.Sw * CD
                x  += V*dt
                dVdt = (1/V)*W*(Pa - Pr)
                V += dVdt * dt
                W_F_used += (Pa/ac_obj.prop_eff) * ac_obj.SFC
                W  -= (Pa/ac_obj.prop_eff) * ac_obj.SFC
                t += dt
            LD_Max = np.sqrt((np.pi*ac_obj.A*ac_obj.e)/(2*ac_obj.CD0))
            TOD_dtt = LD_Max * h_cruise
            x_start_des = target_dist - TOD_dtt
            p, rho = atm_parameters(atm_obj, h_cruise)[0], atm_parameters(atm_obj, h_cruise)[2]
            print(p, rho, "cruise, dropregion = None")
            while x < x_start_des:
                if V_cruise == None:                                                        # optimal case
                    CL = CL_cr_opt
                else: 
                    CL = 2*W*atm_obj.g/(rho*ac_obj.Sw*V_cruise**2)
                CD = dragpolar(ac_obj, h)
                Pr = 1/2 * rho * V**3 * ac_obj.Sw * CD
                FF = (Pr/ac_obj.prop_eff) * ac_obj.SFC
                x += V*dt
                W_F_used += FF * dt
                W -= FF * ac_obj.SFC
                t += dt
            FF_idle = 1/3600                                                                # 1 kg of fuel per hour in idle
            while h > 15.0:
                CL = CL_des_opt
                V  = 2*W*atm_obj.g/(rho*ac_obj.Sw*CL)
                CD = dragpolar(ac_obj, CL)
                Pr  = 1/2 * rho * V**3 * ac_obj.Sw * CD
                ROD = -Pr / (W*atm_obj.g)
                gamma = ROD/V
                h  += ROD
                p, rho = atm_parameters(atm_obj, h_cruise)[0], atm_parameters(atm_obj, h_cruise)[2]
                print(p, rho, "descent, dropregion = None")
                x += V*np.cos(gamma)*dt
                W_F_used += FF_idle * dt
                W -= FF_idle * dt
                t += dt
    else:
        target_dist = Range - dropregion * 1000
        target_dist = Range / n_drops
    
        if 0.0 <= target_dist <= 10000.0:
            h_cruise = 1000 * 0.3048                                                              # m
        elif 10000.0 <= target_dist <= 50000.0:
            h_cruise = 1000 + (target_dist - 10000.0) * 0.1
        elif 50000.0 <= target_dist <= 100000.0:
            h_cruise = 5000 + (target_dist - 50000.0) * 0.05
        elif 100000.0 <= target_dist <= 200000.0:
            h_cruise = 7500 + (target_dist - 100000.0) * 0.025
        else:
            h_cruise = cruiseheight
        h_cruise = round(h_cruise/500) * 500 * 0.3048
        for i in range(n_drops):
            while h < h_cruise:
                Pa = ac_obj.power * ac_obj.prop_eff * p/(atm_obj.p0) * hp_to_watt           # 100% power setting
                CL = 2*W*atm_obj.g/(rho*ac_obj.Sw*V**2)
                CD = dragpolar(ac_obj, CL)
                Pr = 1/2 * rho * V**3 * ac_obj.Sw * CD
                ROC= (Pa - Pr)/(W*atm_obj.g)
                gamma = ROC / V
                x  += V * np.cos(gamma) * dt
                h  += ROC * dt
                W_F_used += (Pa/ac_obj.prop_eff) * ac_obj.SFC
                W  -= (Pa/ac_obj.prop_eff) * ac_obj.SFC
                p, rho = atm_parameters(atm_obj, h)[0], atm_parameters(atm_obj, h)[2]
                print(p, rho, "optimum climb, dropregion defined")
                t += dt
            while V < V_cruise:
                Pa = ac_obj.power * ac_obj.prop_eff * p/(atm_obj.p0) * hp_to_watt
                CL = 2*W*atm_obj.g/(rho*ac_obj.Sw*V**2)
                CD = dragpolar(ac_obj, CL)
                Pr = 1/2 * rho * V**3 * ac_obj.Sw * CD
                x  += V*dt
                dVdt = (1/V)*W*(Pa - Pr)
                V += dVdt * dt
                W_F_used += (Pa/ac_obj.prop_eff) * ac_obj.SFC
                W  -= (Pa/ac_obj.prop_eff) * ac_obj.SFC
                t += dt
            LD_Max = np.sqrt((np.pi*ac_obj.A*ac_obj.e)/(2*ac_obj.CD0))
            TOD_dtt = LD_Max * h_cruise
            x_start_des = target_dist - TOD_dtt
            p, rho = atm_parameters(atm_obj, h_cruise)[0], atm_parameters(atm_obj, h_cruise)[2]
            while x < x_start_des:
                if V_cruise == None:                                                        # optimal case
                    CL = CL_cr_opt
                else: 
                    CL = 2*W*atm_obj.g/(rho*ac_obj.Sw*V_cruise**2)
                CD = dragpolar(ac_obj, h)
                Pr = 1/2 * rho * V**3 * ac_obj.Sw * CD
                FF = (Pr/ac_obj.prop_eff) * ac_obj.SFC
                x += V*dt
                W_F_used += FF * dt
                W -= FF * ac_obj.SFC
                t += dt
            FF_idle = 1/3600                                                                # 1 kg of fuel per hour in idle
            while h > 15.0:
                CL = CL_des_opt
                V  = 2*W*atm_obj.g/(rho*ac_obj.Sw*CL)
                CD = dragpolar(ac_obj, CL)
                Pr  = 1/2 * rho * V**3 * ac_obj.Sw * CD
                ROD = -Pr / (W*atm_obj.g)
                gamma = ROD/V
                h  += ROD
                p, rho = atm_parameters(atm_obj, h_cruise)[0], atm_parameters(atm_obj, h_cruise)[2]
                print(p, rho)
                x += V*np.cos(gamma)*dt
                W_F_used += FF_idle * dt
                W -= FF_idle * dt
                t += dt
        # now for the drops in the dropregion
        d_interdrop = dropregion*1000 / n_drops
        if 0.0 <= d_interdrop <= 10000.0:
            h_interdrop = 1000 * 0.3048                                                              # m
        elif 10000.0 <= d_interdrop <= 50000.0:
            h_interdrop = 1000 + (d_interdrop - 10000.0) * 0.1
        elif 50000.0 <= d_interdrop <= 100000.0:
            h_interdrop = 5000 + (d_interdrop - 50000.0) * 0.05
        elif 100000.0 <= d_interdrop <= 200000.0:
            h_interdrop = 7500 + (d_interdrop - 100000.0) * 0.025
        h_interdrop = round(h_interdrop/500) * 500 * 0.3048
        for i in range(n_drops):
            x_inter = 0.0
            while h < h_interdrop:
                Pa = ac_obj.power * ac_obj.prop_eff * p/(atm_obj.p0) * hp_to_watt           # 100% power setting
                CL = 2*W*atm_obj.g/(rho*ac_obj.Sw*V**2)
                CD = dragpolar(ac_obj, CL)
                Pr = 1/2 * rho * V**3 * ac_obj.Sw * CD
                ROC= (Pa - Pr)/(W*atm_obj.g)
                gamma = ROC / V
                x  += V * np.cos(gamma) * dt
                x_inter  += V * np.cos(gamma) * dt
                h  += ROC * dt
                W_F_used += (Pa/ac_obj.prop_eff) * ac_obj.SFC
                W  -= (Pa/ac_obj.prop_eff) * ac_obj.SFC
                p, rho = atm_parameters(atm_obj, h)[0], atm_parameters(atm_obj, h)[2]
                print(p, rho)
                t += dt
            while V < V_cruise:
                Pa = ac_obj.power * ac_obj.prop_eff * p/(atm_obj.p0) * hp_to_watt
                CL = 2*W*atm_obj.g/(rho*ac_obj.Sw*V**2)
                CD = dragpolar(ac_obj, CL)
                Pr = 1/2 * rho * V**3 * ac_obj.Sw * CD
                x  += V*dt
                x_inter  += V * np.cos(gamma) * dt
                dVdt = (1/V)*W*(Pa - Pr)
                V += dVdt * dt
                W_F_used += (Pa/ac_obj.prop_eff) * ac_obj.SFC
                W  -= (Pa/ac_obj.prop_eff) * ac_obj.SFC
                t += dt
            LD_Max = np.sqrt((np.pi*ac_obj.A*ac_obj.e)/(2*ac_obj.CD0))
            TOD_dtt = LD_Max * h_interdrop
            x_inter_start_des = d_interdrop - TOD_dtt
            p, rho = atm_parameters(atm_obj, h_cruise)[0], atm_parameters(atm_obj, h_cruise)[2]
            print(p, rho)
            while x_inter < x_inter_start_des:
                if V_cruise == None:                                                        # optimal case
                    CL = CL_cr_opt
                else: 
                    CL = 2*W*atm_obj.g/(rho*ac_obj.Sw*V_cruise**2)
                CD = dragpolar(ac_obj, h)
                Pr = 1/2 * rho * V**3 * ac_obj.Sw * CD
                FF = (Pr/ac_obj.prop_eff) * ac_obj.SFC
                x += V*dt
                W_F_used += FF * dt
                W -= FF * ac_obj.SFC
                t += dt
            FF_idle = 1/3600                                                                # 1 kg of fuel per hour in idle
            while h > 15.0:
                CL = CL_des_opt
                V  = 2*W*atm_obj.g/(rho*ac_obj.Sw*CL)
                CD = dragpolar(ac_obj, CL)
                Pr  = 1/2 * rho * V**3 * ac_obj.Sw * CD
                ROD = -Pr / (W*atm_obj.g)
                gamma = ROD/V
                h  += ROD
                p, rho = atm_parameters(atm_obj, h_cruise)[0], atm_parameters(atm_obj, h_cruise)[2]
                print(p, rho)
                x += V*np.cos(gamma)*dt
                W_F_used += FF_idle * dt
                W -= FF_idle * dt
                t += dt
    # Return to base
    if 0.0 <= Range <= 10000.0:
            h_cruise = 1000 * 0.3048                                                              # m
    elif 10000.0 <= Range <= 50000.0:
        h_cruise = 1000 + (Range - 10000.0) * 0.1
    elif 50000.0 <= Range <= 100000.0:
        h_cruise = 5000 + (Range - 50000.0) * 0.05
    elif 100000.0 <= Range <= 200000.0:
        h_cruise = 7500 + (Range - 100000.0) * 0.025
    else:
        h_cruise = cruiseheight
    h_cruise = round(h_cruise/500) * 500 * 0.3048
    while h < h_cruise:
        Pa = ac_obj.power * ac_obj.prop_eff * p/(atm_obj.p0) * hp_to_watt           # 100% power setting
        CL = 2*W*atm_obj.g/(rho*ac_obj.Sw*V**2)
        CD = dragpolar(ac_obj, CL)
        Pr = 1/2 * rho * V**3 * ac_obj.Sw * CD
        ROC= (Pa - Pr)/(W*atm_obj.g)
        gamma = ROC / V
        x  += V * np.cos(gamma) * dt
        h  += ROC * dt
        W_F_used += (Pa/ac_obj.prop_eff) * ac_obj.SFC
        W  -= (Pa/ac_obj.prop_eff) * ac_obj.SFC
        p, rho = atm_parameters(atm_obj, h)[0], atm_parameters(atm_obj, h)[2]
        print(p, rho)
        t += dt
    while V < V_cruise:
        Pa = ac_obj.power * ac_obj.prop_eff * p/(atm_obj.p0) * hp_to_watt
        CL = 2*W*atm_obj.g/(rho*ac_obj.Sw*V**2)
        CD = dragpolar(ac_obj, CL)
        Pr = 1/2 * rho * V**3 * ac_obj.Sw * CD
        x  += V*dt
        dVdt = (1/V)*W*(Pa - Pr)
        V += dVdt * dt
        W_F_used += (Pa/ac_obj.prop_eff) * ac_obj.SFC
        W  -= (Pa/ac_obj.prop_eff) * ac_obj.SFC
        t += dt
    LD_Max = np.sqrt((np.pi*ac_obj.A*ac_obj.e)/(2*ac_obj.CD0))
    TOD_dtt = LD_Max * h_cruise
    x_start_des = target_dist - TOD_dtt
    p, rho = atm_parameters(atm_obj, h_cruise)[0], atm_parameters(atm_obj, h_cruise)[2]
    print(p, rho)
    while x < x_start_des:
        if V_cruise == None:                                                        # optimal case
            CL = CL_cr_opt
        else: 
            CL = 2*W*atm_obj.g/(rho*ac_obj.Sw*V_cruise**2)
        CD = dragpolar(ac_obj, h)
        Pr = 1/2 * rho * V**3 * ac_obj.Sw * CD
        FF = (Pr/ac_obj.prop_eff) * ac_obj.SFC
        x += V*dt
        W_F_used += FF * dt
        W -= FF * ac_obj.SFC
        t += dt
    FF_idle = 1/3600                                                                # 1 kg of fuel per hour in idle
    while h > 15.0:
        CL = CL_des_opt
        V  = 2*W*atm_obj.g/(rho*ac_obj.Sw*CL)
        CD = dragpolar(ac_obj, CL)
        Pr  = 1/2 * rho * V**3 * ac_obj.Sw * CD
        ROD = -Pr / (W*atm_obj.g)
        gamma = ROD/V
        h  += ROD
        p, rho = atm_parameters(atm_obj, h_cruise)[0], atm_parameters(atm_obj, h_cruise)[2]
        x += V*np.cos(gamma)*dt
        W_F_used += FF_idle * dt
        W -= FF_idle * dt
        t += dt
    W_final  = W * ac_obj.WfinalW10
    W_F_used += W - W_final
    # print(f"================================= Sortie Summary =================================")
    # print(f"Boxes delivered: {n_boxes} boxes / {n_boxes*ac_obj.boxweight} [kg]")
    # print(f"Distance flown: {np.round(x/1000, 2)} [km]")
    # print(f"Fuel used: {W_F_used} [kg]")
    # print(f"==================================================================================")
    return W_F_used

fuelusesortie(aircraft, atm, 70, 1, 12, 3048)