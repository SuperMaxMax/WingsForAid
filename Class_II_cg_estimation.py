import matplotlib.pyplot as plt

def cg_calc(obj):
    # --- Wing group
    # Wing
    if obj.sweep_angle == 0:
        wing_cg = 0.4 * obj.rootchord                 # 40% of root chord plus Leading Edge location
    else:
        wing_cg = 0.4 * obj.rootchord                 # to be done later, depends on spar locations (table 8-15 Torenbeek)

    # Control surfaces
    control_surfaces_cg = obj.x_lemac + obj.MAC_length  # guess for now
    
    W_wing_gr = obj.W_w + obj.W_sc
    x_wcg = (wing_cg*obj.W_w + control_surfaces_cg*obj.W_sc)/(W_wing_gr)

    # --- Fuselage group # propellor to be done
    # Fuselage and engine
    prop_correction = 0.06                      # correction for propeller weight
    if obj.engine_pos == 'tractor':
        fus_cg = 0.45 * obj.l_f                 # educated guess
        engine_cg = 0.327 - prop_correction     # based on Rotax 912is (.g. or rotax 912is is at 327 mm, total length is 665.1 mm)
    elif obj.engine_pos == 'pusher':
        fus_cg = 0.55 * obj.l_f                 # educated guess
        engine_cg = obj.l_f - (0.6651-0.327) + prop_correction # based on Rotax 912is
    elif obj.engine_pos == 'fuselage':
        fus_cg = 0.53 * obj.l_f                 # educated guess
        engine_cg = 0.8 * obj.l_f               # educated guess

    # Tail
    if obj.boom == True:
        tail_cg = obj.l_f + 0.9*obj.l_f_boom    # educated guess
    else:
        tail_cg = 0.9*obj.l_f                   # for now at 0.9 of fuselage length

    # Equipment
    eq_cg = 0.5*obj.l_f                        # educated guess

    # Nacelle
    # at engine cg

    # Undercarriage
    # For now: cg assumed to be at aircraft cg -> not taken into account for X_FCG, but is part of OEW

    W_fus_gr = obj.W_fus + obj.W_pg + obj.W_t + obj.W_eq + obj.W_n + obj.W_uc
    X_FCG = (fus_cg*obj.W_fus + engine_cg*obj.W_pg + tail_cg*obj.W_t + eq_cg*obj.W_eq + engine_cg*obj.W_n)/(W_fus_gr - obj.W_uc)

    # X_LEMAC and xc_OEW
    xc_OEW = obj.xc_OEW_p*obj.MAC_length
    X_LEMAC = X_FCG + obj.MAC_length * ((x_wcg/obj.MAC_length)*(W_wing_gr/W_fus_gr)-(xc_OEW)*(1+W_wing_gr/W_fus_gr))

    # Final CG
    W_OEW = W_wing_gr+W_fus_gr
    X_OEW = X_LEMAC + xc_OEW
    print(f"W_OEW = {W_OEW} N, X_OEW = {X_OEW} m")

    # Fuel
    W_fuel_wi = obj.W_F
    X_fuel_wi = X_LEMAC + 0.5*obj.MAC_length
    print(f"W_fuel_wi = {W_fuel_wi} N, X_fuel_wi = {X_fuel_wi} m")

    # Payload
    if obj.engine_pos == 'tractor':
        dist_front = 0.6651 + 0.6  # [m]
    elif obj.engine_pos == 'pusher':
        dist_front = 0.6
    elif obj.engine_pos == 'fuselage':
        dist_front = 0.6
    
    # 2 boxes in front
    W_2box_f = 1/6*obj.W_PL
    X_2box_f = dist_front + 0.60/2
    # 4 boxes in the front
    W_4box_f = 1/3*obj.W_PL
    X_4box_f = dist_front + 2*0.60/2
    # 2 boxes in back
    W_2box_b = 1/6*obj.W_PL
    X_2box_b = dist_front + 5*0.6 + 0.6/2
    # 4 boxes in back
    W_4box_b = 1/3*obj.W_PL
    X_4box_b = dist_front + 4*0.6 + 2*0.6/2
    # all boxes
    W_allbox = obj.W_PL
    X_allbox = dist_front + 6*0.6/2

    # Calculate points to plot
    # OEW + fuel
    W_OEW_fuel_frac = (W_OEW + W_fuel_wi)/obj.W_TO
    X_OEW_fuel = (W_OEW*X_OEW + W_fuel_wi*X_fuel_wi)/(W_OEW + W_fuel_wi)
    print(f"W_OEW_fuel_frac = {W_OEW_fuel_frac}, X_OEW_fuel = {X_OEW_fuel}")
    # OEW + fuel + 2 boxes in front
    W_OEW_fuel_2box_f_frac = W_OEW_fuel_frac + W_2box_f/obj.W_TO
    X_OEW_fuel_2box_f = (W_OEW_fuel_frac*X_OEW_fuel + W_2box_f*X_2box_f)/(W_OEW_fuel_frac + W_2box_f)
    print(f"W_OEW_fuel_2box_f_frac = {W_OEW_fuel_2box_f_frac}, X_OEW_fuel_2box_f = {X_OEW_fuel_2box_f}")
    # OEW + fuel + 4 boxes in front
    W_OEW_fuel_4box_f_frac = W_OEW_fuel_frac + W_4box_f/obj.W_TO
    X_OEW_fuel_4box_f = (W_OEW_fuel_frac*X_OEW_fuel + W_4box_f*X_4box_f)/(W_OEW_fuel_frac + W_4box_f)
    print(f"W_OEW_fuel_4box_f_frac = {W_OEW_fuel_4box_f_frac}, X_OEW_fuel_4box_f = {X_OEW_fuel_4box_f}")
    # OEW + fuel + 2 boxes in back
    W_OEW_fuel_2box_b_frac = W_OEW_fuel_frac + W_2box_b/obj.W_TO
    X_OEW_fuel_2box_b = (W_OEW_fuel_frac*X_OEW_fuel + W_2box_b*X_2box_b)/(W_OEW_fuel_frac + W_2box_b)
    print(f"W_OEW_fuel_2box_b_frac = {W_OEW_fuel_2box_b_frac}, X_OEW_fuel_2box_b = {X_OEW_fuel_2box_b}")
    # OEW + fuel + 4 boxes in back
    W_OEW_fuel_4box_b_frac = W_OEW_fuel_frac + W_4box_b/obj.W_TO
    X_OEW_fuel_4box_b = (W_OEW_fuel_frac*X_OEW_fuel + W_4box_b*X_4box_b)/(W_OEW_fuel_frac + W_4box_b)
    print(f"W_OEW_fuel_4box_b_frac = {W_OEW_fuel_4box_b_frac}, X_OEW_fuel_4box_b = {X_OEW_fuel_4box_b}")
    # OEW + fuel + all boxes
    W_OEW_fuel_allbox_frac = W_OEW_fuel_frac + W_allbox/obj.W_TO
    X_OEW_fuel_allbox = (W_OEW_fuel_frac*X_OEW_fuel + W_allbox*X_allbox)/(W_OEW_fuel_frac + W_allbox)
    print(f"W_OEW_fuel_allbox_frac = {W_OEW_fuel_allbox_frac}, X_OEW_fuel_allbox = {X_OEW_fuel_allbox}")

    # Plot each point
    xs = [X_OEW, X_fuel_wi, X_2box_f, X_4box_f, X_2box_b, X_4box_b, X_allbox]
    w_fracs = [W_OEW_fuel_frac, W_OEW_fuel_frac, W_OEW_fuel_2box_f_frac, W_OEW_fuel_4box_f_frac, W_OEW_fuel_2box_b_frac, W_OEW_fuel_4box_b_frac, W_OEW_fuel_allbox_frac]
    labels = ['OEW', 'OEW + Fuel', 'OEW + Fuel + 2 boxes front', 'OEW + Fuel + 4 boxes front', 'OEW + Fuel + 2 boxes back', 'OEW + Fuel + 4 boxes back', 'OEW + Fuel + all boxes']
    # plot points with labels
    for x, w, label in zip(xs, w_fracs, labels):
        plt.scatter(x, w, label=label)
    plt.xlabel('X_cg [m]')
    plt.ylabel('Mass fraction [-]')
    plt.grid()
    plt.legend()
    plt.title(f'Mass fraction vs X_cg for {obj.name}')
    plt.show()

    # Save most forward and most aft and fully loaded c.g. in object
    obj.X_cg_fwd = min(xs)
    obj.X_cg_aft = max(xs)
    obj.X_cg_range = obj.X_cg_aft - obj.X_cg_fwd
    obj.X_cg_full = X_OEW_fuel_allbox