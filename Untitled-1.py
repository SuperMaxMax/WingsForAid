print("hello peeps")

# for i in range(len(CL_max_clean)):
    #     WS_stall_max = stallWS(V_s_min, rho, CL_max_clean[i])
    #     WS_stall_max = np.full(len(WP), WS_stall_max)
    #     WS_stall_max = np.vstack((WS_stall_max, WP))
    #     if i == 0:
    #         WS_stall = WS_stall_max
    #     else:
    #         WS_stall = np.vstack((WS_stall, WS_stall_max))

# WP_TO   = np.empty(0)
    # for j in range(len(CL_TO)):
    #     WP_takeoff  = (TOP_req/WS)*CL_TO[j]*sigma
    #     WP_takeoff  = np.vstack((WP_takeoff, WS))
    #     if j == 0:
    #         WP_TO   = WP_takeoff
    #     else:
    #         WP_TO   = np.vstack((WP_TO, WP_takeoff))

# WS_landing = np.empty(0)
    # for k in range(len(CL_LDG)): 
    #     WS_ldg = (CL_LDG[i]*rho*(LDG_dist/0.5915))/(2*f)
    #     WS_ldg = np.full(len(WP), WS_ldg)
    #     WS_ldg = np.vstack((WS_ldg, WP))
    #     if k == 0:
    #         WS_landing = WS_ldg
    #     else:
    #         WS_landing = np.vstack((WS_landing, WS_ldg))

# for z in range(int(len(WS_stall)/2)):
    #     if WS_stall[2*z][0] < WS_landing[2*z][0]:
    #         stall_limit += 1
    #     else:
    #         landing_limit += 1
    # if stall_limit > landing_limit:
    #     plot_stall_WS = True
    # elif stall_limit == landing_limit:
    #     plot_both_WS = True
    # else:
    #     plot_LDG_WS = True