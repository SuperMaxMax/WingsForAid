from parameters import UAV

concept = UAV('naam', engine_pos='tractor', boom = True)


def VC_lim_low(obj):
    return 33 * (obj.WS * 0.020885) ** 0.5

def VD_lim_low(obj):
    return 1.4 * VC_lim_low(obj)

def VA_lim_low(obj):
    return 


