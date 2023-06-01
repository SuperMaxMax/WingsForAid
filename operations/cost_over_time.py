import sys
sys.path.append("..")

# Start your import below this
import matplotlib.pyplot as plt
import math as m

class costs:
    def __init__(self):
        self.timemax    = 100    # years
        self.startAC    = 40    # aircraft
        self.ACextramin = 30     # aircraft per year
        self.ACextramax = 31    # aircraft per year
        self.dead       = 5.85688845510E-04 # chance of crach/sortie/ac
        self.Smaint     = 4.961363E-02  # chance of small maintenance/sortie/ac
        self.Lmaint     = 5.774406E-04   # chance of large maintenance/sortie/ac
        self.price      = 1.6   # price we ask per kg of paylaod delivered
        self.maxac      = 300   # maximum number of aircraft 
        self.ops_ext    = 58.4  # approximate percentage of EXT OPS costs (excl. manufact. and EOL)
        self.ops_gb     = 8.3   # approximate percentage of GB OPS costs (excl. manufact. and EOL)
        self.ops_sortie = 33.3  # approximate percentage of SORTIE OPS costs (excl. manufact. and EOL)
        self.EOL        = 5000  # per dead ac
        self.sortie     = 233000    #cost per sortie
        self.gb_cost    = 57870 # cost of one GB (30 ac)
        self.cost_ac    = 29640 # cost per ac
        self.acperop    = 30    # necessary ac per operation
        
        self.payl_box   = 20    # kg of payload per box
        self.box_w      = 3     # weight of box in kg
        self.boxes      = 12    # amount of boxes per sortie
        self.kgprice    = 4   # euros per kg

        self.sortie_pd  = 3.4148976
        self.days_perop = 28
        self.ops_year   = 4     # amount of operations (non-simulanious per year)
        #so if you have this as 4, and 2 simult. possible, that is 8 ops per year


# Per 30 ac at GB 27 are expected to be operational at any time on average
plot = True
c = costs()

def op(tot_ac):
            # ac_at_ops = tot_ac
            # sorties = ac_at_ops / c.sortie_pd
            # rev = sorties * c.boxes * c.kgprice
            deadac = tot_ac * c.dead *c.sortie_pd * c.days_perop
            # maint = ac_at_ops * c.Lmaint + ac_at_ops * c.Smaint
            # ac_ops_end  = ac_at_ops - m.ceil(deadac)
            print(deadac)
            return deadac #,rev

def rev(ops_year):
     sorties_p_op = c.sortie_pd * c.days_perop
     rev_p_op = c.kgprice * c.payl_box * c.boxes * sorties_p_op
     print(rev_p_op)
     rev_tot = rev_p_op *  ops_year
     print(rev_tot)
     return rev_tot

def calc(extra_ac):
    init_cost = c.cost_ac * c.startAC
    tot_ac = c.startAC
    totalcost = [init_cost]
    costtime = [init_cost]
    revenue = [0]
    totalrev = [0]
    ac = []
    # acend = []

    for i in range(c.timemax):
        ac.append(tot_ac)
        if tot_ac < 100:
            ops_parallel = 1
        else:
             ops_parallel = 2

        ops_series = min((tot_ac // 15),2)
        ops_year = ops_parallel * ops_series
        ac_at_ops = c.acperop * ops_year

        deadac = op(ac_at_ops)
        print(deadac)
        reven = rev(ops_year)

        revenue.append(reven)

        tot_ac = tot_ac - deadac
        # ac_togo = c.maxac - tot_ac

        if tot_ac < c.maxac:
            manu_cost = extra_ac * c.cost_ac
            # ac_togo = ac_togo - extra_ac
            tot_ac = tot_ac + extra_ac
            totalcost.append(manu_cost)
        
        costtime.append(sum(totalcost))
        totalrev.append(sum(revenue))
    return(costtime, totalrev,ac)

total = []
if plot:
    for i in range(c.ACextramin,(c.ACextramax+1)):

        costovertime, revovertime, actot = calc(i) 
        # plt.plot(range(c.timemax),actot,label='aircraft')
        plt.plot(range(c.timemax+1),costovertime, label=i)
        plt.plot(range(c.timemax+1),revovertime, label = 'revenue')
        for j in range(len(revovertime)):
            total.append(( revovertime[j] + costovertime[j] ))
        # plt.plot(range((c.timemax+1)*2), total, label = 'total')
                
    plt.legend()                                                      
    plt.show()                                                                                        




