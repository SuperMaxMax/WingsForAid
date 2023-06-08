import sys
sys.path.append("..")

# Start your import below this
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math as m

class costs:
    def __init__(self):
        self.timemax    = 11    # years
        self.startAC    = 30    # aircraft
        self.ACextramin = 30     # aircraft per year
        self.ACextramax = 31    # aircraft per year
        self.dead       = 5.85688845510E-04 # chance of crach/sortie/ac
        self.Smaint     = 4.961363E-02  # chance of small maintenance/sortie/ac
        self.Lmaint     = 5.774406E-05   # chance of large maintenance/sortie/ac
        self.price      = 1.6   # price we ask per kg of paylaod delivered
        self.maxac      = 300   # maximum number of aircraft 
        # self.ops_ext    = 58.4  # approximate percentage of EXT OPS costs (excl. manufact. and EOL)
        # self.ops_gb     = 8.3   # approximate percentage of GB OPS costs (excl. manufact. and EOL)
        # self.ops_sortie = 33.3  # approximate percentage of SORTIE OPS costs (excl. manufact. and EOL)
        self.EOL        = 250  # per dead ac
        # self.sortie     = 233000    #cost per sortie
        # self.gb_cost    = 57870 # cost of one GB (30 ac)
        self.cost_ac    = 29640 # cost per ac
        self.acperop    = 30    # necessary ac per operation
        
        self.payl_box   = 20    # kg of payload per box
        self.box_w      = 3     # weight of box in kg
        self.boxes      = 12    # amount of boxes per sortie
        self.kgprice    = [1.25,1.3,1.35]
        # [1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2]   # euros per kg
        # [1.25,1.3,1.35]

        self.sortie_pd  = 3.4148976
        self.days_perop = 28
        self.ops_year   = 4     # amount of operations (non-simulanious per year)
        #so if you have this as 4, and 2 simult. possible, that is 8 ops per year

        # self.maint_p_ac = 0.0004   # increase in costs/kg per aircraft (excl manufacturing)
        # self.cost_40ac  = 1.085631948   # cost but only OSP (excl manufacturing)

        self.cost_30ac  = 0.842259464  # cost per kg for the first 30 aircraft 
        self.extra_1ac = 0.008995361 / 30   # extra cost per aircraft, including manufacturing costs



# Per 30 ac at GB 27 are expected to be operational at any time on average
changing_ac = True
changing_price = True
roi_calc = False
c = costs()

def op(tot_ac):
            deadac = tot_ac * c.dead *c.sortie_pd * c.days_perop
            return deadac 

def inc(ops_year,price0):
     sorties_p_op = c.sortie_pd * c.days_perop
     rev_p_op = price0 * c.payl_box * c.boxes * sorties_p_op
     rev_tot = rev_p_op *  ops_year
     return rev_tot

def maint(ac_tot,ops_year):
    money = c.cost_30ac + c.extra_1ac * (ac_tot - 30) #euro per kg aid
    sorties_p_op = c.sortie_pd * c.days_perop
    per_op = money * c.payl_box * c.boxes * sorties_p_op #euros per op
    o = per_op * ops_year # euros per year
    return o


def calc(extra_ac,price):
    init_cost = c.cost_ac * c.startAC
    tot_ac = c.startAC
    totalcost = [init_cost]
    costtime = [init_cost]
    maintcost = []
    income = []
    totalinc = []
    ac = []

    for i in range(c.timemax):
        ac.append(tot_ac)

        ops_parallel = tot_ac // 50

        ops_series = min((tot_ac // 15),2)
        ops_year = ops_parallel * ops_series

        a = maint(tot_ac,ops_year)
        maintcost.append(a)
        incom = inc(ops_year,price)
        income.append(incom)

        ac_at_ops = c.acperop * ops_year
        deadac = op(ac_at_ops)
        # print("deadac = ",deadac)
        tot_ac = tot_ac - deadac
        ac_togo = c.maxac - tot_ac

        if tot_ac < c.maxac:
            manu_cost = extra_ac * c.cost_ac
            tot_ac = tot_ac + extra_ac
            totalcost.append(manu_cost)
        
        costtime.append(sum(totalcost))
        totalinc.append(sum(income))
        
    return(costtime, income,ac,maintcost)

################## PLOT STUFF DOWN BELOW ##############################

if changing_price:
    for i in c.kgprice:
        extra_ac = 40
        costovertime, income, actot, maintenance = calc(extra_ac,i)
        revenue = []
        roi = []
        totalrev = []
        for j in range(len(income)):
            if income[j] != 0:
                revenue.append(income[j]-maintenance[j])
                totalrev.append(sum(revenue))
                if income[j] != maintenance[j]:
                    roi.append(income[j]/maintenance[j])
                else:
                    roi.append(0)
        # plt.plot(range(1,c.timemax,1),revenue, label = i)
        plt.plot(range(1,c.timemax,1),totalrev, label = i)
        # plt.plot(range(1,c.timemax,1),roi, label = i)
    # plt.grid(True, which='both') 
    left, bottom, width, height = (0.9, 1.3, 9.2, 0.2)
    rect=mpatches.Rectangle((left,bottom),width,height, 
                        #fill=False,
                        alpha=0.1,
                       facecolor="green")
    plt.gca().add_patch(rect)
    # plt.plot(range(1,c.timemax,1),[1]*10, color='k', linestyle=':')
    # plt.axhline(y=1, color='k')
    # plt.axhline(y=1.5, color='k')
    # Reversing the order of legend labels
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles[::-1], labels[::-1])
    # plt.legend(loc='right')
    plt.show()


total = []
x = [10,20,30,40,50]
if changing_ac:
    for i in x:
    # for i in range(c.ACextramin,(c.ACextramax+1)):
        price = 1.1
        costovertime, income, actot, maintenance = calc(i,price)
        revenue = []
        for j in range(len(income)):
             revenue.append(income[j]-maintenance[j])
        # print(income - maintenance)
        # plt.plot(range(c.timemax),actot,label=i)
        # plt.plot(range(c.timemax+1),costovertime, label=i)
        # plt.plot(range(c.timemax+1),income, label = i)
        # plt.plot(range(c.timemax+1),maintenance, label = i)
        plt.plot(range(c.timemax),revenue, label = i)
        # for j in range(len(revovertime)):
        #     total.append(( revovertime[j] + costovertime[j] ))
        # plt.plot(range((c.timemax+1)*2), total, label = 'total')
                
    plt.legend()                                                      
    plt.show()                                                                                        




