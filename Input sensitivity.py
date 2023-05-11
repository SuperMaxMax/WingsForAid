import numpy as np

def random_value(value, perc):
    """
    Input a value and the percentage randomness
    """
    if type(value) == np.ndarray:
        random_array = np.random.uniform(-1, 1, value.shape)
        random_array = (1+random_array*perc)*value
        random_array = (1-sum(random_array))/len(random_array) + random_array
        return random_array
    else:
        return (1+np.random.uniform(-1, 1)*perc)*value

def total_weights(x, perc):
    """
    Make a vector containing the total weights used in trade-off
    """
    
    tot_weights = []

    for i in range(x):
        # Everything cost side
        cost = random_value(0.9, perc)

        cost_top_r = random_value(cost_top, perc)
        cost_perf_r = random_value(cost_perf, perc)
        cost_oper_r = random_value(cost_oper, perc)
        cost_manu_r = random_value(cost_manu, perc)

        # Everything environment side
        env = 1-cost

        env_top_r = random_value(env_top, perc)
        env_perf_r = random_value(env_perf, perc)
        env_oper_r = random_value(env_oper, perc)
        env_manu_r = random_value(env_manu, perc)
    
        # Total
        cost_tot = cost*cost_top_r
        cost_tot_fin = np.hstack((cost_tot[0]*cost_perf_r, cost_tot[1], cost_tot[2]*cost_oper_r, cost_tot[3]*cost_manu_r))

        env_tot = env*env_top_r
        env_tot_fin = np.hstack((env_tot[0]*env_perf_r, env_tot[1], env_tot[2]*env_oper_r, env_tot[3]*env_manu_r))

        tot_weights.append(cost_tot_fin+env_tot_fin)
    
    return tot_weights

# Define the weights
cost_top = np.array([0.35, 0.1, 0.35, 0.2])
cost_perf = np.array([0.4, 0.2, 0.1, 0.1, 0.2])
cost_oper = np.array([0.3, 0.1, 0.15, 0.14, 0.01, 0.1, 0.2])
cost_manu = np.array([0.8, 0.2])

env_top = np.array([0.3, 0.15, 0.25, 0.3])
env_perf = np.array([0.4, 0.05, 0.05, 0.15, 0.35])
env_oper = np.array([0.25, 0.05, 0.1, 0.15, 0.1, 0.15, 0.2])
env_manu = np.array([0.7, 0.3])

input = total_weights(3, 0)
print(input)
