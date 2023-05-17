import numpy as np

### Scores ###
pull_fuselage = np.array([3, 3, 3, 3, 3])
pull_boom = np.array([2, 3, 3, 3, 4])
push_fuselage = np.array([1, 2, 1, 4, 4])
push_boom = np.array([1, 2, 1, 4, 5])
top_fuselage = np.array([2, 3, 1, 2, 2])
top_boom = np.array([2, 3, 1, 2, 2])

original_weights = np.array([0.15, 0.2, 0.25, 0.1, 0.3])
iterations = 10000

#define storage of all weights
all_random_winner = np.zeros((6,iterations))

for i in range(len(original_weights)):
    for j in range(iterations):    
        #defining new weight array that will be outputted
        new_weight = np.zeros((1,5))[0]

        #random value
        mu = original_weights[i] #mean
        sigma = mu*0.2 #standard deviation
        random_weight = np.random.normal(mu, sigma, 1)[0]
        new_weight[i] = random_weight
        
        #scalling other weights
        difference = original_weights[i] - random_weight
        proportion = original_weights[i-1] + original_weights[i-2] + original_weights[i-3] + original_weights[i-4]

        new_weight[i-1] = original_weights[i-1] + difference * original_weights[i-1]/proportion
        new_weight[i-2] = original_weights[i-2] + difference * original_weights[i-2]/proportion
        new_weight[i-3] = original_weights[i-3] + difference * original_weights[i-3]/proportion
        new_weight[i-4] = original_weights[i-4] + difference * original_weights[i-4]/proportion

        scores = []
        scores.append(np.dot(np.transpose(pull_fuselage), new_weight))
        scores.append(np.dot(np.transpose(pull_boom), new_weight)) 
        scores.append(np.dot(np.transpose(push_fuselage), new_weight))
        scores.append(np.dot(np.transpose(push_boom), new_weight))
        scores.append(np.dot(np.transpose(top_fuselage), new_weight))
        scores.append(np.dot(np.transpose(top_boom), new_weight))

        winner = np.array(scores).argmax() + 1

        all_random_winner[i][j] = winner

print('-----------------------------------------------------------')
print('the overall distribution is:')
print(f'pull fuselage wins {np.count_nonzero(all_random_winner == 1)} out of {5 * iterations} times, which is {np.count_nonzero(all_random_winner == 1)/50000 * 100}%')
print(f'pull boom wins {np.count_nonzero(all_random_winner == 2)} out of {5 * iterations} times, which is {np.count_nonzero(all_random_winner == 2)/50000 * 100}%')
print(f'push fuselage wins {np.count_nonzero(all_random_winner == 3)} out of {5 * iterations} times, which is {np.count_nonzero(all_random_winner == 3)/50000 * 100}%')
print(f'push boom wins {np.count_nonzero(all_random_winner == 4)} out of {5 * iterations} times, which is {np.count_nonzero(all_random_winner == 4)/50000 * 100}%')
print(f'top fuselage wins {np.count_nonzero(all_random_winner == 5)} out of {5 * iterations} times, which is {np.count_nonzero(all_random_winner == 5)/50000 * 100}%')
print(f'top boom wins {np.count_nonzero(all_random_winner == 6)} out of {5 * iterations} times, which is {np.count_nonzero(all_random_winner == 6)/50000 * 100}%')
print('-----------------------------------------------------------')