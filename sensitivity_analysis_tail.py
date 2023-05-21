import numpy as np

### scores ###
conventional = np.array([2, 2, 2, 2])
v_tail = np.array([1, 2, 2, 1])
inverted_v_tail = np.array([3, 2, 2, 1])
y_tail = np.array([3, 1, 1, 2])
h_tail = np.array([3, 1, 1, 2])
flying_tailplane = np.array([2, 3, 2, 3])

original_weights = np.array([0.3, 0.15, 0.35, 0.2])
iterations = 10000

#define storage of all weights
all_random_winner = np.zeros((6,iterations))

for i in range(len(original_weights)):
    for j in range(iterations):    
        #defining new weight array that will be outputted
        new_weight = np.zeros((1,4))[0]

        #random value
        mu = original_weights[i] #mean
        sigma = mu*0.2 #standard deviation
        random_weight = np.random.normal(mu, sigma, 1)[0]
        new_weight[i] = random_weight
        
        #scalling other weights
        difference = original_weights[i] - random_weight
        proportion = original_weights[i-1] + original_weights[i-2] + original_weights[i-3]

        new_weight[i-1] = original_weights[i-1] + difference * original_weights[i-1]/proportion
        new_weight[i-2] = original_weights[i-2] + difference * original_weights[i-2]/proportion
        new_weight[i-3] = original_weights[i-3] + difference * original_weights[i-3]/proportion

        scores = []
        scores.append(np.dot(np.transpose(conventional), new_weight))
        scores.append(np.dot(np.transpose(v_tail), new_weight)) 
        scores.append(np.dot(np.transpose(inverted_v_tail), new_weight))
        scores.append(np.dot(np.transpose(y_tail), new_weight))
        scores.append(np.dot(np.transpose(h_tail), new_weight))
        scores.append(np.dot(np.transpose(flying_tailplane), new_weight))

        winner = np.array(scores).argmax() + 1

        all_random_winner[i][j] = winner

print('-----------------------------------------------------------')
print('the overall distribution is:')
print(f'convential wins {np.count_nonzero(all_random_winner == 1)} out of {4 * iterations} times, which is {np.count_nonzero(all_random_winner == 1)/40000 * 100}%')
print(f'v tail wins {np.count_nonzero(all_random_winner == 2)} out of {4 * iterations} times, which is {np.count_nonzero(all_random_winner == 2)/40000 * 100}%')
print(f'inverted v tail wins {np.count_nonzero(all_random_winner == 3)} out of {4 * iterations} times, which is {np.count_nonzero(all_random_winner == 3)/40000 * 100}%')
print(f'y tail wins {np.count_nonzero(all_random_winner == 4)} out of {4 * iterations} times, which is {np.count_nonzero(all_random_winner == 4)/40000 * 100}%')
print(f'h tail wins {np.count_nonzero(all_random_winner == 5)} out of {4 * iterations} times, which is {np.count_nonzero(all_random_winner == 5)/40000 * 100}%')
print(f'flying tailplane wins {np.count_nonzero(all_random_winner == 6)} out of {4 * iterations} times, which is {np.count_nonzero(all_random_winner == 6)/40000 * 100}%')
print('-----------------------------------------------------------')