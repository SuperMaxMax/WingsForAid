import numpy as np

### scores ###
#[mass, L/D, complexity, C&S, operation]
non_brace = np.array([3, 3, 3, 3, 3])
light_brace = np.array([5, 2, 4, 3, 2])
lift_gen_brace = np.array([4, 5, 1, 2, 1])

original_weights = np.array([0.275, 0.2, 0.275, 0.1, 0.15])
iterations = 10000

#define storage of all weights
all_random_winner = np.zeros((5,iterations))

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

        score_non_brace = np.dot(np.transpose(non_brace), new_weight)
        score_light_brace = np.dot(np.transpose(light_brace), new_weight) 
        score_lift_gen_brace = np.dot(np.transpose(lift_gen_brace), new_weight) 

        if score_non_brace > score_light_brace and score_non_brace > score_light_brace:
            winner = 1
        elif score_light_brace > score_lift_gen_brace and score_light_brace > score_non_brace: 
            winner = 2 
        elif score_lift_gen_brace > score_non_brace and score_lift_gen_brace > score_light_brace:
            winner = 3
        else:
            winner = 4 #no clear winner

        all_random_winner[i][j] = winner

print('By changing the mass weight the winners are distributed as follows:')
print('non brace wins', np.count_nonzero(all_random_winner[0,:] == 1), f'out of {iterations} times')
print('light brace wins', np.count_nonzero(all_random_winner[0,:] == 2), f'out of {iterations} times')
print('lift generator brace wins', np.count_nonzero(all_random_winner[0,:] == 3), f'out of {iterations} times')
print('there is no clear winner', np.count_nonzero(all_random_winner[0,:] == 4), f'out of {iterations} times')
print('--------------------------------------------------')

print('By changing the L/D weight the winners are distributed as follows:')
print('non brace wins', np.count_nonzero(all_random_winner[1,:] == 1), f'out of {iterations} times')
print('light brace wins', np.count_nonzero(all_random_winner[1,:] == 2), f'out of {iterations} times')
print('lift generator brace wins', np.count_nonzero(all_random_winner[1,:] == 3), f'out of {iterations} times')
print('there is no clear winner', np.count_nonzero(all_random_winner[1,:] == 4), f'out of {iterations} times')
print('--------------------------------------------------')

print('By changing the complexity weight the winners are distributed as follows:')
print('non brace wins', np.count_nonzero(all_random_winner[2,:] == 1), f'out of {iterations} times')
print('light brace wins', np.count_nonzero(all_random_winner[2,:] == 2), f'out of {iterations} times')
print('lift generator brace wins', np.count_nonzero(all_random_winner[2,:] == 3), f'out of {iterations} times')
print('there is no clear winner', np.count_nonzero(all_random_winner[2,:] == 4), f'out of {iterations} times')
print('--------------------------------------------------')

print('By changing the C&S weight the winners are distributed as follows:')
print('non brace wins', np.count_nonzero(all_random_winner[3,:] == 1), f'out of {iterations} times')
print('light brace wins', np.count_nonzero(all_random_winner[3,:] == 2), f'out of {iterations} times')
print('lift generator brace wins', np.count_nonzero(all_random_winner[3,:] == 3), f'out of {iterations} times')
print('there is no clear winner', np.count_nonzero(all_random_winner[3,:] == 4), f'out of {iterations} times')
print('--------------------------------------------------')

print('By changing the operation weight the winners are distributed as follows:')
print('non brace wins', np.count_nonzero(all_random_winner[4,:] == 1), f'out of {iterations} times')
print('light brace wins', np.count_nonzero(all_random_winner[4,:] == 2), f'out of {iterations} times')
print('lift generator brace wins', np.count_nonzero(all_random_winner[4,:] == 3), f'out of {iterations} times')
print('there is no clear winner', np.count_nonzero(all_random_winner[4,:] == 4), f'out of {iterations} times')
print('--------------------------------------------------')

print('the overall distribution is:')
print(f'non brace wins {np.count_nonzero(all_random_winner == 1)} out of {5 * iterations} times, which is {np.count_nonzero(all_random_winner == 1)/50000 * 100}%')
print(f'light brace wins {np.count_nonzero(all_random_winner == 2)} out of {5 * iterations} times, which is {np.count_nonzero(all_random_winner == 2)/50000 * 100}%')
print(f'lift generator brace wins {np.count_nonzero(all_random_winner == 3)} out of {5 * iterations} times, which is {np.count_nonzero(all_random_winner == 3)/50000 * 100}%')
print(f'there is no clear winner {np.count_nonzero(all_random_winner == 4)} out of {5 * iterations} times, which is {np.count_nonzero(all_random_winner == 4)/50000 * 100}%')
print('--------------------------------------------------')