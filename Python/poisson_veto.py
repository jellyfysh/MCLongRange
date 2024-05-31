import random
import math

N = 5
q_probs = [random.random() for i in range(N)]
intensities = [-math.log(1 - q) for q in q_probs]
print('Rejection probabilities: ' + str(q_probs))

n_trials = 10 ** 7

# poisson-veto
histo_naive = {i: 0 for i in range(N)}
for n in range(n_trials):
    S_veto = []
    for i in range(N):
        lambda_i = intensities[i]
        t = -math.log(random.random()) / lambda_i
        if t < 1.0:
            S_veto.append(i)
            histo_naive[i] += 1
naive_probs = [histo_naive[i] / n_trials for i in range(N)]
print('Experimental probabilities (naive): ' + str(naive_probs))

# constructing Walker tables
lambda_tot = sum(intensities)
lambda_probs = [lambda_i / lambda_tot for lambda_i in intensities]
prob_table = [N * lambda_prob for lambda_prob in lambda_probs]
alias_table = [k for k in range(N)]
over = []
under = []
for s in range(N):
    if prob_table[s] > 1.0:
        over.append(s)
    elif prob_table[s] < 1.0:
        under.append(s)
while over and under:
    i = random.choice(over)
    j = random.choice(under)
    alias_table[j] = i
    under.remove(j)
    prob_table[i] += prob_table[j] - 1.0
    if prob_table[i] <= 1.0:
        over.remove(i)
    if prob_table[i] < 1.0:
        under.append(i)

# poisson-veto(patch)
histo_patch = {i: 0 for i in range(N)}
for n in range(n_trials):
    S_veto = set()
    t = 0.0
    while True:
        t += (-math.log(random.random()) / lambda_tot)
        if t > 1.0:
            break
        j = random.choice(range(N))
        if random.random() < prob_table[j]:
            S_veto.add(j)
        else:
            S_veto.add(alias_table[j])
    for j in S_veto:
        histo_patch[j] += 1
patch_probs = [histo_patch[i] / n_trials for i in range(N)]
print('Experimental probabilities (patch): ' + str(patch_probs))
