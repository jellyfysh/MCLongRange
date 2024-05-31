import random
import math
import os.path
from functions import per_dist, u_lj

N = 100
density = 0.05
L = math.sqrt(N / density)
delta = 1.0

filename = 'InitialConfs/N' + str(N) + '_rho' + str(density) + '.data'
if os.path.isfile(filename):
    with open(filename, 'r') as file:
        conf = [list(c) for c in eval(file.readline())]
else:
    conf = [[random.uniform(0.0, L), random.uniform(0.0, L)] for i in range(N)]

# u_test = sum([u_lj(per_dist(conf[i], conf[j], L)) for i in range(N) for j in range(i, N) if j != i])

n_samples = 10 ** 2
u_evals = 0
distance = 0.0
acc = 0
pair_correlations = []
for sample in range(n_samples):
    print(sample)
    part = random.choice(conf)
    [del_x, del_y] = [random.uniform(-delta, delta), random.uniform(-delta, delta)]
    new_part = [(part[0] + del_x) % L, (part[1] + del_y) % L]
    delta_u = sum([u_lj(per_dist(new_part, c, L)) - u_lj(per_dist(part, c, L)) for c in conf if c != part])
    u_evals += (2 * (N - 1))
    if random.uniform(0.0, 1.0) < math.exp(min(0.0, -delta_u)):
        acc += 1
        distance += math.sqrt(del_x ** 2 + del_y ** 2)
        # u_test += delta_u
        part[:] = new_part
    pair_part = random.sample(conf, 2)
    pair_correlations.append(per_dist(pair_part[0], pair_part[1], L))

# u_final = sum([u_lj(per_dist(conf[i], conf[j], L)) for i in range(N) for j in range(i, N) if j != i])

# print("Sanity check: |u_test - u_final| = " + str(abs(u_test - u_final)))
print("Acceptance rate: " + str(acc / n_samples))
print("Evaluations/length: " + str(u_evals / distance))

with open(filename, "w") as file:
    file.write(str(conf))
# with open('ScalingData/MetropolisScaling.data', "a") as file:
#     file.write(str([N, u_evals / distance]) + '\n')
# with open('PairCorrelationData/Metropolis' + str(N) + 'Correlation.data', "w") as file:
#     file.write(str(pair_correlations))
