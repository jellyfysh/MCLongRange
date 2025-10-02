# MCLongRange - repository accompanying the manuscript "Markov-chain sampling for long-range systems without
# evaluating the energy" by Gabriele Tartero & Werner Krauth - https://github.com/jellyfysh/MCLongRange
# Copyright (C) 2024 The JeLLyFysh organization
#
# MCLongRange is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
# version (see the LICENSE file).
#
# This program implements the factorized Metropolis algorithm for a two dimensional Lennard-Jones system.
# It corresponds to Algorithm 2 (factorized-metropolis) in the manuscript.
#
# To produce output files, uncomment the last lines of the script.
#
import random
import math
import os.path
from functions import per_dist, u_lj

N = 100
density = 0.05
L = math.sqrt(N / density)
delta = 1.0

# initial configuration
filename = 'InitialConfs/N' + str(N) + '_rho' + str(density) + '.data'
if os.path.isfile(filename):
    with open(filename, 'r') as file:
        conf = [list(c) for c in eval(file.readline())]
else:
    conf = [[random.uniform(0.0, L), random.uniform(0.0, L)] for i in range(N)]

u_test = sum([u_lj(per_dist(conf[i], conf[j], L)) for i in range(N) for j in range(i, N) if j != i])

# sampling
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
    delta_u = 0.0
    for c in conf:
        if c == part:
            continue
        u_evals += 2
        delta_u_pair = u_lj(per_dist(new_part, c, L)) - u_lj(per_dist(part, c, L))
        metr_fil = math.exp(-delta_u_pair) if delta_u_pair > 0.0 else 1.0
        delta_u += delta_u_pair
        if random.uniform(0.0, 1.0) > metr_fil:
            break
    else:
        acc += 1
        distance += math.sqrt(del_x ** 2 + del_y ** 2)
        u_test += delta_u
        part[:] = new_part
    pair_part = random.sample(conf, 2)
    pair_correlations.append(per_dist(pair_part[0], pair_part[1], L))

u_final = sum([u_lj(per_dist(conf[i], conf[j], L)) for i in range(N) for j in range(i, N) if j != i])

print("Sanity check: |u_test - u_final| = " + str(abs(u_test - u_final)))
print("Acceptance rate: " + str(acc / n_samples))
print("Evaluations/distance: " + str(u_evals / distance))

with open(filename, "w") as file:
    file.write(str(conf))

# to produce a file containing scaling data for this algorithm, uncomment the following two lines and run the program
# several times with different values of N
# with open('ScalingData/FactorizedMetropolisScaling_rho' + str(density) + '.data', "a") as file:
#     file.write(str([N, u_evals / distance]) + '\n')

# to produce a file containing pair-correlation data for this algorithm, uncomment the following two lines
# with open('PairCorrelationData/FactorizedMetropolis_' + str(N) + '_rho' + str(density) + '_Correlation.data', "w") as file:
#     file.write(str(pair_correlations))
