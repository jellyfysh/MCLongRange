# MCLongRange - repository accompanying the manuscript "Markov-chain sampling for long-range systems without
# evaluating the energy" by Gabriele Tartero & Werner Krauth - https://github.com/jellyfysh/MCLongRange
# Copyright (C) 2024 The JeLLyFysh organization
#
# MCLongRange is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
# version (see the LICENSE file).
#
# This program implements the multi-time-step Metropolis algorithm for a two dimensional Lennard-Jones system.
# It corresponds to Algorithm 1 (multi-step-metropolis) in the manuscript.
#
# To produce output files, uncomment the last lines of the script.
#
import math
import random
import os.path
from functions import u_lj, per_dist, sigma

N = 100
density = 0.05
L = math.sqrt(N / density)
delta = 1.0
r_c = 1.3 * sigma
n_short = 4

# initial configuration
filename = 'InitialConfs/N' + str(N) + '_rho' + str(density) + '.data'
if os.path.isfile(filename):
    with open(filename, 'r') as file:
        conf = [list(c) for c in eval(file.readline())]
else:
    conf = [[random.uniform(0.0, L), random.uniform(0.0, L)] for i in range(N)]

u_short_test = sum([u_lj(per_dist(conf[i], conf[j], L)) for i in range(N) for j in range(i, N)
                    if 0 < per_dist(conf[i], conf[j], L) < r_c])
u_long_test = sum([u_lj(per_dist(conf[i], conf[j], L)) for i in range(N) for j in range(i, N)
                   if per_dist(conf[i], conf[j], L) >= r_c])

# sampling
n_samples = 10 ** 2
u_evals = 0
distance = 0.0
short_acc = 0
long_acc = 0
pair_correlations = []
for sample in range(n_samples):
    print(sample)
    Y_old = {}  # keeping track of the particles that are displaced with short-range moves
    # short-range steps
    delta_u_short = 0.0
    dist_short = 0.0
    for n_s in range(n_short):
        i = random.randint(0, N - 1)
        part = conf[i]
        [del_x, del_y] = [random.uniform(-delta, delta), random.uniform(-delta, delta)]
        new_part = ((part[0] + del_x) % L, (part[1] + del_y) % L)
        u_short_old = 0.0
        u_short_new = 0.0
        for j in list(range(i)) + list(range(i + 1, N)):
            if per_dist(conf[i], conf[j], L) < r_c:
                u_short_old += u_lj(per_dist(conf[i], conf[j], L))
                u_evals += 1
            if per_dist(new_part, conf[j], L) < r_c:
                u_evals += 1
                u_short_new += u_lj(per_dist(new_part, conf[j], L))
        metr_fil_short = math.exp(-(u_short_new - u_short_old)) if u_short_new - u_short_old > 0.0 else 1.0
        if random.uniform(0.0, 1.0) < metr_fil_short:
            short_acc += 1
            dist_short += math.sqrt(del_x ** 2 + del_y ** 2)
            delta_u_short += (u_short_new - u_short_old)
            conf[i] = new_part
            if i not in Y_old:
                Y_old[i] = part  # part has not been moved yet
    # long-range decision
    delta_u_long = 0.0
    for i in Y_old:
        for j in list(range(i)) + list(range(i + 1, N)):
            alpha = 0.5 if j in Y_old else 1.0
            if per_dist(conf[i], conf[j], L) >= r_c:
                delta_u_long += alpha * u_lj(per_dist(conf[i], conf[j], L))
                u_evals += 1
            target = Y_old[j] if j in Y_old else conf[j]
            if per_dist(Y_old[i], target, L) >= r_c:
                delta_u_long -= alpha * u_lj(per_dist(Y_old[i], target, L))
                u_evals += 1
    metr_fil_long = math.exp(-delta_u_long) if delta_u_long > 0.0 else 1.0
    if random.uniform(0.0, 1.0) > metr_fil_long:
        for i in Y_old:
            conf[i] = Y_old[i]
    else:
        long_acc += 1
        distance += dist_short
        u_short_test += delta_u_short
        u_long_test += delta_u_long
    pair_part = random.sample(conf, 2)
    pair_correlations.append(per_dist(pair_part[0], pair_part[1], L))

u_short_final = sum([u_lj(per_dist(conf[i], conf[j], L)) for i in range(N) for j in range(i, N)
                     if 0 < per_dist(conf[i], conf[j], L) < r_c])
u_long_final = sum([u_lj(per_dist(conf[i], conf[j], L)) for i in range(N) for j in range(i, N)
                    if per_dist(conf[i], conf[j], L) >= r_c])

print("Sanity check: |u_short_test - u_short_final| = " + str(abs(u_short_test - u_short_final)) +
      ", |u_long_test - u_long_final| = " + str(abs(u_long_test - u_long_final)))
print("Short-range acceptance rate: " + str(short_acc / (n_samples * n_short)))
print("Long-range acceptance rate: " + str(long_acc / n_samples))
print("Evaluations/distance: " + str(u_evals / distance))

with open(filename, "w") as file:
    file.write(str(conf))

# to produce a file containing scaling data for this algorithm, uncomment the following two lines and run the program
# several times with different values of N
# with open('ScalingData/MultiStep' + str(n_short) + 'Scaling.data', "a") as file:
#     file.write(str([N, u_evals / distance]) + '\n')

# to produce a file containing pair-correlation data for this algorithm, uncomment the following two lines
# with open('PairCorrelationData/MultiStep' + str(N) + '_' + str(n_short) + 'Correlation.data', "w") as file:
#     file.write(str(pair_correlations))
