# MCLongRange - repository accompanying the manuscript "Markov-chain sampling for long-range systems without
# evaluating the energy" by Gabriele Tartero & Werner Krauth - https://github.com/jellyfysh/MCLongRange
# Copyright (C) 2024 The JeLLyFysh organization
#
# MCLongRange is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
# version (see the LICENSE file).
#
# This program implements the non-reversible cell-veto algorithm for a two dimensional Lennard-Jones system.
# It corresponds to the event-chain version of Algorithm 4 (cell-veto(patch)) in the manuscript, with the
# set S_veto sampled using Walker's algorithm (illustrated in Figure 8).
# Event-chain Monte Carlo algorithms are not explicitly discussed in the manuscript. For a full explanation, see
# references [32], [36], [38].
#
# To produce output files, uncomment the last lines of the script.
#
import random
import math
import os.path
from functions import du_lj, per_dist, displacement_advanced

N = 100
density = 0.05
L = math.sqrt(N / density)
directions = [[1.0, 0.0], [0.0, 1.0]]
chain_length = 2.0

# generating the differential cell rates
n = 20
n_cells = n ** 2
cell_size = L / n
neighbor_indices = [(i % n) + n * (j % n) for i in [-1, 0, 1] for j in [-1, 0, 1]]
n_trials = 10 ** 3
diff_rates = []
for k in range(n):
    print("Row number " + str(k))
    for j in range(n):
        if j + n * k in neighbor_indices:
            diff_rates.append(0.0)
        else:
            x_active = [[random.uniform(0, cell_size), random.uniform(0, cell_size)] for n in range(n_trials)]
            x_target = [[random.uniform(j * cell_size, (j + 1) * cell_size),
                         random.uniform(k * cell_size, (k + 1) * cell_size)] for m in range(n_trials)]
            distances = [per_dist(x, y, L) for x, y in zip(x_active, x_target)]
            if 5.0 * max([abs(du_lj(r)) for r in distances]) > 0.5:
                diff_rates.append(0.0)
                neighbor_indices.append(j + n * k)
            else:
                diff_rates.append(5.0 * max([abs(du_lj(r)) for r in distances]))

# generating Walker tables
Q = sum(diff_rates)
q_probs = [q / Q for q in diff_rates]
n_q = len(q_probs)
prob_table = [n_q * q for q in q_probs]
alias_table = [k for k in range(n_q)]
over = []
under = []
for s in range(n_q):
    if prob_table[s] > 1.0:
        over.append(s)
    elif prob_table[s] < 1.0:
        under.append(s)
count = 0
while over and under:
    print(count)
    count += 1
    i = random.choice(over)
    j = random.choice(under)
    alias_table[j] = i
    under.remove(j)
    prob_table[i] += prob_table[j] - 1.0
    if prob_table[i] <= 1.0:
        over.remove(i)
    if prob_table[i] < 1.0:
        under.append(i)

# initial configuration
filename = 'InitialConfs/N' + str(N) + '_rho' + str(density) + '.data'
if os.path.isfile(filename):
    with open(filename, 'r') as file:
        conf = [list(c) for c in eval(file.readline())]
        particle_cells = {cell: [] for cell in range(n_cells)}
        surplus_cells = []
        for c in conf:
            cell = int(c[0] / cell_size) + n * int(c[1] / cell_size)
            particle_cells[cell].append(c)
            if len(particle_cells[cell]) > 1 and cell not in surplus_cells:
                surplus_cells.append(cell)
else:
    conf = []
    particle_cells = {cell: [] for cell in range(n_cells)}
    r = int(n_cells / N)
    for j in range(N):
        [n_x, n_y] = [(r * j) % n, int((r * j) / n)]
        x = random.uniform(n_x * cell_size, (n_x + 1) * cell_size)
        y = random.uniform(n_y * cell_size, (n_y + 1) * cell_size)
        conf.append([x, y])
        particle_cells[(r * j)].append([x, y])
    surplus_cells = []

# sampling
n_samples = 10 ** 2
u_evals = 0
distance = 0.0
pair_correlations = []
for sample in range(n_samples):
    print(sample)
    sigma = random.choice(directions)
    part = random.choice(conf)
    current_length = 0.0
    while True:
        active_cell = int(part[0] / cell_size) + n * int(part[1] / cell_size)
        # identifying neighbor cells and neighbor particles
        neighbor_cells = [(active_cell % n + k % n) % n + n * ((int(active_cell / n) + int(k / n)) % n)
                          for k in neighbor_indices]
        neighbor_particles = []
        for n_cell in neighbor_cells:
            possible_neighbors = [particle for particle in particle_cells[n_cell] if particle != part]
            if possible_neighbors:
                neighbor_particles.extend(possible_neighbors[0:])
        # identifying surplus particles
        surplus_particles = []
        for s_cell in surplus_cells:
            if s_cell not in neighbor_cells:
                surplus_particles.extend(particle_cells[s_cell][1:])
        extra_particles = neighbor_particles + surplus_particles
        # treating neighbor and surplus particles
        extra_u_hats = [-math.log(random.random()) for i in range(len(extra_particles))]
        extra_targets = [[y, L - x] for [x, y] in extra_particles] if sigma == [0.0, 1.0] else extra_particles
        active = part if sigma == [1.0, 0.0] else [part[1], L - part[0]]
        displacements = [displacement_advanced(active, c, u_hat, L) for c, u_hat in zip(extra_targets, extra_u_hats)]
        u_evals += len(extra_targets)
        # treating target particles
        FakeRejection = False
        if len(extra_particles) < len(conf) - 1:  # making sure that target particles actually exist
            x = active[0]
            act_cell = int(active[0] / cell_size) + n * int(active[1] / cell_size)
            act_coord = [act_cell % n, int(act_cell / n)]
            delta_s_max = per_dist([(act_coord[0] + 1) * cell_size], [x], L)
            delta_s_max = cell_size if delta_s_max == 0.0 else delta_s_max
            delta_s = -math.log(random.random()) / Q
            u_evals += 1
            if delta_s > delta_s_max:
                displacements.append(delta_s_max)
                FakeRejection = True
            else:
                displacements.append(delta_s)
        # finding the next displacement
        delta_s = min(displacements)
        # if the total length is reached, the chain ends
        if current_length + delta_s > chain_length:
            distance += (chain_length - current_length)
            new_part = [(part[j] + sigma[j] * (chain_length - current_length)) % L for j in range(2)]
            new_cell = int(new_part[0] / cell_size) + n * int(new_part[1] / cell_size)
            particle_cells[active_cell].remove(part)
            particle_cells[new_cell].append(new_part)
            part[:] = new_part
            # updating the list of surplus particles
            if new_cell != active_cell:
                if len(particle_cells[active_cell]) == 1:
                    surplus_cells.remove(active_cell)
                if len(particle_cells[new_cell]) == 2:
                    surplus_cells.append(new_cell)
            break
        # if the total length has not been reached yet, the chain continues
        distance += delta_s
        current_length += delta_s
        new_part = [(part[j] + sigma[j] * delta_s) % L for j in range(2)]
        new_cell = int(new_part[0] / cell_size) + n * int(new_part[1] / cell_size)
        particle_cells[active_cell].remove(part)
        particle_cells[new_cell].append(new_part)
        part[:] = new_part
        # updating the list of surplus particles
        if new_cell != active_cell:
            if len(particle_cells[active_cell]) == 1:
                surplus_cells.remove(active_cell)
            if len(particle_cells[new_cell]) == 2:
                surplus_cells.append(new_cell)
        # selecting the new active particle
        k = displacements.index(delta_s)
        if k < len(extra_particles):
            j = conf.index(extra_particles[k])
            part = conf[j]
        elif k == len(extra_particles) and not FakeRejection:
            t_index = random.choice(range(n_q))
            if random.random() > prob_table[t_index]:
                t_index = alias_table[t_index]
            t_cell = (active_cell % n + t_index % n) % n + n * ((int(active_cell / n) + int(t_index / n)) % n)
            if particle_cells[t_cell]:
                target_particle = particle_cells[t_cell][0]
                if random.random() < du_lj(per_dist(part, target_particle, L)) / diff_rates[t_index]:
                    j = conf.index(target_particle)
                    part = conf[j]
    pair_part = random.sample(conf, 2)
    pair_correlations.append(per_dist(pair_part[0], pair_part[1], L))

print("Evaluations/distance: " + str(u_evals / distance))

with open(filename, "w") as file:
    file.write(str(conf))

# to produce a file containing scaling data for this algorithm, uncomment the following two lines and run the program
# several times with different values of N
# with open('ScalingData/ECCellVetoScaling.data', "a") as file:
#     file.write(str([N, u_evals / distance]) + '\n')

# to produce a file containing pair-correlation data for this algorithm, uncomment the following two lines
# with open('PairCorrelationData/ECCellVeto' + str(N) + 'Correlation.data', "w") as file:
#     file.write(str(pair_correlations))
