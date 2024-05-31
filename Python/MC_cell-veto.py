import random
import math
import os.path
from functions import u_lj, per_dist

N = 100
density = 0.05
L = math.sqrt(N / density)
delta = 1.0

# generating the cell rates
n = 20
n_cells = n ** 2
cell_size = L / n
neighbor_indices = [(i % n) + n * (j % n) for i in [-1, 0, 1] for j in [-1, 0, 1]]
n_trials = 10 ** 3
cell_rates = []
for k in range(n):
    print("Row number " + str(k))
    for j in range(n):
        if j + n * k in neighbor_indices:
            cell_rates.append(0.0)
        else:
            x_active = [[random.uniform(0, cell_size), random.uniform(0, cell_size)] for n in range(n_trials)]
            new_x = [[(x[0] + random.uniform(-delta, delta)) % L, (x[1] + random.uniform(-delta, delta)) % L]
                     for x in x_active]
            x_target = [[random.uniform(j * cell_size, (j + 1) * cell_size),
                        random.uniform(k * cell_size, (k + 1) * cell_size)] for m in range(n_trials)]
            delta_u = max([u_lj(per_dist(new_x[i], x_target[i], L)) - u_lj(per_dist(x_active[i], x_target[i], L))
                           for i in range(n_trials)])
            if 1.0 - math.exp(-5.0 * delta_u) == 1.0:
                cell_rates.append(0.0)
                neighbor_indices.append(j + n * k)
            else:
                cell_rates.append(1.0 - math.exp(-5.0 * delta_u))

# generating Walker tables
intensities = [-math.log(1.0 - q) for q in cell_rates]
P = sum(intensities)
p_probs = [p / P for p in intensities]
n_p = len(p_probs)
prob_table = [n_p * p for p in p_probs]
alias_table = [k for k in range(n_p)]
over = []
under = []
for s in range(n_p):
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
        particle_cells[r * j].append([x, y])
    surplus_cells = []

# sampling valid configurations
n_samples = 10 ** 2
u_evals = 0
distance = 0.0
pair_correlations = []
for sample in range(n_samples):
    print(sample)
    part = random.choice(conf)
    [del_x, del_y] = [random.uniform(-delta, delta), random.uniform(-delta, delta)]
    new_part = [(part[0] + del_x) % L, (part[1] + del_y) % L]
    active_cell = int(part[0] / cell_size) + n * int(part[1] / cell_size)
    new_cell = int(new_part[0] / cell_size) + n * int(new_part[1] / cell_size)
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
    # treating surplus and neighbor particles
    for s in surplus_particles + neighbor_particles:
        u_evals += 2
        delta_u = u_lj(per_dist(new_part, s, L)) - u_lj(per_dist(part, s, L))
        metr_fil = math.exp(-delta_u) if delta_u > 0.0 else 1.0
        if random.random() > metr_fil:
            break
    else:
        # identifying the set of target cells
        powerset = set()
        t = 0.0
        while True:
            t -= math.log(random.uniform(0.0, 1.0)) / P
            if t > 1.0:
                break
            i = random.choice(range(n_p))
            if random.random() < prob_table[i]:
                powerset.add(i)
            else:
                powerset.add(alias_table[i])
        target_cells = {(active_cell % n + j % n) % n + n * ((int(active_cell / n) + int(j / n)) % n): cell_rates[j]
                        for j in powerset}
        # identifying target particles
        target_particles = []
        for t_cell in target_cells.keys():
            if particle_cells[t_cell]:
                target_particles.append(particle_cells[t_cell][0])
        # deciding whether the former rejections were real or not
        for t in target_particles:
            u_evals += 2
            t_cell = int(t[0] / cell_size) + n * int(t[1] / cell_size)
            delta_u = u_lj(per_dist(new_part, t, L)) - u_lj(per_dist(part, t, L))
            fil = math.exp(-delta_u) if delta_u > 0.0 else 1.0
            if random.random() < (1.0 - fil) / target_cells[t_cell]:
                break
        else:
            # if there are no real rejections, update the active particle's position and associate it to its cell
            distance += math.sqrt(del_x ** 2 + del_y ** 2)
            particle_cells[active_cell].remove(part)
            particle_cells[new_cell].append(new_part)
            part[:] = new_part
            # updating the list of surplus particles
            if new_cell != active_cell:
                if len(particle_cells[active_cell]) == 1:
                    surplus_cells.remove(active_cell)
                if len(particle_cells[new_cell]) == 2:
                    surplus_cells.append(new_cell)
    pair_part = random.sample(conf, 2)
    pair_correlations.append(per_dist(pair_part[0], pair_part[1], L))

print("Evaluations/length: " + str(u_evals / distance))

with open(filename, "w") as file:
    file.write(str(conf))
# with open('ScalingData/MCCellVetoScaling.data', "a") as file:
#     file.write(str([N, u_evals / distance]) + '\n')
# with open('PairCorrelationData/MCCellVeto' + str(N) + 'Correlation.data', "w") as file:
#     file.write(str(pair_correlations))
