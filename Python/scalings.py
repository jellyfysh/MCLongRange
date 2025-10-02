# MCLongRange - repository accompanying the manuscript "Markov-chain sampling for long-range systems without
# evaluating the energy" by Gabriele Tartero & Werner Krauth - https://github.com/jellyfysh/MCLongRange
# Copyright (C) 2024 The JeLLyFysh organization
#
# MCLongRange is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
# version (see the LICENSE file).
#
# This program plots the scaling data, stored in the folder ScalingData/, for different Markov chain algorithms.
# It can be used to reproduce Figure 10 of the manuscript.
#
import glob
from matplotlib import pyplot
import scienceplots

pyplot.style.use("science")

density = 0.05

title = r"$\rho = $" + str(density)

fig, ax = pyplot.subplots(figsize=(7, 5))
ax.set_xlabel("Number of particles", fontsize=20)
ax.set_ylabel("Energy evaluations per displacement", fontsize=20)
ax.set_xscale("log")
ax.set_yscale("log")
markers = ["d", "o", "v", "^", "X"]
for file in glob.glob('ScalingData/*.data'):
    j = glob.glob('ScalingData/*.data').index(file)
    label = file.split("_")[0][12:-7]
    data = [eval(line.rstrip()) for line in open(file, 'r')]
    data.sort()
    n_particles = [d[0] for d in data]
    scaling_data = [d[1] for d in data]
    ax.plot(n_particles, scaling_data, label=label, linestyle='-', marker=markers[j], markersize=10)
pyplot.legend(fontsize=16)
ax.set_title(title, fontsize=18)
pyplot.show()
