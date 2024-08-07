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

pyplot.xlabel("Number of particles", fontsize=20)
pyplot.ylabel("Energy evaluations per displacement", fontsize=20)
pyplot.xscale("log")
pyplot.yscale("log")
markers = ["d", "o", "v", "^", "X"]
for file in glob.glob('ScalingData/*.data'):
    j = glob.glob('ScalingData/*.data').index(file)
    label = file.split("data")[0][12:-8]
    data = [eval(line.rstrip()) for line in open(file, 'r')]
    data.sort()
    n_particles = [d[0] for d in data]
    scaling_data = [d[1] for d in data]
    pyplot.plot(n_particles, scaling_data, label=label, linestyle='-', marker=markers[j], markersize=10)
pyplot.legend(fontsize=18)
pyplot.show()