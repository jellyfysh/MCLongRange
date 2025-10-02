# MCLongRange - repository accompanying the manuscript "Markov-chain sampling for long-range systems without
# evaluating the energy" by Gabriele Tartero & Werner Krauth - https://github.com/jellyfysh/MCLongRange
# Copyright (C) 2024 The JeLLyFysh organization
#
# MCLongRange is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
# version (see the LICENSE file).
#
# This program produces an histogram of pair-correlation data, stored in the folder PairCorrelationData/, for
# different Markov chain algorithms.
# To draw a cumulative histogram, set the "cumulative" parameter to "True" on line 24
#
import numpy as np
from matplotlib import pyplot
import glob
import scienceplots

pyplot.style.use("science")

N = 100
density = 0.05

title = r"$N = $" + str(N) + r"$, \, \rho = $" + str(density)

fig, ax = pyplot.subplots(figsize=(7, 5))
ax.set_xlabel(r"$|\Delta r|$", fontsize=20)
ax.set_ylabel(r"$\pi(|\Delta r|)$", fontsize=20)
for file in glob.glob("PairCorrelationData/*_N" + str(N) + "_rho" + str(density) + "*_Correlation.data"):
    label = file.split("_")[0][20:]
    data = eval(open(file).readline())
    bins_x = np.linspace(min(data), max(data), 100)
    ax.hist(data, bins_x, label=label, density=True, cumulative=False, histtype="step")
ax.legend(fontsize=16)
ax.set_title(title, fontsize=18)
pyplot.show()
