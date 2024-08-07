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

pyplot.rcParams["text.usetex"] = True

pyplot.xlabel(r"$|\Delta r|$", fontsize=20)
pyplot.ylabel(r"$\pi(|\Delta r|)$", fontsize=20)
for file in glob.glob("PairCorrelationData/*.data"):
    label = file.split("data")[0][20:-12]
    data = eval(open(file).readline())
    bins_x = np.linspace(min(data), max(data), 100)
    pyplot.hist(data, bins_x, label=label, density=True, cumulative=False, histtype="step")
pyplot.legend(fontsize=18)
pyplot.show()