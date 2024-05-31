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
    pyplot.hist(data, bins_x, label=label, density=True, cumulative=True, histtype="step")
pyplot.legend(fontsize=18)
pyplot.show()
