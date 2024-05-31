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
    n_particles = [d[0] for d in data]
    scaling_data = [d[1] for d in data]
    pyplot.plot(n_particles, scaling_data, label=label, linestyle='-', marker=markers[j], markersize=10)
pyplot.legend(fontsize=18)
pyplot.show()
