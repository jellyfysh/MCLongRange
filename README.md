## MCLongRange
This repository accompanies the article 
[Markov-chain sampling for long-range systems without evaluating the energy](https://doi.org/10.48550/arXiv.2406.00184),
by Gabriele Tartero and Werner Krauth. It contains all the supplementary 
material discussed in Appendix A.

### Python programs
The [Python](https://github.com/jellyfysh/MCLongRange/tree/master/Python) directory
contains five simulations of a system of $N$ particles interacting via a Lennard-Jones 
potential in a periodic box, with $\beta = 1$: 
* [metropolis.py](https://github.com/jellyfysh/MCLongRange/blob/master/Python/metropolis.py),
  the standard (non-factorized) Metropolis algorithm (discussed in Section II.A);
* [multi-step_metropolis.py](https://github.com/jellyfysh/MCLongRange/blob/master/Python/multi-step_metropolis.py),
  the multi-time-step Metropolis (Algorithm 1 (multi-step-metropolis) in the manuscript);
* [factorized_metropolis.py](https://github.com/jellyfysh/MCLongRange/blob/master/Python/factorized_metropolis.py),
  the factorized Metropolis algorithm (Algorithm 2 (factorized-metropolis) in the manuscript);
* [MC_cell-veto.py](https://github.com/jellyfysh/MCLongRange/blob/master/Python/MC_cell-veto.py),
  the reversible version of the cell-veto algorithm (Algorithm 4 (cell-veto(patch)) in the manuscript, with the 
  set $\mathcal{S}_{\text{veto}}$ sampled using Algorithm 6 (poisson-veto(patch)));
* [EC_cell-veto.py](https://github.com/jellyfysh/MCLongRange/blob/master/Python/EC_cell-veto.py),
  the non-reversible version of the cell-veto algorithm (Algorithm 4 (cell-veto(patch)) in the manuscript, with the 
  set $\mathcal{S}_{\text{veto}}$ sampled using Walker's algorithm).

Moreover, the script [posson_veto.py](https://github.com/jellyfysh/MCLongRange/blob/master/Python/poisson_veto.py) implements Algorithms 5 (poisson-veto) and 6 (poisson-veto(patch)) in a pedagogical way, and compares their outcomes.

All these algorithms rely on [functions.py](https://github.com/jellyfysh/MCLongRange/blob/master/Python/functions.py),
which contains some basic functions and parameters. For given values of $N$ and $\rho$ (the system density),
the final configuration is stored in 
[InitialConf/](https://github.com/jellyfysh/MCLongRange/tree/master/Python/InitialConfs), and then used as the starting configuration of the following run. In this way, it is sufficient 
to reach the equilibrium
only once. The consistency of the implementations is checked by comparing the histograms for 
the pair-correlation function. The experimental data are stored in 
[PairCorrelationData/](https://github.com/jellyfysh/MCLongRange/tree/master/Python/PairCorrelationData)
and the histograms are produced with the script 
[pair_correlation.py](https://github.com/jellyfysh/MCLongRange/blob/master/Python/pair_correlation.py). 
Some programs also contain internal sanity checks, based on the computation of energy variations or acceptance rates. 
As for the scaling, the data are stored in [ScalingData/](https://github.com/jellyfysh/MCLongRange/tree/master/Python/ScalingData) and plotted with [scalings.py](https://github.com/jellyfysh/MCLongRange/blob/master/Python/scalings.py).
This program can be used to reproduce Figure 10 of the manuscript.

All programs can be executed with any Python3 implementation 
(e.g., standard [CPython](https://www.python.org/) or 
[PyPy3](https://www.pypy.org/)). The scripts
[pair_correlation.py](https://github.com/jellyfysh/MCLongRange/blob/master/Python/pair_correlation.py)
and [scalings.py](https://github.com/jellyfysh/MCLongRange/blob/master/Python/scalings.py)
rely both on [NumPy](https://numpy.org/) and 
[Matplotlib](https://matplotlib.org/) to produce the plots, while
all other programs do not have any further requirements.

### Authors
The authors of this project are:
* Gabriele Tartero 
([gabriele.tartero@phys.ens.fr](mailto:gabriele.tartero@phys.ens.fr));
* Werner Krauth ([werner.krauth@ens.fr](mailto:werner.krauth@ens.fr)).

For any question about the MCLongRange software package, or the related
paper, please raise an issue here on GitHub or contact us via e-mail.

### License
This project is licensed under the GNU General Public License, 
version 3 (see the 
[LICENSE](https://github.com/jellyfysh/MCMCNutshell/blob/master/LICENSE) 
file).


