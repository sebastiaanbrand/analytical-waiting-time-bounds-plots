# Improved analytical bounds on delivery times of long-distance entanglement
This repository contains the code required to reproduce the plots in "[Improved analytical bounds on delivery times of long-distance entanglement](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.105.012608)" ([arXiv](https://arxiv.org/abs/2103.11454)).

## Files overview
The files are structured as follows:

* `bounds_plots.py` contains the main code to actually generate the plots, as well as the implementation of the formulas in the paper.

The following files

* `repeater_chain_analyzer.py` 
* `plotting_tools.py`
* `probability_tools.py`
* `repeater_chain_analyzer.py`
* `werner_tools.py`

contain the implementation of the numerical numerical calculations from [previous work](https://ieeexplore.ieee.org/abstract/document/8972391) ([arXiv](https://arxiv.org/abs/1912.07688)), and are originally from [this repository](https://github.com/sebastiaanbrand/waiting-time-quantum-repeater-chains).

## Dependencies
The code is written for `Python 3` and requires the following packages:
```
numpy
matplotlib
scipy
```

## Reproducing the plots
The plots contain comparisons between numerical and analytical results. When generating the plots, the numerical results are first calculated and written to csv files in `results/`. This can take some time.

### Figure 4
The plot in Fig. 4 can be generated with
```
$ python bounds_plots.py mean
```
Running the numerical calculation can take up to ~40 minutes (depending on hardare). The plot is written to `results/mean_bounds_ratios_plot.pdf`.

### Figure 5
The plots in Fig. 5 can be generated with
```
$ python bouds_plots.py tail
```
Running the numerical calculation can take up to ~15 minutes. The plots are written to  `results/tail_bounds_plot_pswap0.2.pdf` and `results/tail_bounds_plot_pswap0.5.pdf`.
