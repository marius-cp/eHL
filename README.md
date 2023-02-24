Replication material ‘A safe Hosmer-Lemeshow test’
================

## Citation

A. Henzi, M. Puke, T. Dimitriadis, and J. Ziegel (2023). A safe
Hosmer-Lemeshow test. *Preprint*,
[*arXiv:2203.00426*](https://doi.org/10.48550/arXiv.2203.00426).

## Simulation readme

- The simulation code can be found in the simulation folder. To execute
  the simulation, please run the `sim_HLstyle_oos.R` script. The results
  from each replication will be saved in an rds file. In order to carry
  out the simulation, the script relies on the functions contained in
  the `nle.R`, `eHL.R`, and `oracle_eHL.R` files located in the main
  folder. Those files do the following:

  - `nle.R`: The `nle.R` function facilitates the generation of data for
    the simulation study of the Hosmer-Lemeshow type calibration test.
    Using this function, one can calculate the parameters $\beta$, which
    are described in section 3 of the document.

  - `eHL.R`: The `eHL.R` function is responsible for calculating the
    e-value by utilizing the sample split and the isotonic/monotonic
    recalibrated probabilities. It is based on Algorithm 2, which is
    described in Section 2.5 of the document.

  - `oracle_eHL.R`: Calculates the optimal (hypothetical) e-value.

- To obtain the results presented in Table 1 of the document, please
  refer to the `table.R` script. Similarly, to generate Figure 1, please
  run the `plots_simulation_main.R` script.

- The `sim_instability_appendix.R` script generates the results
  presented in Figure 4 of the document. To access the results from each
  run, please refer to the corresponding `rds` file.

## Application readme

- In application, we analyze (re-)calibration of probability predictions
  for the binary event of credit card defaults in Taiwan in the
  year 2005. The results can be replicated using the file `credit.R`.

- The data analyzed here was obtained from the [UCI Machine Learning
  Repository](https://archive.ics.uci.edu/ml/datasets/default+of+credit+card+clients).
  To reproduce the study’s findings, it is necessary to access and
  download the dataset from the aforementioned webpage.

## Computational requirements

The software versions that were used to run these analyses are

- R 4.2.1
  - `readxl` (1.4.1)
  - `dplyr` (1.0.10)
  - `tidyr` (1.2.1)
  - `ggplot2` (3.4.0)
  - `ResourceSelection` (0.3-5)
  - `kableExtra` (1.3.4)
  - `foreach`(1.5.2)
  - `tidyverse` (1.3.2)
  - `doParallel` (1.0.17)

The time-intensive simulations were performed on 100 cores of the
cluster of the [data
laboratory](https://cfh.uni-hohenheim.de/en/compute-server) of the
University of Hohenheim. See chapter `Simulations` for details. We are
grateful for the provided computing power.

## References

Dua, D. and Graff, C. (2019). UCI Machine Learning Repository
\[<http://archive.ics.uci.edu/ml>\]. Irvine, CA: University of
California, School of Information and Computer Science.

R Core Team (2023). R: A language and environment for statistical
computing. R Foundation for Statistical Computing, Vienna, Austria. URL
<https://www.R-project.org/>.

Yeh, I. C., & Lien, C. H. (2009). The comparisons of data mining
techniques for the predictive accuracy of probability of default of
credit card clients. Expert Systems with Applications, 36(2), 2473-2480.
