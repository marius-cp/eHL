Code for technical report ‘A safe Hosmer-Lemeshow test’
================

### Citation

T. Dimitriadis, A. Henzi, M. Puke, and J. Ziegel. A safe Hosmer-Lemeshow
test. *Preprint*,
[*arXiv:2203.00426*](https://doi.org/10.48550/arXiv.2203.00426), 2022.

### eHL Simulation

-   The code for the simulation is in `sim_quad_misspec.R`. The
    simulation script uses the functions in `nle.R`, `eHL.R`, and
    `oracle_eHL.R`.

-   `nle.R`: Data generating processes for simulation studies of
    Hosmer-Lemeshow type goodness-of-fit test follow a particular
    (traditional) pattern. This function allows to calculate the
    parameters *β* described in section 3.1.

-   `eHL.R`: Calculates the e-value based on the sample split and the
    isotonic/monotonic recalibrated probabilities.

-   `oracle_eHL.R`: Calculates the optimal (hypothetical) e-value

-   Fig 1: created by the `Fig_1-2.R`

-   Fig 2: created by the `Fig_1-2.R`

-   Fig 3: created by the `Fig_3.R`
