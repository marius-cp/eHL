# eHL
Code for technical report "A safe Hosmer-Lemeshow test". 


### Simulation

* The code for the simulation is in `sim_quad_misspec.R`. The simulation script uses  the functions in `nle.R`, `eHL.R`, and `oracle_eHL.R`. 

* `nle.R`: data generating processes for simulation studies of Hosmer-Lemeshow type goodness-of-fit test follow a particular (traditional) pattern. This function allows to calculate the parameters $\beta$ described in section 3.1.

* `eHL.R`: calculates the e-value based on the sample split and the isotonic/monotonic recalibrated probabilities. 

* `oracle_eHL.R`: calculates the optimal (hypothetical) e-value 


* Fig 1 created by the   `Fig_1-2.R`

<img src="plots/Fig_1.pdf" />


* Fig 2: created by the   `Fig_1-2.R`


<img src="plots/Fig_2.pdf" />


*  Fig 3 created by the   `Fig_3.R`

<img src="plots/Fig_3.pdf" />
