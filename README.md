# Decomposition of Differentials in Health Expectancies from Multistate Life Tables
***Tianyu Shen, Tim Riffe, Collin F. Payne and Vladimir Canudas-Romo***

This is a repository for our paper in [*Demography*](https://doi.org/10.1215/00703370-11058373).

Details of the mathematical derivations and the demographic interpretation of the results are included in our paper.

The code is available in "**code.R**".

Our example data are authors' calculated from the HRS. They are stored in folder "**RMLE**". Inside this folder, there are two files:

- "**BASELINE.csv**" includes the initial health structure of the population. The columns are explained below

**state**: health state (with 1 and 2, healthy and unhealthy respectively)

**ragender**: gender

**pro**: proportion of the population in that health state (they are rescaled to 1 by sex and iteration to calculated the HLE for male and female separately)

**iter**: bootstrap iteration number

- "**PROB.csv**" includes the transition probabilities by age. The columns are explained below

**pre_state**: the initial state

**ragender**: gender

**age**: age

**iter**: bootstrap iteration number

**A**: probability to "Healthy" given the initial state

**L**: probability to "Unhealthy" given the initial state

**H**: probability to "Death" given the initial state

In the R output: Column "A" represent the healthy life expectancy, Column "L" the unhealthy life expectancy and Column "H" the death expectancy (time in death). We didn’t report this value because it is only part of the calculation and not meaningful. We are computing the temporary life expectancy, so in the example, the maximum life expectancy is 51 (i.e., 55-105). Within the 51 years, female would be dead for 24.33 years. The total life expectancy is 51-24.33=26.68=22.47+4.21. The only possible way to interpret this value is the “Year of Life Lost” (YYL), but it is not very useful for anyone to know how long one spends dead. Unless you have a state space without absorbing state, in which all the column has its meaning.

For other enquiries related to the paper please email tianyu.shen@anu.edu.au
