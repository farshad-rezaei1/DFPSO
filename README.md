# DFPSO
In any meta-heuristic algorithm, each search agent must move to the high-fitness areas in the search space while preserving its diversity. At first glance, there is no relationship between fitness and diversity, as two key factors to be considered in selecting a guide for the solutions. In other words, each of these factors must be evaluated in its specific and independent way. Since the independent ways to evaluate the fitness and diversity usually make any meta-heuristic consider these factors disproportionately to choose the guides, the solutions’ movements may be unbalanced. In this project, a novel version of the Particle Swarm Optimization (PSO) algorithm, named Dual Fitness PSO (DFPSO) is proposed. In this algorithm, not only fitness and diversity of the particles are properly evaluated, but also the abilities to evaluate these features are integrated to avoid the abovementioned problem in determining the global guide particles. 

For further information about this algorithm, please refer to the following reference:

Rezaei, F., Safavi, H.R. Sustainable Conjunctive Water Use Modeling Using Dual Fitness Particle Swarm Optimization Algorithm. Water Resour Manage 36, 989–1006 (2022). https://doi.org/10.1007/s11269-022-03064-w. 

This article introduces the DFPSO and its unique characteristics discriminating it from the other meta-heuristics, in detail. Then the performance of the algorithm is tested on a variaty of the standard benchmark functions and compared with a set of other well-known PSO variants as well as a set of popular independent meta-heuristic optimization algorithms. Finally, the DFPSO is successfully applied to solve a real-world high-constrained engineering optimization problem. The results revealed high competence of the proposal to handle a vast range of the difficulties optimization problems may be engaged with, and efficiently solve such problems. 

Please cite this article upon using this source code.
