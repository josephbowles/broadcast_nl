## Code to accompany: [Device-independent and semi-device-independent entanglement certification in broadcast Bell scenarios](https://arxiv.org/abs/2111.06358)
#### Flavien Hirsch, Emanuel-Cristian Boghiu, Pei-Sheng Lin, Marco Túlio Quintino, and Joseph Bowles

This is a repository for the codes used to calculate the numerical results presented in the article [Device-independent and semi-device-independent entanglement certification in broadcast Bell scenarios](https://arxiv.org/abs/2111.06358).

 MATLAB code requires:
- [CVX](http://cvxr.com/) - a free MATLAB toolbox for rapid prototyping of optimization problems.
- [YALMIP](https://github.com/yalmip) - a free MATLAB toolbox for optimization modeling.
- [QETLAB](http://www.qetlab.com/) - a free MATLAB toolbox for quantum entanglement theory.
- [MOSEK](https://www.mosek.com/) - a solver for linear and semidefinite programs with free academic licenses; the code can be perfectly adapted to use free and open sources alternatives such as [CVXOPT](https://cvxopt.org/).

The MATLAB code of this repository contains:

- [PrepareAndRun_HeuristicSearchBroadcast.m](https://github.com/josephbowles/broadcast_nl/blob/main/PrepareAndRun_HeuristicSearchBroadcast.m):
Script which is read to perform a Heuristic Search method to find the optimal visibility of a given quantum state in the broadcast nonlocality scenario. This script has several adjustable parameters such as, setting a target state and the relative noise, number of parties, number of input per parties, number of outputs.

- [PrepareAndRun_SteeringHeuristicSearchBroadcast2Bobs.m](https://github.com/josephbowles/broadcast_nl/blob/main/PrepareAndRun_SteeringHeuristicSearchBroadcast2Bobs.m):
Script which is read to perform a Heuristic Search method to find the optimal visibility of a given quantum state in the broadcast steering scenario with 2 untrusted Bobs. This script has various adjustable parameters such as, setting a target state and the relative noise, number of input per parties, number of outputs.

- [PrepareAndRun_SteeringHeuristicSearchBroadcast3Bobs.m](https://github.com/josephbowles/broadcast_nl/blob/main/PrepareAndRun_SteeringHeuristicSearchBroadcast3Bobs.m):
Script which is read to perform a Heuristic Search method to find the optimal visibility of a given quantum state in the broadcast steering scenario with 3 untrusted Bobs. This script allows the user to set a target state and a relative noise state.

- Various subroutines and useful functions which are commented in a way to allow the users to adapt it for similar problems.

- [Scripts Sec3Sec4](https://github.com/ecboghiu/broadcast_nl/tree/main/Scripts%20Sec3Sec4): 
This is a folder with various scripts reproducing results from Sections 3 and 4. [promote_chained_1A2B.m](https://github.com/ecboghiu/broadcast_nl/blob/main/Scripts%20Sec3Sec4/promote_chained_1A2B.m) and [promote_elegant_1A2B.m](https://github.com/ecboghiu/broadcast_nl/blob/main/Scripts%20Sec3Sec4/promote_elegant_1A2B.m), using the Symbolic Math Toolbox™, calculate the coefficients of a Bell inequality promoted to the broadcast scenario from Fig. 1a) through the construction given in the paper. [PromotedElegant_visibility.m](https://github.com/ecboghiu/broadcast_nl/blob/main/Scripts%20Sec3Sec4/PromotedElegant_visibility.m) and [PromotedChained_visibility.m](https://github.com/ecboghiu/broadcast_nl/blob/main/Scripts%20Sec3Sec4/PromotedChained_visibility.m) optimize said inequalities over measurements and channel to find how robust they are to white noise when the noiseless quantum state is a maximally entangled two-qubit state. [PromotedIneqs_visibility_4party.m](https://github.com/ecboghiu/broadcast_nl/blob/main/Scripts%20Sec3Sec4/PromotedIneqs_visibility_4party.m) finds the robustness to white noise of the inequalities promoted to the 4-partite broadcast scenario from Fig 1b). [detectionefficiency_07355.m](https://github.com/ecboghiu/broadcast_nl/blob/main/Scripts%20Sec3Sec4/detectionefficiency_07355.m) verifies the detection efficiency threshold of 0.7355. [NS2activation_POVMstate_ineq16.m](https://github.com/ecboghiu/broadcast_nl/blob/main/Scripts%20Sec3Sec4/NS2activation_POVMstate_ineq16.m) verifies the broadcast activation of NS genuine network nonlocality of a bipartite state with a local model for all general measurements.

- Scripts to verify the numerical results for broadcast steering and device-independent entanglement certification, found in the directories BroadcastSteering and entanglement_detection. For the entanglement detection results, the code RebuildPPTResult.m verifies the numerical claim in the paper. 

