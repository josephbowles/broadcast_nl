## Code to accompany: [Device-independent and semi-device-independent entanglement certification in broadcast Bell scenarios](https://arxiv.org/abs/2111.xxxx)
#### Emanuel-Cristian Boghiu, Flavien Hirsch, Pei-ShengLin, Marco Túlio Quintino, and Joseph Bowles

This is a repository for the code used to calculate the numerical results presented in the article [Device-independent and semi-device-independent entanglement certification in broadcast Bell scenarios](https://arxiv.org/abs/2111.xxxx).

 MATLAB code requires:
- [cvx](http://cvxr.com/) - a free MATLAB toolbox for rapid prototyping of optimization problems.
- [YALMIP](https://github.com/yalmip) - a free MATLAB toolbox for optimization modeling.
- [QETLAB](http://www.qetlab.com/) - a free MATLAB toolbox for quantum entanglement theory.

The MATLAB code of this repository contains:

- [PrepareAndRun_HeuristicSearchBroadcast.m](https://github.com/josephbowles/broadcast_nl/blob/main/PrepareAndRun_HeuristicSearchBroadcast.m):
Script which is read to perform a Heuristic Search method to find the optimal visibility of a given quantum state in the broadcast nonlocal scenario. This script has several adjustable parameters such as, setting a target state and the relative noise, number of parties, number of input per parties, number of outputs.

- [PrepareAndRun_SteeringHeuristicSearchBroadcast2Bobs.m](https://github.com/josephbowles/broadcast_nl/blob/main/PrepareAndRun_SteeringHeuristicSearchBroadcast2Bobs.m):
Script which is read to perform a Heuristic Search method to find the optimal visibility of a given quantum state in the broadcast steering scenario with 2 untrusted Bobs. This script has various adjustable parameters such as, setting a target state and the relative noise, number of input per parties, number of outputs.

- [PrepareAndRun_SteeringHeuristicSearchBroadcast3Bobs.m](https://github.com/josephbowles/broadcast_nl/blob/main/PrepareAndRun_SteeringHeuristicSearchBroadcast3Bobs.m):
Script which is read to perform a Heuristic Search method to find the optimal visibility of a given quantum state in the broadcast steering scenario with 3 untrusted Bobs. This script allows the user to set a target state and a relative noise state.

- Various subroutines and useful functions which are commented in a way to allow the users to adapt it for similar problems.

- The folder scripts_Sec3A3CSec4 which contains the code and required variables to verify the results presented in Section III


