function res = eval_eta(bellcoeffs, ini_povms, channel, eta, lamA, lamB, lamC, ins, outs)
% auxiliary function to find the value of the score of the bell inequality
% when some of the detectors fail with probability 1-eta
oldPOVMs = ini_povms;
newPOVMs = eta_POVMs(oldPOVMs,1,lamA,eta,ins,outs);
newPOVMs = eta_POVMs(newPOVMs,2,lamB,eta,ins,outs);
newPOVMs = eta_POVMs(newPOVMs,3,lamC,eta,ins,outs);

p_entangled = ProbMultidimArray(final_state(NoisyWernerState(0), channel), newPOVMs, ins, outs);
res = sum(bellcoeffs.*p_entangled, 'all');
end
