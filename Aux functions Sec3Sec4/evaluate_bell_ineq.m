function out = evaluate_bell_ineq(bellcoeffs, belloffset, finalstate, povms, ins, outs)
    %assert(mod(length(size(bellcoeffs)),2)==0,"There should be as many inputs as outputs.");
    probarray = ProbMultidimArray(finalstate, povms, ins, outs);
	aux = probarray.*bellcoeffs;
    out = sum(aux(:)) + belloffset;
end

