function probbellexpression = ToProjectorNotation(corrbellexpression, ins, outs, FLAG_Use01obsInsteadOfCorrelator)
%ToProbabilityNotationTerm Summary of this function goes here
%   Detailed explanation goes here

expression = expand(corrbellexpression);

[C,T] = coeffs(expression);

nrterms = length(T);

summ = 0;
for termidx=1:nrterms
    summ = summ + C(termidx)*ToProjectorNotationTerm(T(termidx),ins, outs, FLAG_Use01obsInsteadOfCorrelator);
end

probbellexpression = expand(summ);

end

