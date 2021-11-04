function [summ, scaling] = dispBellCoeffsCorrelators(bellcoeffs,ins,outs) 
    cartproductOUT = ind2subv(outs, 1:prod(outs(:)));
    cartproductIN  = ind2subv(ins,  1:prod(ins(:)));
    
    summ = 0;
    for i1 = 1:size(cartproductOUT,1)
        outputs = cartproductOUT(i1,:);
        for i2 = 1:size(cartproductIN,1)
            settings = cartproductIN(i2,:);
            indexes = num2cell([settings, outputs]);
            summ = summ + bellcoeffs(indexes{:}) * ToCorrelatorNotationWithInOut(settings,outputs);
        end
    end
    summ = expand(summ);
    
    
    % clean the terms
    %tolerance = 1e-7;
    %[C,T] = coeffs(summ);
    %C(abs(C)<tolerance)=0;
    %summ = dot(C,T);
    
    %scaling = 1;
    % normalize by smallest one
    %[C,T] = coeffs(summ);  % redo the split so there are no 0 coeffs to divide by
    %scaling = 1.0/min(abs(C(:)));
    %C = scaling*C;
    %summ = dot(C,T);
    
    %summ = vpa(summ,3); % vpa to simplify integer fractions
end

