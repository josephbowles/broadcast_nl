function flag=checkThatProbSumsToOne(prob,ins,outs)
flag = true;
tol  = 1e-7;

cartproductOUT = ind2subv(outs, 1:prod(outs(:)));
cartproductIN  = ind2subv(ins,  1:prod(ins(:)));
for i2 = 1:length(cartproductIN)
    inputs = cartproductIN(i2,:);
    summ = 0;
    for i1 = 1:length(cartproductOUT)
        outputs = cartproductOUT(i1,:);
        indexes = num2cell([inputs,outputs]);
        if prob(indexes{:})<-tol
            warning("Probability not positive. prob(%d %d %d %d %d %d)=%g", indexes{:}, prob(indexes{:}));
        end
        summ = summ + prob(indexes{:});
    end
    if abs(summ-1)>tol
        warning("Something went wrong. sum_outputs(prob) = %g, (x,y,z)=",summ);
        disp(inputs);
        flag = false;
        return
    end
end        
                        
end

