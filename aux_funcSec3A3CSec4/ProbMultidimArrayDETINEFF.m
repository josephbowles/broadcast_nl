function probability_ndarray = ProbMultidimArrayDETINEFF(state,povms, etas, ins, outs)
    nrparties = length(ins);
    
    
    noisy_povms = cell(nrparties, max(ins), max(outs));%povms;
    for p=1:nrparties
        for x=1:ins(p)
           for a=1:(outs(p)+1)
              if a ~= outs(p)+1
                  noisy_povms{p,x,a} = etas(p)*povms{p}{x}{a};
              else
                  noisy_povms{p,x,outs(p)+1} = (1-etas(p))*eye(outs(p));
              end
           end
        end
    end
    
    outs_plus_one = outs + 1;
    
    % Define a zero multidim where to place the probability distribution
    dims = num2cell([ins,outs_plus_one]);
    probability_ndarray = zeros(dims{:});
        
    aux = [ins, outs_plus_one];
    allinputoutputcombinations = ind2subv(aux, 1:prod(aux(:)));
    for slice=1:size(allinputoutputcombinations,1)
        ins_slice = num2cell(allinputoutputcombinations(slice,1:nrparties));
        outs_slice = num2cell(allinputoutputcombinations(slice,nrparties+1:end));

        tensor = 1;
        for p=1:nrparties
             tensor = kron(tensor,noisy_povms{p,ins_slice{p},outs_slice{p}});
        end

        probability_ndarray(ins_slice{:}, outs_slice{:}) = real(trace(tensor*state));
    end
end