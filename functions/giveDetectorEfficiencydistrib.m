function probability_ndarray = giveDetectorEfficiencydistrib(input_prob, efficiencies, ins, outs)
% This function, given an input array indexed as p(x,y,z,a,b,c) where
% p==input_prob, gives the distribution P(x,y,z,a',b',c') where
% a'\in{{a},FAIL}, b' \in {{b}, FAIL}, c' \in {{c},FAIL}
% i.e., we add an extra output which correspodns to detector failure.
% This is calculated as follows for a 2 party distribution for simplicity:
% P(a,b|x,y)=eta^2 p(a,b|x,y)
% P(FAIL,b|x,y) = P(FAIL|x,y)*P(b NOT FAIL|x,y)*p(b|y) = eta*(1-eta)*p(b|y)
% P(a, FAIL|x,y) = eta*(1-eta)*p(a|x)
% P(FAIL, FAIL|x,y)=(1-eta)^2
% For other projects that use this definition see arXiv:1812.10017v3
% The generalisation to 3 parties is straightforwards. This is what this
% code does.


outs_plus_one = outs + 1; % Add the failure mode. a in {0,1,2,...,d} and {d+1} will be detector failure

nr_parties = length(outs);

assert(length(efficiencies)==nr_parties,"Incorrect input");

dims = [ins, outs_plus_one];
dims_cell = num2cell(dims);
probability_ndarray = zeros(dims_cell{:});

prob_fail = 1 - efficiencies;
unos = 1:nr_parties;

aux = [ins, outs_plus_one];
allinputoutputcombinations = ind2subv(aux, 1:prod(aux(:)));
for slice=1:size(allinputoutputcombinations,1)
    ins_slice = num2cell(allinputoutputcombinations(slice,1:nr_parties));
    outs_slice = num2cell(allinputoutputcombinations(slice,nr_parties+1:end));
    num_outs_slice = [outs_slice{:}];
    
    % given this slice, get the index of the parties whose detectors fail
    parties_whose_detectors_fail = num_outs_slice==outs_plus_one;
    parties_whose_detectors_work = ~parties_whose_detectors_fail;
                

    
    if any(parties_whose_detectors_fail)
        % in this case there is at least one detector failur
        % i.e., in P(a'b'c'|xyz) at least one of the outputs corresponds to
        % FAIL. now we need to identify which parties fail and write the
        % according expression
        indices_where_it_fails = unos(parties_whose_detectors_fail);
        indices_where_it_detects = unos(parties_whose_detectors_work);
    
        product_of_probs_of_fail = prob_fail(indices_where_it_fails);
        product_of_probs_of_fail = prod(product_of_probs_of_fail(:));
        
        product_of_probs_of_detection = efficiencies(indices_where_it_detects);
        product_of_probs_of_detection = prod(product_of_probs_of_detection(:));
        
        %disp([efficiencies, efficiencies(indices_where_it_detects)])
        %disp([prob_fail(indices_where_it_fails)])
        %disp([ins_slice{:}, outs_slice{:}, product_of_probs_of_detection, product_of_probs_of_fail, prob_fail(indices_where_it_fails), ])
        
        % assume shape(prob) = [ 3 3 3 2 2 2] and that we have failure for b.  
        % then P(a FAIL c|xyz) = etaA (1-eta_B) etaC p_{2party} (a c |x z)
        % when we do sum(prob, [2]) the final shape of prob goes from 
        % [3 3 3 2 2 2 ] to [3 1 3 2 2 2]. The point is that in our
        % selection of outputs we need to set output=1 for all outputs that
        % fail to detect
        % i.e., p_{2party} (a c |x z) is an array that gets indexed as 
        % p2party(x,y,za,1,c) and I need to add the '1' for the failed
        aux_outs_slice = [outs_slice{:}];
        aux_outs_slice(indices_where_it_fails) = 1; 
        
        aux_outs_slice = num2cell(aux_outs_slice);
        marginal_over_failed = sum(input_prob, nr_parties + indices_where_it_fails);

        probability_ndarray(ins_slice{:}, outs_slice{:}) = ... 
                product_of_probs_of_fail * product_of_probs_of_detection * ...
                marginal_over_failed(ins_slice{:}, aux_outs_slice{:});
    else
        % in this case no failure, P(a'b'c'|xyz)=eta_A*eta_B*eta_C *
        % p(a'b'c'|xyz)
        probability_ndarray(ins_slice{:},outs_slice{:}) = prod(efficiencies(:))*input_prob(ins_slice{:},outs_slice{:});
    end
        %if all(parties_whose_detectors_work==0)
       %error("ASD"); 
        %end
end

end
