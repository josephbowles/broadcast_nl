function opt = optimizer_povm_party(partyidx, ins, outs)
    fprintf("Calculating YALMIP's 'optimizer' for the POVM of party %d. It can take a while.\n", partyidx);
    nrparties = length(ins);

    dimA = outs(1);
    dimB1 = outs(2);
    dimB2 = outs(3);
    
    final_state_dim = dimA*dimB1*dimB2;
    final_state = sdpvar(final_state_dim,final_state_dim,'hermitian','complex');
    ia_state = give_ia_state(final_state);
    
    allbutone = logical(ones(nrparties,1));
    allbutone(partyidx) = false;
    prod_outs_without_partyidx = prod(outs(allbutone));
    
    povms_party = cell(ins(partyidx), outs(partyidx));
    partial_products = cell(ins(partyidx), outs(partyidx));
    %ia_partial_products = cell(ins(partyidx), outs(partyidx));
    for x=1:ins(partyidx)
        for a=1:outs(partyidx)
            povms_party{x,a} = sdpvar(outs(partyidx),outs(partyidx),'hermitian','complex');
            partial_products{x,a} = sdpvar(prod_outs_without_partyidx,prod_outs_without_partyidx,'full','complex');  % the combination with the bell coefficients can be arbitrary
            %ia_partial_products{x,a} = partial_products{x,a}(find(triu(ones(prod_outs_without_partyidx)))); 
        end
    end

    constraints_povm = [];
    for x=1:ins(partyidx)
        summ = 0;
        for a=1:outs(partyidx)
            summ = summ + povms_party{x,a};
            constraints_povm = [constraints_povm,  povms_party{x,a} >= 0];
        end
        constraints_povm = [constraints_povm, summ == eye(outs(partyidx))];
    end

    if partyidx == 2
        swapop = Tensor(SwapOperator([outs(2) outs(1)]), eye(outs(3)));
    end
    
    BellOperator = 0;
    for x=1:ins(partyidx)
       for a=1:outs(partyidx)
           if partyidx == 1  % Alice
               BellOperator = BellOperator + kron(povms_party{x,a}, partial_products{x,a});
           elseif partyidx == 2  % Bob
               % The idea is that we input 'partial_products' which
               % will be products of Alice times Charlie but we want to
               % tensor this with Bob. I do Bob x Alice x Charlie and then
               % use the 'swapop' to permute to Alice x Bob x Charlie
               BellOperator = BellOperator + swapop * kron(povms_party{x,a}, partial_products{x,a}) * swapop';
           elseif partyidx == 3  % Charlie
               BellOperator = BellOperator + kron(partial_products{x,a}, povms_party{x,a});
           else
              error('not supported for more than 3 parties'); 
           end
       end
    end
    
    objective = real(trace(final_state*BellOperator));
    
    opt_povm_party = optimizer(constraints_povm, -objective, ...
                         sdpsettings('solver','mosek'), ...
                        [{ia_state}, partial_products(:)'], ...
                        [{objective}, povms_party(:)']);

    opt = opt_povm_party;
end