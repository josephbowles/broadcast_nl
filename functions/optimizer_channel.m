function out_optimizer = optimizer_Channel(ins, outs, inputdimspace)
    fprintf("Calculating YALMIP's 'optimizer' for optimization over channels. It can take a while depending on the size of the problem.\n");
    
    nrparties = length(ins);

    dimA = outs(1);
    dimB = inputdimspace;
    dimB1 = outs(2);
    dimB2 = outs(3);
    idB = eye(dimB);
    idB1 = eye(dimB1);
    idB2 = eye(dimB2);
    idB1B2 = Tensor(idB1,idB2);

    outputdimspace = dimB1*dimB2;
    choidim = inputdimspace*outputdimspace;
    
    state_dim = dimA*dimB;
    state = sdpvar(state_dim,state_dim,'hermitian','complex');
    ia_state = give_ia_state(state);
    
    choi = sdpvar(choidim,choidim,'hermitian','complex');
    ia_choi = give_ia_state(choi);
 
    state_small = state; %ini_state(alpha);
    state_small_reshaped = reshape(state_small, dimA,dimB,dimA,dimB);

    biggerstate = Tensor(state.', eye(dimB1*dimB2) );
    biggerchannel = kron(eye(dimA), choi);
    output_state = PartialTrace(biggerchannel*biggerstate, 2, [dimA,dimB,dimB1,dimB2]);

    constraints_choi = cell(2);
    constraints_choi{1} = choi >= 0;
    constraints_choi{2} = PartialTrace(choi, 2, [dimB, dimB1*dimB2]) == idB;

    BellOperator = sdpvar(dimA*dimB1*dimB2,dimA*dimB1*dimB2,'full','complex');
    
    objective = real(trace(BellOperator*output_state));
   
    out_optimizer = optimizer([constraints_choi{:}], -objective, ...
                         sdpsettings('solver','mosek'), ...
                        [{ia_state}, {BellOperator}], ...
                        {choi, objective});
                    
    % WARNING: do NOT use 'yalmip('clear')' as that kills the variables to which
    % 'optimizer' makes reference
end