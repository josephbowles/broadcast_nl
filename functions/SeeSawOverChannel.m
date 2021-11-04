function [newChoiMap,finalObj,problemStatus] = SeeSawOverChannel(state, bellcoeffs, povms)
    dimA = 2;
    dimB = 2;
    dimB1 = 2;
    dimB2 = 2;
    idB = eye(dimB);
    idB1 = eye(dimB1);
    idB2 = eye(dimB2);
    idB1B2 = Tensor(idB1,idB2);
    
    assert(mod(length(size(bellcoeffs)),2)==0,"There should be as many inputs as outputs.");
    dims = size(bellcoeffs);
    nrparties = length(dims)/2;
    ins = dims(1:nrparties);
    outs = dims(nrparties+1:end);

    inputdimspace = 2;
    outputdimspace = 4;
    choidim = inputdimspace*outputdimspace;
    
    choi = sdpvar(choidim,choidim,'hermitian','complex');

    summ = 0;
    for x=1:ins(1)
        for y=1:ins(2)
            for z=1:ins(3)
                for a=1:outs(1)
                    for b=1:outs(2)
                        for c=1:outs(3)
                            term = Tensor(povms{1}{x}{a}, ...
                                            povms{2}{y}{b},...
                                            povms{3}{z}{c});
                            summ = summ + bellcoeffs(x,y,z,a,b,c)*term;
                        end
                    end
                end
            end
        end        
    end
    
%     state_small = ini_state(alpha);
%     state_small_reshaped = reshape(state_small, 2,2,2,2);
%     
%     % TODO change this to a cleaner form, check that it gives the same
%     state = 0;
%     for i=1:dimA
%         for j=1:dimB
%             for k=1:dimA
%                 for l=1:dimB
%                     %scnd = PartialTrace( choi * Tensor( ketbra(j,l,dimB).', idB1B2),...
%                     %                    1, [dimB,dimB1*dimB2] );
%                     state = state + state_small_reshaped(i,j,k,l) * ...
%                                     kron(ketbra(i,k,dimA), ApplyMap(ketbra(j,l,dimB), choi));
%                 end
%             end
%         end
%     end
    
    biggerstate = kron( state.', eye(dimA*dimB1*dimB2) );
    Phi = auxPHI(dimA);
    biggerchannel = kron(Phi*Phi',choi);
    % after previous line tensor spaces are A x A x B x B1 x B2, we need
    % to swap the 2nd and 3rd
    swapop=Tensor(eye(dimA),SwapOperator(2),eye(dimB1*dimB2));
    biggerchannel = swapop * biggerchannel * swapop';
    output_state = PartialTrace(biggerchannel*biggerstate, [1,2], [dimA,dimB,dimA,dimB1,dimB2]);
    %state = AppleMap(inistate, biggerchannel);  % either this or the
                                                % previous line work well
    
    constraints = [choi >= 0];
    constraints = [constraints, PartialTrace(choi, 2, [dimB, dimB1*dimB2]) == idB];
    
    objective = real(trace(output_state*summ));
    
    sol = optimize(constraints, -objective, ...
                    sdpsettings('solver','mosek','verbose', 0, ...
                                'dualize', 0, 'showprogress', 0,...
                                'savesolverinput', 0,'savesolveroutput', 0, ...
                                'debug', 0, 'warning',0));
    if sol.problem ~= 0
       disp(sol);
       error('Check what problem there is.'); 
    end
    newChoiMap = value(choi);
    finalObj = value(objective);
    problemStatus = sol.problem;
    
    %list = whos;for i = 1:length(list);if strcmp(list(i).class,'sdpvar');clear(list(i).name);end;end
    %yalmip("clear");
end
