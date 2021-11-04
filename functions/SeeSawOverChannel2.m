function [newChoiMap,finalObj,problemStatus] = SeeSawOverChannel2(state, bellcoeffs, povms, channels, position, ins, outs)
    % This is a modification of the unnumbered function for the 4 party scenario
    % with 2 alices and 2 bobs.

    dimA = 2;
    dimB = 2;
    dimB1 = 2;
    dimB2 = 2;
    idB = eye(dimB);
    idB1 = eye(dimB1);
    idB2 = eye(dimB2);
    idB1B2 = Tensor(idB1,idB2);
    

    inputdimspace = 2;
    outputdimspace = 4;
    choidim = inputdimspace*outputdimspace;
    
    choi = sdpvar(choidim,choidim,'hermitian','complex');

    summ = 0;
    for x=1:ins(1)
        for y=1:ins(2)
            for z=1:ins(3)
                for w=1:ins(4)
                    for a=1:outs(1)
                        for b=1:outs(2)
                            for c=1:outs(3)
                                for d=1:outs(4)
                            term =   Tensor(povms{1}{x}{a}, ...
                                            povms{2}{y}{b}, ...
                                            povms{3}{z}{c}, ...
                                            povms{4}{w}{d});
                            summ = summ + bellcoeffs(x,y,z,w,a,b,c,d)*term;
                                end
                            end
                        end
                    end
                end
            end
        end        
    end
    
    if position=='A'
        output_state = final_state2(state, choi, channels{2});
    elseif position=='B'
        output_state = final_state2(state, channels{1}, choi);
    else
        error('bad position argument');
    end
    
    constraints = [choi >= 0];
    constraints = [constraints, PartialTrace(choi, 2, [dimB, dimB1*dimB2]) == idB];
    
    objective = real(trace(output_state*summ));
    
    sol = optimize(constraints, -objective, ...
                    sdpsettings('solver','mosek','verbose', 0));
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
