function [newPovms,finalAlpha,problemStatus] = SeeSawOverASingleParty2(partyidx, state, bellcoeffs, povms, ins, outs)
    % This is a modification of the unnumbered function for the 4 party scenario
    % with 2 alices and 2 bobs.

    nrparties = length(ins);

    %% Define the SDP variables.
    chosenPartyMeasurements = {{}};
    positivityconstraints  = [];
    dimofProjector = outs(partyidx);
    for x=1:ins(partyidx)
        for a=1:outs(partyidx)
            chosenPartyMeasurements{x}{a} = sdpvar(dimofProjector,dimofProjector,'hermitian','complex');
            positivityconstraints = [positivityconstraints, ...
                chosenPartyMeasurements{x}{a} >= 0];
        end
    end
    
    %% POVM constraints
    povmconstraints = [];
    for x=1:ins(partyidx)
        summ = 0;
        for a=1:outs(partyidx)
            summ = summ + chosenPartyMeasurements{x}{a};
        end
        povmconstraints = [povmconstraints, summ == eye(dimofProjector)];
    end

%     %% Define the objective function.
%     coords_structure = [ins, outs];
%     coords = ind2subv(coords_structure, 1:prod(coords_structure(:)));
%     summ = 0;
%     for row=1:size(coords,1)
%        settings = num2cell(coords(row,1:nrparties));
%        outputs = num2cell(coords(row,nrparties+1:end));
%        
%        term = 1;
%        for party=1:nrparties
%             if party == partyidx
%                 term = kron(term, povms{party}{settings{party}}{outputs{party}});
%             else
%                 term = kron(term, chosenPartyMeasurements{settings{party}}{outputs{party}});
%             end
%        end
%        summ = summ + bellcoeffs(settings{:},outputs{:}) * term;
%     end
summ = 0;
    for x=1:ins(1)
        for y=1:ins(2)
            for z=1:ins(3)
                for w = 1:ins(4)
                    for a=1:outs(1)
                        for b=1:outs(2)
                            for c=1:outs(3)
                                for d=1:outs(4)
                            if partyidx == 1
                               term = Tensor(chosenPartyMeasurements{x}{a}, ...
                                            povms{2}{y}{b},...
                                            povms{3}{z}{c},...
                                            povms{4}{w}{d});
                            elseif partyidx == 2
                               term = Tensor(povms{1}{x}{a}, ...
                                            chosenPartyMeasurements{y}{b},...
                                            povms{3}{z}{c},...
                                            povms{4}{w}{d});
                            elseif partyidx == 3
                                term = Tensor(povms{1}{x}{a}, ...
                                            povms{2}{y}{b},...
                                            chosenPartyMeasurements{z}{c},...
                                            povms{4}{w}{d});
                                        
                            elseif partyidx == 4
                                term = Tensor(povms{1}{x}{a}, ...
                                            povms{2}{y}{b},...
                                            povms{3}{z}{c},...
                                            chosenPartyMeasurements{w}{d});
                            else
                                disp('error');
                            end
                            summ = summ + term*bellcoeffs(x,y,z,w,a,b,c,d);
                                end
                            end
                        end
                    end
                end        
            end
        end
    end
    objective = real(trace(state*summ));
    
    %% Solve the problem
    constraints = [positivityconstraints, povmconstraints];
    sol = optimize(constraints, -objective, ...
                    sdpsettings('solver','mosek', 'verbose', 0));

    if sol.problem ~= 0
       disp(sol);
       warning('Check what problem there is.'); 
    end
    %% Return output
    newPovms = povms;
    for x=1:ins(partyidx)
        for a=1:outs(partyidx)
            newPovms{partyidx}{x}{a} = value(chosenPartyMeasurements{x}{a});
        end
    end
    finalAlpha = value(objective);
    problemStatus = sol.problem;
    
    %list = whos;for i = 1:length(list);if strcmp(list(i).class,'sdpvar');clear(list(i).name);end;end
    %yalmip("clear");
end