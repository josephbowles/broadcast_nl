function probability_ndarray2 = CalcProbArray(state,povms, inputs_per_party, outputs_per_party)
    % This is used somewhere for the optimization over the 4 party
    % broadcast scenario with 2 alices and 2 bobs.

    nrparties = length(inputs_per_party);
    
    % Define a zero multidim where to place the probability distribution
    dims = num2cell([inputs_per_party,outputs_per_party]);
    probability_ndarray2 = zeros(dims{:});
     
    aux = [inputs_per_party, outputs_per_party];
    allinputoutputcombinations = ind2subv(aux, 1:prod(aux(:)));
    for slice=1:size(allinputoutputcombinations,1)
        ins = num2cell(allinputoutputcombinations(slice,1:nrparties));
        outs = num2cell(allinputoutputcombinations(slice,nrparties+1:end));
        %disp([ins{:}, outs{:}]);
        tensor = 1;
        for p=1:nrparties
             tensor = Tensor(tensor,povms{p}{ins{p}}{outs{p}});
        end

        probability_ndarray2(ins{:},outs{:}) = real(trace(tensor*state));
    end
 

% ins = inputs_per_party;
% outs = outputs_per_party;
% 
% probability_ndarray = zeros([ins,outs]);
% for x=1:ins(1)
%     for y=1:ins(2)
%         for z=1:ins(3)
%             for w=1:ins(4)
%                 for a=1:outs(1)
%                     for b=1:outs(2)
%                         for c=1:outs(3)
%                             for d=1:outs(4)
%                                     probability_ndarray(x,y,z,w,a,b,c,d) = real(trace(state*Tensor(povms{1}{x}{a}, ...
%                                                                          povms{2}{y}{b}, ...
%                                                                          povms{3}{z}{c}, ...
%                                                                          povms{4}{w}{d})));
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end




end