function probability_ndarray = ProbMultidimArray(state,povms, inputs_per_party, outputs_per_party)
    % povms is called as povms{party}{input}{output}
% %     nrparties = max(size(povms));
% %     inputs_per_party = zeros(1,nrparties);
% %     outputs_per_party = zeros(1,nrparties);
% %     for p=1:nrparties
% %         inputs_per_party(p) = max(size(povms{p}));
% %         outputs_per_party(p) = max(size(povms{p}{1}));
% %     end
    nrparties = length(inputs_per_party);
    
    % Define a zero multidim where to place the probability distribution
    dims = num2cell([inputs_per_party,outputs_per_party]);
    probability_ndarray = zeros(dims{:});
        
    aux = [inputs_per_party, outputs_per_party];
    allinputoutputcombinations = ind2subv(aux, 1:prod(aux(:)));
    for slice=1:size(allinputoutputcombinations,1)
        ins = num2cell(allinputoutputcombinations(slice,1:nrparties));
        outs = num2cell(allinputoutputcombinations(slice,nrparties+1:end));

        tensor = 1;
        for p=1:nrparties
             tensor = kron(tensor,povms{p}{ins{p}}{outs{p}});
        end

        probability_ndarray(ins{:},outs{:}) = real(trace(tensor*state));
    end
    %probability_ndarray = clean(probability_ndarray, 1e-6);
end
% 
% function p = prob(state,Pproj,coords)
%     nrparties = max(size(coords))/2;
%     assert(mod(max(size(coords)),2)==0,"'coords' must have an even number of components");
%     settings = coords(1:nrparties);
%     outputs = coords(nrparties+1:end);
%     tensor = 1;
%     for p=1:nrparties
%         tensor = kron(tensor,Pproj{p}{settings(p)}{outputs(p)});
%     end
%     p = real(trace(state*tensor));
% end