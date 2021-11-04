function finaloutput= FromProjectorNotationToProbability(correlToProj, ins, outs)
max_nr_of_parties = length(ins);
[C,T] = coeffs(correlToProj);
nrterms = length(T);
probs = cell(nrterms,1);
summ = 0;
for i=1:nrterms
    auxout = GivePartySettingsOutputsFromProjectorProducts(T(i));
    if ~isequal(1:max_nr_of_parties,auxout(:,1)')
        disp(1:max_nr_of_parties)
        disp(auxout(:,1))
        error('the factors should be ordered!');
    end
    settings = auxout(:,2);
    outputs = auxout(:,3);
    probsym = GiveProbSymVarFromInputsOutputs(settings,outputs);
    probs{i} = probsym;
    summ = summ + C(i) * probsym;
end
finaloutput = expand(summ);
end

function out = GivePartyIdxAndInputFromSingleProjector(symvar)
% example: input Pi_2_3_1
% output: [2,3,1]
name = char(symvar); % ex: sym var A0 --> 'A0'
partyIdx = str2double(name(4));%-'A'+1;
partyInput = str2double(name(6));
partyOutput = str2double(name(8)); % assuming single digits
out(1) = partyIdx;
out(2) = partyInput;
out(3) = partyOutput;
end

function out = GivePartySettingsOutputsFromProjectorProducts(symvar)
theFactors = factor(symvar);
nrFactors = length(theFactors);
out = zeros(nrFactors,3);
for i=1:nrFactors
    auxoutput = GivePartyIdxAndInputFromSingleProjector(theFactors(i));
    out(i,1) = auxoutput(1);%partyIdx;
    out(i,2) = auxoutput(2);%partyInput;
    out(i,3) = auxoutput(3);%partyOutput;
end
end

function symvar = GiveProbSymVarFromInputsOutputs(settings,outputs)
str = 'p_';
for i_out = 1:length(outputs)
   str = strcat(str, string( outputs(i_out) ));  
end
str = strcat(str,'_');
for i_in = 1:length(settings)
   str = strcat(str, string( settings(i_in) )); 
end
symvar = sym(str);
end