function finaloutput = ToProjectorNotationTerm(symcorrelator, ins, outs, FLAG_Use01obsInsteadOfCorrelator)
if isa(symcorrelator,'double')
   finaloutput = symcorrelator;
   return;
end
try  
    % if the input is a sym number instead of just a double then converting to double will work and we return this
    % else we catch this error and proceed with the code as intended
    finaloutput = double(symcorrelator);
    return;
catch 
% if double(symcorrelator) gave an error then it is a monomial and we
% proceed as follows
max_nr_of_parties = length(ins);

factorarray = factor(symcorrelator);
lengthfactorarray = length(factorarray);

factorsInIndexFormat = zeros(lengthfactorarray,2); % row is the factor, first column 
                                      % is the party (1->'A', 2->'B', etc.)
                                      % and the second column is the 
                                      % input
                                     
for i=1:lengthfactorarray
    auxoutput = GivePartyIdxAndInputFromCorrSymName(factorarray(i)); 
    factorsInIndexFormat(i,1) = auxoutput(1); % partyidx
    factorsInIndexFormat(i,2) = auxoutput(2); % partyinput
end
% sort in increasing order so 'A' is first, 'B' second and 'C' third, etc.
factorsInIndexFormat = sortrows(factorsInIndexFormat); % this sorts by the first row in increasing order


% now we need to identify the missing parties
% we assume that all parties are indexed from 1 to max_nr_of_parties
missing_parties = setdiff(1:max_nr_of_parties,factorsInIndexFormat(:,1)');
present_parties = factorsInIndexFormat(:,1);


correlToProj = 1;
for partyidx = 1:max_nr_of_parties
    if ismember(partyidx, present_parties)
        IndexForPartyIdx = dot(factorsInIndexFormat(:,1)==partyidx, ...
                               1:size(factorsInIndexFormat,1));
        partyinput = factorsInIndexFormat(IndexForPartyIdx, 2);
        if FLAG_Use01obsInsteadOfCorrelator
            % if outputs are 0,1 then in papers such as https://arxiv.org/pdf/1006.3032.pdf
            % they don't use correlators but Ai = 0*\Pi_1 + 1*\Pi_2
            correlToProj = correlToProj * ( sym(strcat('Pi_',string(partyidx), ...
                                                 '_',string(partyinput), ...
                                                 '_','1')) );
        else
            correlToProj = correlToProj * ( sym(strcat('Pi_',string(partyidx), ...
                                     '_',string(partyinput), ...
                                     '_','1')) ...
                                     - ...
                                     sym(strcat('Pi_',string(partyidx), ...
                                                 '_',string(partyinput), ...
                                                 '_','2')));
        end

    elseif ismember(partyidx, missing_parties)
        partyinput = 1;
        %partyinput = randi(ins(partyidx));
        %sum over outputs
        sumoveroutputs = 0;
        for output=1:outs(partyidx)
            sumoveroutputs = sumoveroutputs + sym(strcat('Pi_',string(partyidx), ...
                                                           '_',string(partyinput), ...
                                                           '_',string(output)));
        end
        correlToProj = correlToProj * sumoveroutputs;
    end
end
finaloutput = expand(correlToProj);
end
end

function out = GivePartyIdxAndInputFromCorrSymName(symvar)
% input something like A0 or B1 or C0
% example: input B1
% output: [2,1], substitute A->1, B->2, C->3
out=zeros(1,2);
name = char(symvar); % ex: sym var A0 --> 'A0'
partyIdx = name(1)-'A'+1;
partyInput = str2double(name(2));
out(1) = partyIdx;
out(2) = partyInput;
end

