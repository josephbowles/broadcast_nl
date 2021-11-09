function OutputP = Transform_NPartite_Correlators(Type,Scenario,PTarget,Normalization)
% This code can transform multipartite probability distributions in CG form/
% probabilities into probabilites/ CG form.
%
% 'InputType' should be either 'CG' or 'full'
% 'Scenario' shoud in the form of [Oa|x=0, Oa|x=1, Oa|x=2, ...
%                                  Ob|y=0, Ob|y=1, Ob|y=2, ...
%                                  Oc|z=0, Oc|z=1, Oc|z=2  ...]
% 'InputTarget' for the 'CG' form should be
% P(a-1,b-1,c-1,...|x+1,y+1,z+1,...)
% 'InputTarget' for the 'full' form should be in the form of probabilities
% P(a,b,c,...|x,y,z,...)
assert(isequal(unique(Scenario(find(vec(Scenario)))),2))
if nargin ==3
    Normalization = 1;
end
Num_Parties = size(Scenario,1);
Num_Inputs_EachParty = zeros(1,Num_Parties);
for i =1:Num_Parties
    Num_Inputs_EachParty(i) = length(find(Scenario(i,:)~=0));
end
Max_NumOutcomes_EachParty = max(Scenario');
if isequal(Type,'full')
    assert(isequal(size(PTarget),[Max_NumOutcomes_EachParty,Num_Inputs_EachParty]),'Dimension of the given target mismatches with the specified Bell scenario.')
    OutputP = -9*ones(Num_Inputs_EachParty+1);
    E_xyz_FillingIndex = cell(1,Num_Parties);
    PTarget_abc_FillingIndex = cell(1,Num_Parties); PTarget_xyz_FillingIndex = cell(1,Num_Parties);
    % Full probabilities part
    for i =1:prod(Num_Inputs_EachParty)
        [PTarget_xyz_FillingIndex{:}] = ind2sub(Num_Inputs_EachParty,i);
        TmpValue = 0;
        for j=1:2^Num_Parties
            [PTarget_abc_FillingIndex{:}] = ind2sub(2*ones(1,Num_Parties),j);
            TmpValue = TmpValue + (-1)^sum(cell2mat(PTarget_abc_FillingIndex)-1)*PTarget(PTarget_abc_FillingIndex{:},PTarget_xyz_FillingIndex{:});
        end
        OutputP(PTarget_xyz_FillingIndex{:}) = TmpValue;
    end
    % Marginal probabilities part, of course, the nonsignaling constraints are
    % assumed to hold true from one party to another.
    for k_PartiesMarginal=1:Num_Parties-1
        AllPossibleRelatedParties_for_MarginalTerms = nchoosek(1:Num_Parties,k_PartiesMarginal);
        for ith_PartiesCombinationsOfMarginalTerms =1:size(AllPossibleRelatedParties_for_MarginalTerms,1)
            RelatedParties_for_MarginalTerms = AllPossibleRelatedParties_for_MarginalTerms(ith_PartiesCombinationsOfMarginalTerms,:);
            UnRelatedParties_for_MarginalTerms = setdiff(1:Num_Parties,RelatedParties_for_MarginalTerms);
            for j =1:length(UnRelatedParties_for_MarginalTerms)
                E_xyz_FillingIndex{UnRelatedParties_for_MarginalTerms(j)} = Num_Inputs_EachParty(UnRelatedParties_for_MarginalTerms(j))+1;
                PTarget_abc_FillingIndex{UnRelatedParties_for_MarginalTerms(j)} = ':';
                PTarget_xyz_FillingIndex{UnRelatedParties_for_MarginalTerms(j)} = 1;
            end
            Inptus_of_RelatedParties = Num_Inputs_EachParty(RelatedParties_for_MarginalTerms);
            Outputs_of_RelatedParties = 2*ones(1,length(RelatedParties_for_MarginalTerms));
            for m_choice_of_inputs = 1:prod(Inptus_of_RelatedParties)
                [E_xyz_FillingIndex{[RelatedParties_for_MarginalTerms]}] = ind2sub(Inptus_of_RelatedParties,m_choice_of_inputs);
                [PTarget_xyz_FillingIndex{[RelatedParties_for_MarginalTerms]}] = ind2sub(Inptus_of_RelatedParties,m_choice_of_inputs);
                TmpValue = 0;
                for z = 1:prod(Outputs_of_RelatedParties)
                    [PTarget_abc_FillingIndex{[RelatedParties_for_MarginalTerms]}] = ind2sub(Outputs_of_RelatedParties,z);
                    TmpValue = TmpValue + (-1)^(sum([PTarget_abc_FillingIndex{[RelatedParties_for_MarginalTerms]}]-1))*sum(vec(PTarget(PTarget_abc_FillingIndex{:},PTarget_xyz_FillingIndex{:})));
                end
                OutputP(E_xyz_FillingIndex{:}) = TmpValue;
            end
        end
    end
    OutputP(end)=Normalization;
elseif isequal(Type,'Corr')
    Output_abc_FillingIndex = cell(1,Num_Parties); Output_xyz_FillingIndex = cell(1,Num_Parties);
    Input_xyz_FillingIndex = cell(1,Num_Parties);
    OutputP = zeros([Max_NumOutcomes_EachParty,Num_Inputs_EachParty]);
    Default_xyz_FillingIndex = cell(1,Num_Parties); Default_abc_FillingIndex = cell(1,Num_Parties);
    for i =1:Num_Parties
        Default_xyz_FillingIndex{i} = Num_Inputs_EachParty(i)+1;
        Default_abc_FillingIndex{i} = 1;
    end
    for zth_input = 1:prod(Num_Inputs_EachParty)
        [Input_xyz_FillingIndex{:}] = ind2sub(Num_Inputs_EachParty,zth_input);
        [Output_xyz_FillingIndex{:}] = ind2sub(Num_Inputs_EachParty,zth_input);
        for jth_output = 1:prod(ones(1,Num_Parties)*2)
            [Output_abc_FillingIndex{:}] = ind2sub(ones(1,Num_Parties)*2,jth_output);
%             Output_abc_FillingIndex
%             Output_xyz_FillingIndex
            TmpValue = 1;
            for k_PartiesMarginal=1:Num_Parties
                AllPossibleRelatedParties_for_MarginalTerms = nchoosek(1:Num_Parties,k_PartiesMarginal);
                for ith_PartiesCombinationsOfMarginalTerms =1:size(AllPossibleRelatedParties_for_MarginalTerms,1)
                    TmpInput_xyz_FillingIndex= Input_xyz_FillingIndex;
                    TmpOutput_abc_FillingIndex= Output_abc_FillingIndex;
                    RelatedParties_for_MarginalTerms = AllPossibleRelatedParties_for_MarginalTerms(ith_PartiesCombinationsOfMarginalTerms,:);
                    UnRelatedParties_for_MarginalTerms = setdiff(1:Num_Parties,RelatedParties_for_MarginalTerms);
                    if ~isempty(UnRelatedParties_for_MarginalTerms)
                        for kkk = 1:length(UnRelatedParties_for_MarginalTerms)
                        TmpInput_xyz_FillingIndex{[UnRelatedParties_for_MarginalTerms(kkk)]} = Default_xyz_FillingIndex{[UnRelatedParties_for_MarginalTerms(kkk)]};
                        TmpOutput_abc_FillingIndex{[UnRelatedParties_for_MarginalTerms(kkk)]} = Default_abc_FillingIndex{[UnRelatedParties_for_MarginalTerms(kkk)]};
                        end
                    end
%                     TmpInput_xyz_FillingIndex
%                     (-1)^sum(cell2mat(TmpOutput_abc_FillingIndex)-1)
                    TmpValue = TmpValue + (-1)^sum(cell2mat(TmpOutput_abc_FillingIndex)-1)*PTarget(TmpInput_xyz_FillingIndex{:});
                end
            end
            OutputP(Output_abc_FillingIndex{:},Output_xyz_FillingIndex{:}) = 1/2^(Num_Parties)*TmpValue;
        end
        
    end
    
    
    
end
end
