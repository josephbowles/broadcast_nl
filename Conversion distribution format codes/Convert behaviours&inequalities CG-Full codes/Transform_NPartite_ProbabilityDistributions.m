function POutput = Transform_NPartite_ProbabilityDistributions(Type,Scenario,PTarget,Normalization)
% This code can transform multipartite probability distributions in full
% form which is a tensor with all probabilities p(a,b,x,y), or p(a,b,c,x,y,z)
% into the probabilities on the CG (Collins-Gisin, https://arxiv.org/abs/quant-ph/0306129)form (which is more compact)
% Also, if the probabilities are in CG form, this code transforms it into
% the full probability form

% 'InputType' should be either 'CG' or 'full'
% 'Scenario' shoud in the form of [Oa|x=0, Oa|x=1, Ob|x=2, ...
%                                  Ob|y=0, Ob|y=1, Ob|y=2, ...
%                                  Oc|z=0, Oc|z=1, Oc|z=2  ...]
% 'InputTarget' for the 'CG' form should be 
% P(a-1,b-1,c-1,...|x+1,y+1,z+1,...)
% 'InputTarget' for the 'full' form should be in the form of probabilities
% P(a,b,c,...|x,y,z,...)
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
    POutput = -9*ones([Max_NumOutcomes_EachParty-1 Num_Inputs_EachParty+1]);
    CG_abc_FillingIndex = cell(1,Num_Parties); CG_xyz_FillingIndex = cell(1,Num_Parties);
    PTarget_abc_FillingIndex = cell(1,Num_Parties); PTarget_xyz_FillingIndex = cell(1,Num_Parties);
    % Full probabilities part
    for i =1:prod(Num_Inputs_EachParty)
        [PTarget_xyz_FillingIndex{:}] = ind2sub(Num_Inputs_EachParty,i);
        for j=1:Num_Parties
            PTarget_abc_FillingIndex{j} = 1:Scenario(j,PTarget_xyz_FillingIndex{j})-1;            
        end
        POutput(PTarget_abc_FillingIndex{:},PTarget_xyz_FillingIndex{:}) = PTarget(PTarget_abc_FillingIndex{:},PTarget_xyz_FillingIndex{:});        
    end
    % Marginal probabilities part, of course, the nonsignaling constraints are
    % assumed to hold true from one party to another.
    for k_PartiesMarginal=1:Num_Parties-1
        AllPossibleRelatedParties_for_MarginalTerms = nchoosek(1:Num_Parties,k_PartiesMarginal);
        for ith_PartiesCombinationsOfMarginalTerms =1:size(AllPossibleRelatedParties_for_MarginalTerms,1)
            RelatedParties_for_MarginalTerms = AllPossibleRelatedParties_for_MarginalTerms(ith_PartiesCombinationsOfMarginalTerms,:);
            UnRelatedParties_for_MarginalTerms = setdiff(1:Num_Parties,RelatedParties_for_MarginalTerms);
            for j =1:length(UnRelatedParties_for_MarginalTerms)
                CG_abc_FillingIndex{UnRelatedParties_for_MarginalTerms(j)} = 1;
                CG_xyz_FillingIndex{UnRelatedParties_for_MarginalTerms(j)} = Num_Inputs_EachParty(UnRelatedParties_for_MarginalTerms(j))+1;
                PTarget_abc_FillingIndex{UnRelatedParties_for_MarginalTerms(j)} = ':';
                PTarget_xyz_FillingIndex{UnRelatedParties_for_MarginalTerms(j)} = 1;
            end
            Inptus_of_RelatedParties = Num_Inputs_EachParty(RelatedParties_for_MarginalTerms);
            Outputs_of_RelatedParties = Scenario(RelatedParties_for_MarginalTerms,:);
            TmpOutputs_of_RelatedParties_GivenInputs = zeros(1,length(RelatedParties_for_MarginalTerms));
            for m_choice_of_inputs = 1:prod(Inptus_of_RelatedParties)
                [CG_xyz_FillingIndex{[RelatedParties_for_MarginalTerms]}] = ind2sub(Inptus_of_RelatedParties,m_choice_of_inputs);
                [PTarget_xyz_FillingIndex{[RelatedParties_for_MarginalTerms]}] = ind2sub(Inptus_of_RelatedParties,m_choice_of_inputs);
                for i = 1:length(RelatedParties_for_MarginalTerms)
                    TmpOutputs_of_RelatedParties_GivenInputs(i) = Outputs_of_RelatedParties(i,CG_xyz_FillingIndex{RelatedParties_for_MarginalTerms(i)});
                end
                for z = 1:prod(TmpOutputs_of_RelatedParties_GivenInputs-1)
                    [CG_abc_FillingIndex{[RelatedParties_for_MarginalTerms]}] = ind2sub(TmpOutputs_of_RelatedParties_GivenInputs-1,z);
                    [PTarget_abc_FillingIndex{[RelatedParties_for_MarginalTerms]}] = ind2sub(TmpOutputs_of_RelatedParties_GivenInputs-1,z);
                    POutput(CG_abc_FillingIndex{:},CG_xyz_FillingIndex{:}) = sum(vec(PTarget(PTarget_abc_FillingIndex{:},PTarget_xyz_FillingIndex{:})));
                end
            end
        end
    end
    for i =1:Num_Parties
        CG_abc_FillingIndex{i} = 1;
        CG_xyz_FillingIndex{i} = Num_Inputs_EachParty(i)+1;
    end
    POutput(find(POutput<0))=0;
    POutput(CG_abc_FillingIndex{:},CG_xyz_FillingIndex{:})=Normalization;
elseif isequal(Type,'CG')
    assert(isequal(size(PTarget),[Max_NumOutcomes_EachParty-1,Num_Inputs_EachParty+1]),'Dimension of the given target mismatches with the specified Bell scenario.')    
    POutput = zeros([Max_NumOutcomes_EachParty,Num_Inputs_EachParty+1]);
    Output_abc_FillingIndex = cell(1,Num_Parties); Output_xyz_FillingIndex = cell(1,Num_Parties);
    Input_abc_FillingIndex = cell(1,Num_Parties);  Input_xyz_FillingIndex = cell(1,Num_Parties);
    Subtracted_Input_abc_FillingIndex = cell(1,Num_Parties);
    % Copying data from the given input
    for i =1:prod(Num_Inputs_EachParty+1)
        %         [Output_xyz_FillingIndex{:}] = ind2sub(Num_Inputs_EachParty,i);
        [Input_xyz_FillingIndex{:}] = ind2sub(Num_Inputs_EachParty+1,i);
        for j =1:Num_Parties
            if Input_xyz_FillingIndex{j}~=Num_Inputs_EachParty(j)+1
                Input_abc_FillingIndex{j} = 1:Scenario(j,Input_xyz_FillingIndex{j})-1;
            else
                Input_abc_FillingIndex{j} = 1;
            end
        end
        POutput(Input_abc_FillingIndex{:},Input_xyz_FillingIndex{:}) = PTarget(Input_abc_FillingIndex{:},Input_xyz_FillingIndex{:});
    end
    % Marginal probabilities part, of course, the nonsignaling constraints are
    % assumed to hold true from one party to another.
    for k_PartiesMarginal=1:Num_Parties
        AllPossibleRelatedParties_for_MarginalTerms = nchoosek(1:Num_Parties,k_PartiesMarginal);
        for ith_PartiesCombinationsOfMarginalTerms =1:size(AllPossibleRelatedParties_for_MarginalTerms,1)
            RelatedParties_for_MarginalTerms = AllPossibleRelatedParties_for_MarginalTerms(ith_PartiesCombinationsOfMarginalTerms,:);
            UnRelatedParties_for_MarginalTerms = setdiff(1:Num_Parties,RelatedParties_for_MarginalTerms);
            % Fix the inputs and outputs of the parties which are not
            % related to the marginal terms under consideration
            for j =1:length(UnRelatedParties_for_MarginalTerms)
                Input_abc_FillingIndex{UnRelatedParties_for_MarginalTerms(j)} = 1;
                Input_xyz_FillingIndex{UnRelatedParties_for_MarginalTerms(j)} = Num_Inputs_EachParty(UnRelatedParties_for_MarginalTerms(j))+1;
                Output_abc_FillingIndex{UnRelatedParties_for_MarginalTerms(j)} = 1;
                Output_xyz_FillingIndex{UnRelatedParties_for_MarginalTerms(j)} = Num_Inputs_EachParty(UnRelatedParties_for_MarginalTerms(j))+1;
                Subtracted_Input_abc_FillingIndex{UnRelatedParties_for_MarginalTerms(j)}= 1;
            end
            Inptus_of_RelatedParties = Num_Inputs_EachParty(RelatedParties_for_MarginalTerms);
            Outputs_of_RelatedParties = Scenario(RelatedParties_for_MarginalTerms,:);
            for m_choice_of_inputs = 1:prod(Inptus_of_RelatedParties)
                [Output_xyz_FillingIndex{[RelatedParties_for_MarginalTerms]}] = ind2sub(Inptus_of_RelatedParties,m_choice_of_inputs);
                for i = 1:length(RelatedParties_for_MarginalTerms)
                    Parties_in_PreviousLayer = setdiff(RelatedParties_for_MarginalTerms,RelatedParties_for_MarginalTerms(i));
                    Output_abc_FillingIndex{RelatedParties_for_MarginalTerms(i)} = Outputs_of_RelatedParties(i,Output_xyz_FillingIndex{RelatedParties_for_MarginalTerms(i)});
                    Input_abc_FillingIndex{RelatedParties_for_MarginalTerms(i)} = 1;
                    [Input_xyz_FillingIndex{[RelatedParties_for_MarginalTerms]}] = ind2sub(Inptus_of_RelatedParties,m_choice_of_inputs);
                    [Input_xyz_FillingIndex{RelatedParties_for_MarginalTerms(i)}] = Num_Inputs_EachParty(RelatedParties_for_MarginalTerms(i))+1;
                    Subtracted_Input_abc_FillingIndex{RelatedParties_for_MarginalTerms(i)} = ':';
                    Outputs_of_Parties_in_PreviousLayer_GivenInputs = zeros(1,length(Parties_in_PreviousLayer));
                    if isempty(Parties_in_PreviousLayer)
                        POutput(Output_abc_FillingIndex{:},Output_xyz_FillingIndex{:}) = POutput(Input_abc_FillingIndex{:},Input_xyz_FillingIndex{:}) -...
                            sum(vec(POutput(Subtracted_Input_abc_FillingIndex{:},Output_xyz_FillingIndex{:})));
                    else
                        for z = 1:length(Parties_in_PreviousLayer)
                            Tmp = Scenario(Parties_in_PreviousLayer(z),:);
                            if Parties_in_PreviousLayer(z)>RelatedParties_for_MarginalTerms(i)
                                Outputs_of_Parties_in_PreviousLayer_GivenInputs(z) = Tmp(1,Output_xyz_FillingIndex{Parties_in_PreviousLayer(z)})-1;
                            else
                                Outputs_of_Parties_in_PreviousLayer_GivenInputs(z) = Tmp(1,Output_xyz_FillingIndex{Parties_in_PreviousLayer(z)});
                            end
                        end
                        for z = 1:prod(Outputs_of_Parties_in_PreviousLayer_GivenInputs)
                            [Output_abc_FillingIndex{[Parties_in_PreviousLayer]}] = ind2sub(Outputs_of_Parties_in_PreviousLayer_GivenInputs,z);
                            [Subtracted_Input_abc_FillingIndex{[Parties_in_PreviousLayer]}] = ind2sub(Outputs_of_Parties_in_PreviousLayer_GivenInputs,z);
                            [Input_abc_FillingIndex{[Parties_in_PreviousLayer]}] = ind2sub(Outputs_of_Parties_in_PreviousLayer_GivenInputs,z);
                            POutput(Output_abc_FillingIndex{:},Output_xyz_FillingIndex{:}) = POutput(Input_abc_FillingIndex{:},Input_xyz_FillingIndex{:}) -...
                                sum(vec(POutput(Subtracted_Input_abc_FillingIndex{:},Output_xyz_FillingIndex{:})));
                        end
                    end
                end
            end
        end
    end
    for i =1:Num_Parties
       Output_abc_FillingIndex{i} = 1:Max_NumOutcomes_EachParty(i);
       Output_xyz_FillingIndex{i} = 1:Num_Inputs_EachParty(i);
    end
    POutput = POutput(Output_abc_FillingIndex{:},Output_xyz_FillingIndex{:});
end
end
