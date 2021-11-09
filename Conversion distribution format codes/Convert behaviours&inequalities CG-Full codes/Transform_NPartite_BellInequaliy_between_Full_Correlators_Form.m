function OutputBellIneq = Transform_NPartite_BellInequaliy_between_Full_Correlators_Form(Type,Scenario,InputBellIneq)

%Transforms a Bell inequality between correlator notation and full notation

%'InputBellineq' should be given in table form


%'Type' specifies the type of 'InputBellIneq, should be either 'full' or
%'Corr', the output will be in the other type 

%'Scenario' is a 2-dimensional matrix of the form:
%[Oa|x=1, Oa|x=2, Oa|x=3, …, Oa|x=N_x,
% Ob|y=1, Ob|y=2, Ob|y=3, …, Ob|y=N_y,
% Oc|z=1, Oc|z=2, Oc|z=3  ...]


assert(isequal(unique(Scenario(find(vec(Scenario)))),2))
Num_Parties = size(Scenario,1);
Num_Inputs_EachParty = zeros(1,Num_Parties);
for i =1:Num_Parties
    Num_Inputs_EachParty(i) = length(find(Scenario(i,:)~=0));
end
Max_NumOutcomes_EachParty = max(Scenario');
Output_abc_FillingIndex = cell(1,Num_Parties); Output_xyz_FillingIndex = cell(1,Num_Parties);
Input_abc_FillingIndex = cell(1,Num_Parties);  Input_xyz_FillingIndex = cell(1,Num_Parties);
if isequal(Type,'full')
    OutputBellIneq = zeros(Num_Inputs_EachParty+1);
    % Marginal probabilities part, of course, the nonsignaling constraints are
    % assumed to hold true from one party to another.
    for k_PartiesMarginal=0:Num_Parties
        AllPossibleRelatedParties_for_MarginalTerms = nchoosek(1:Num_Parties,k_PartiesMarginal);
        for ith_PartiesCombinationsOfMarginalTerms =1:size(AllPossibleRelatedParties_for_MarginalTerms,1)
            RelatedParties_for_MarginalTerms = AllPossibleRelatedParties_for_MarginalTerms(ith_PartiesCombinationsOfMarginalTerms,:);
            UnRelatedParties_for_MarginalTerms = setdiff(1:Num_Parties,RelatedParties_for_MarginalTerms);
            for j =1:length(UnRelatedParties_for_MarginalTerms)
                Output_xyz_FillingIndex{UnRelatedParties_for_MarginalTerms(j)} = Num_Inputs_EachParty(UnRelatedParties_for_MarginalTerms(j))+1;
                Input_xyz_FillingIndex{UnRelatedParties_for_MarginalTerms(j)} = 1:Num_Inputs_EachParty(UnRelatedParties_for_MarginalTerms(j));
                Input_abc_FillingIndex{UnRelatedParties_for_MarginalTerms(j)} = ':';
            end
            NumOutputs_of_RelatedParties = Scenario(RelatedParties_for_MarginalTerms,:);
            
            Num_Inputs_of_RelatedParties = Num_Inputs_EachParty(RelatedParties_for_MarginalTerms);
            TmpNumOutputs_of_RelatedParties_GivenInputs = zeros(1,length(RelatedParties_for_MarginalTerms));
            for m_choice_of_inputs = 1:prod(Num_Inputs_of_RelatedParties)
                [Input_xyz_FillingIndex{[RelatedParties_for_MarginalTerms]}] = ind2sub(Num_Inputs_of_RelatedParties,m_choice_of_inputs);
                [Output_xyz_FillingIndex{[RelatedParties_for_MarginalTerms]}] = ind2sub(Num_Inputs_of_RelatedParties,m_choice_of_inputs);
                for i = 1:length(RelatedParties_for_MarginalTerms)
                    TmpNumOutputs_of_RelatedParties_GivenInputs(i) = NumOutputs_of_RelatedParties(i,Output_xyz_FillingIndex{RelatedParties_for_MarginalTerms(i)});
                end
                TmpValue = 0;
                for z_output_of_RelatedParties = 1:prod(TmpNumOutputs_of_RelatedParties_GivenInputs)
                    [Input_abc_FillingIndex{[RelatedParties_for_MarginalTerms]}] = ind2sub(TmpNumOutputs_of_RelatedParties_GivenInputs,z_output_of_RelatedParties);
                    T = 0;
                    for i = 1:length(RelatedParties_for_MarginalTerms)
                        T = T +(Input_abc_FillingIndex{RelatedParties_for_MarginalTerms(i)}-1);
                    end
                    TmpPlusMinus =(-1)^(T);
                    TmpValue = TmpValue + TmpPlusMinus*sum(vec(InputBellIneq(Input_abc_FillingIndex{:},Input_xyz_FillingIndex{:})));
                end
                OutputBellIneq(Output_xyz_FillingIndex{:}) = TmpValue;
            end
        end
    end
    OutputBellIneq = OutputBellIneq/2^Num_Parties;
elseif isequal(Type,'Corr')
    OutputBellIneq = zeros([Max_NumOutcomes_EachParty,Num_Inputs_EachParty]);
    for k_PartiesMarginal=0:Num_Parties
        AllPossibleRelatedParties_for_MarginalTerms = nchoosek(1:Num_Parties,k_PartiesMarginal);
        for ith_PartiesCombinationsOfMarginalTerms =1:size(AllPossibleRelatedParties_for_MarginalTerms,1)
            RelatedParties_for_MarginalTerms = AllPossibleRelatedParties_for_MarginalTerms(ith_PartiesCombinationsOfMarginalTerms,:);
            UnRelatedParties_for_MarginalTerms = setdiff(1:Num_Parties,RelatedParties_for_MarginalTerms);
            for j =1:length(UnRelatedParties_for_MarginalTerms)
                Output_xyz_FillingIndex{UnRelatedParties_for_MarginalTerms(j)} = 1 ;
                Output_abc_FillingIndex{UnRelatedParties_for_MarginalTerms(j)} = 1:Scenario(UnRelatedParties_for_MarginalTerms(j),1);
                Input_xyz_FillingIndex{UnRelatedParties_for_MarginalTerms(j)} = Num_Inputs_EachParty(UnRelatedParties_for_MarginalTerms(j))+1;
%                 Input_abc_FillingIndex{UnRelatedParties_for_MarginalTerms(j)} = ':';
            end
            NumOutputs_of_RelatedParties = Scenario(RelatedParties_for_MarginalTerms,:);            
            Num_Inputs_of_RelatedParties = Num_Inputs_EachParty(RelatedParties_for_MarginalTerms);
            TmpNumOutputs_of_RelatedParties_GivenInputs = zeros(1,length(RelatedParties_for_MarginalTerms));
            for m_choice_of_inputs = 1:prod(Num_Inputs_of_RelatedParties)
                [Input_xyz_FillingIndex{[RelatedParties_for_MarginalTerms]}] = ind2sub(Num_Inputs_of_RelatedParties,m_choice_of_inputs);
                [Output_xyz_FillingIndex{[RelatedParties_for_MarginalTerms]}] = ind2sub(Num_Inputs_of_RelatedParties,m_choice_of_inputs);
                for i = 1:length(RelatedParties_for_MarginalTerms)
                    TmpNumOutputs_of_RelatedParties_GivenInputs(i) = NumOutputs_of_RelatedParties(i,Output_xyz_FillingIndex{RelatedParties_for_MarginalTerms(i)});
                end
                for z_output_of_RelatedParties = 1:prod(TmpNumOutputs_of_RelatedParties_GivenInputs)
                    [Output_abc_FillingIndex{[RelatedParties_for_MarginalTerms]}] = ind2sub(TmpNumOutputs_of_RelatedParties_GivenInputs,z_output_of_RelatedParties);
                    T = 0;
                    for i = 1:length(RelatedParties_for_MarginalTerms)
                        T = T +(Output_abc_FillingIndex{RelatedParties_for_MarginalTerms(i)}-1);
                    end
                    TmpPlusMinus =(-1)^(T);
                    OutputBellIneq(Output_abc_FillingIndex{:},Output_xyz_FillingIndex{:}) = OutputBellIneq(Output_abc_FillingIndex{:},Output_xyz_FillingIndex{:})  + TmpPlusMinus*InputBellIneq(Input_xyz_FillingIndex{:});
                end
            end
        end
    end
end
end