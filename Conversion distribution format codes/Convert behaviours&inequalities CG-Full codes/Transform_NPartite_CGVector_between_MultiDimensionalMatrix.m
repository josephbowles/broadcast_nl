function POutput = Transform_NPartite_CGVector_between_MultiDimensionalMatrix(Type,Scenario,InputTarget)
Num_Parties = size(Scenario,1);
Num_Inputs_EachParty = zeros(1,Num_Parties);
for i =1:Num_Parties
    Num_Inputs_EachParty(i) = length(find(Scenario(i,:)~=0));
end
Max_NumOutcomes_EachParty = max(Scenario');
if isequal(Type,'CGV')
    Num_Parties = size(Scenario,1);
    Num_Inputs_EachParty = zeros(1,Num_Parties);
    for i =1:Num_Parties
        Num_Inputs_EachParty(i) = length(find(Scenario(i,:)~=0));
    end
    Max_NumOutcomes_EachParty = max(Scenario');
    POutput = zeros([Max_NumOutcomes_EachParty-1,Num_Inputs_EachParty]);
    Output_abc_FillingIndex = cell(1,Num_Parties); Output_xyz_FillingIndex = cell(1,Num_Parties);
    CountedIndex = 1;
    for i = 1:prod(Num_Inputs_EachParty+1)
        [Output_xyz_FillingIndex{:}] = ind2sub(Num_Inputs_EachParty+1,i);
        Outputs_of_Parties_GivenInputs = zeros(1,Num_Parties);
        for j =1:Num_Parties
            if Output_xyz_FillingIndex{j} ~= Num_Inputs_EachParty(j)+1
                Outputs_of_Parties_GivenInputs(j) = Scenario(j,Output_xyz_FillingIndex{j})-1;
            else
                Outputs_of_Parties_GivenInputs(j) = 1;
            end
        end
        for k =1:prod(Outputs_of_Parties_GivenInputs)
            [Output_abc_FillingIndex{:}] = ind2sub(Outputs_of_Parties_GivenInputs,k);
            POutput(Output_abc_FillingIndex{:},Output_xyz_FillingIndex{:}) = InputTarget(CountedIndex);
            CountedIndex = CountedIndex + 1;
        end
    end
    assert(isequal(CountedIndex-1,length(vec(InputTarget))),'Dimension of the input CGVector mismatches with the specified scenario');
elseif isequal(Type,'CGM')
    Input_abc_FillingIndex = cell(1,Num_Parties);
    Input_xyz_FillingIndex = cell(1,Num_Parties);
    Num_output_given_input = zeros(1,Num_Parties);
    POutput =0;
    CountedIndex = 1; 
    for ith_inputs =1:prod(Num_Inputs_EachParty+1);
        [Input_xyz_FillingIndex{:}] = ind2sub(Num_Inputs_EachParty+1,ith_inputs);
        for j =1:Num_Parties
            if Input_xyz_FillingIndex{j} ~= Num_Inputs_EachParty(j)+1
                Num_output_given_input(j) = Scenario(j,Input_xyz_FillingIndex{j})-1;
            else
                Num_output_given_input(j) = 1;
            end
        end
        for kth_output_given_input = 1:prod(Num_output_given_input)
            [Input_abc_FillingIndex{:}] = ind2sub(Num_output_given_input,kth_output_given_input);            
            POutput(CountedIndex) = InputTarget(Input_abc_FillingIndex{:},Input_xyz_FillingIndex{:});
            CountedIndex = CountedIndex+1;
        end
    end
    POutput = POutput';
end
end