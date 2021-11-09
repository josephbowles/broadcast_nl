function p = MaximallyMixedBehaviour(N_parties,Inputs,Outputs)

%In scenario where there is 'N_parties', with inputs in 'Inputs', outputs in 'Outputs', gives the maximally mixed behaviour


%'N_parties' is an integer, representing the number of distant parties
%'Inputs' is a 1 x N_parties vector of integers, containing the respective numbers of inputs of each party. That is, Inputs = [nx ny nz ... ]
%'Outputs' is a 1 x N_parties vector of integers, containing the respective numbers of outputs of each party. That is, Outputs = [oa ob oc ... ]

%If 'Inputs'/'Outputs' is a scalar, attributes the same number of inputs/outputs to all parties

%Notation of the output: p = [ p(000|000) p(001|000) ... p(oa ob oc |000) p(000|001) ... p(oa ob oc| nx ny nz) ]


if length(Inputs) == 1
    
    Inputs=Inputs*ones(1,N_parties);
    
end

if length(Outputs) == 1
    
    Outputs=Outputs*ones(1,N_parties);
    
end


%Global inputs and outputs:

NI = prod(Inputs); 

NO = prod(Outputs); 


p = ones(1,NI*NO) * 1/NO;

end




    