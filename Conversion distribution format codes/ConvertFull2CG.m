
function pCG = ConvertFull2CG(pfull,Inputs,Outputs)

%Converts behaviour 'pfull' to Collin-Gisin notation, in multipartite Bell scenario characterized by 'Inputs' and 'Outputs'


%Inputs:
%'pfull' is a 1 x prod(Inputs)*prod(Outputs) vector, containing the probabilities p(abc..|xyz...), that is, pfull = [ p(000..|000..) p(00..1|000..) ... p(oa ob oc ..|000..) p(000..|00..1) ... p(oa ob oc..| nx ny nz ..) ]
%'Inputs' is a 1 x length(Inputs) vector of integers, containing the respective numbers of inputs of each party. That is, Inputs = [nx ny nz ... ]
%'Outputs' is a 1 x length(Inputs) vector of integers, containing the respective numbers of outputs of each party. That is, Outputs = [oa ob oc ... ]

%If 'Inputs'/'Outputs' is a scalar, attributes the same number of inputs/outputs to all parties


%Outputs:

%'pCG' is a 1 x prod(Outputs-1)*prod(Inputs+1) vector, representing the
%behaviour in Collin-Gisin notation, with convention p( 1:oa-1 1:ob-1 .. |
%1:nx+1 1:ny+1 ..), where the last input corresponds to the identity
%measurement (meaning, tracing out the corresponding party)


assert(isequal(size(Inputs),size(Outputs)),'Vectors Inputs and Outputs should be the same size')


%If Inputs and Outputs are scalar, convert to vector 
if length(Inputs) == 1 
    
    NParties = log(length(pfull))/log(Inputs*Outputs);
    
    Inputs = Inputs*ones(1,NParties);
    Outputs = Outputs*ones(1,NParties);  
    
else

    
NParties = length(Inputs); %number of parties

end


%Create 'Scenario' matrix, Scenario = [Oa|x=0, Oa|x=1, Ob|x=2, ...
%                                      Ob|y=0, Ob|y=1, Ob|y=2, ...
%                                      Oc|z=0, Oc|z=1, Oc|z=2  ...]

Scenario = zeros(NParties,max(Inputs));

for k=1:NParties
    
    Scenario_k = Outputs(k)*ones(1,Inputs(k));
    
    Scenario(k,1:length(Scenario_k)) = Scenario_k;
    
end

pT=convP2Table(pfull,[Inputs Outputs]);%transform behaviour to table form
pT = permute(pT,[NParties+1:2*NParties 1:NParties]); %Swap inputs and outputs so that indexing starts with outputs

C = Transform_NPartite_ProbabilityDistributions('full',Scenario,pT); %transform to Collin-Gisin

pCG = convTable2P(permute(C,[NParties+1:2*NParties 1:NParties]),[Inputs Outputs]); %Swap back inputs and outputs and transform to vector



end
