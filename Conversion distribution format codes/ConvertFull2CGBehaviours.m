
function MpCG = ConvertFull2CGBehaviours(pMatrix,Inputs,Outputs)

%Converts collection of behaviours 'pMatrix' to Collin-Gisin notation, in multipartite Bell scenario characterized by 'Inputs' and 'Outputs'


%Inputs:
%'pMatrix' is a N x prod(Inputs)*prod(Outputs) vector, each line of which contains a behaviour, that is, the probabilities p(abc..|xyz...), that is, pMatrix(k,:) = [ p(000..|000..) p(00..1|000..) ... p(oa ob oc ..|000..) p(000..|00..1) ... p(oa ob oc..| nx ny nz ..) ]
%'Inputs' is a 1 x length(Inputs) vector of integers, containing the respective numbers of inputs of each party. That is, Inputs = [nx ny nz ... ]
%'Outputs' is a 1 x length(Inputs) vector of integers, containing the respective numbers of outputs of each party. That is, Outputs = [oa ob oc ... ]

%If 'Inputs'/'Outputs' is a scalar, attributes the same number of inputs/outputs to all parties


%Outputs:

%'MpCG' is a N x prod(Outputs-1)*prod(Inputs+1) vector, each line of which represents a
%behaviour in Collin-Gisin notation, with convention p( 1:oa-1 1:ob-1 .. |
%1:nx+1 1:ny+1 ..), where the last input corresponds to the identity
%measurement (meaning, tracing out the corresponding party)


MpCG = zeros(size(pMatrix,1),prod(Outputs-1)*prod(Inputs+1));


for k = 1 : size(pMatrix,1) 
    
    MpCG(k,:) = ConvertFull2CG(pMatrix(k,:),Inputs,Outputs);
    
end
    


end
