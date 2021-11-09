function f=local_polytope(n_parties,inputs,outputs)

%Outputs the local vertices in the standard Bell scenario with 'n_parties' parties, and respective numbers of inputs 'inputs', respective numbers of outputs 'outputs' 

%'n_parties' is an integer, representing the number of distant parties
%'inputs' is a 1 x n_parties vector of integers, containing the respective numbers of inputs of each party. That is, inputs = [nx ny nz ... ]
%'outputs' is a 1 x n_parties vector of integers, containing the respective numbers of outputs of each party. That is, outputs = [oa ob oc ... ]

%If 'inputs'/'outputs' are scalars, attribute the same number of inputs/outputs to all parties

%Notation of the output: each line of 'f' is a vector, interpreted as p = [ p(000|000) p(001|000) ... p(oa ob oc |000) p(000|001) ... p(oa ob oc| nx ny nz) ]


%required subfunctions: Matgen_Dax, loop_vector

if length(inputs) == 1
    
    iall=inputs;
    inputs=iall*ones(1,n_parties);
    
end

if length(outputs) == 1
    
    oall=outputs;
    outputs=oall*ones(1,n_parties);
    
end

%number of players

N=n_parties;


%dimension of the probability space

d=prod(outputs)*prod(inputs);


%number of deterministic strategies

L = prod( outputs .^ inputs);


%individual strategies

Ptot=cell(N,1);

for k=1:N
    
    Ptot{k}=Matgen_Dax(outputs(k),inputs(k));
    
end


%Global strategies

P=zeros(L,d);

Idx=loop_vector( outputs .^ inputs) ;

for k=1:L
    
    vec=Idx(k,:);
    
    Paux=Ptot{1};
    V=Paux(:,:,vec(1));
    
    for i=2:N
        
        Paux=Ptot{i};
        Vp=Paux(:,:,vec(i));
        
        V=kron(V,Vp);
        
    end
    
    
    P(k,:)=reshape(V,d,1);
    
end



f=P;
    

end

    
    
