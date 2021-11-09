%Set the parameters for HeuristicSearchBroadcastScript

%This script is meant to set all the variables needed to run a heuristic
%search in the broadcast scenario

%Initial family of states (the family of states is of the form v*rho+(1-v)*rhoNoise)

clear %Clear all previously defined variables to avoid potential conflict
d=2; dims=[d d]; %vector of local dimensions of the initial state
rho = GHZState(d,size(dims,2))*GHZState(d,size(dims,2))'; %entangled part of the state
rhoNoise = eye(prod(dims))/prod(dims); %noisy part of the state

%Final local dimensions (here one sets the desired dimension partition of each final party (i.e. after the Choi maps). That is, dloc{1} = [dA1 dA2 .. ], dloc{2} =  [dB1 dB2 .. ], ...
dloc=cell(1,2); %create the cell that contains vectors of local dimensions
dloc{1} = [d]; %local dimensions of Alice (party 1)
dloc{2} = [d d]; %local dimension of Bob (party 2)
%dloc{3} = [d d]; %local dimension of Charlie (party 3)

%Scenario (here one sets the desired numbers of inputs/outputs of each final party (i.e. after the Choi maps). That is, for party 1, nInputs{1} = [nx1 nx2 ..], for party 2, nInputs{2} =  [ny1 ny2 ..], etc.

% Inputs
nInputs=cell(1,2); %create the cell that contains vectors of local inputs
nInputs{1}=[2]; %local inputs of Alice (party 1)
nInputs{2}=[2 2]; %local inputs of Bob (party 2)
%nInputs{3} = [2 2]; %local inputs of Charlie (party 3)

% Outputs
nOutputs=cell(1,2); %create the cell that contains vectors of local outputs
nOutputs{1}=[2]; %local outputs of Alice (party 1)
nOutputs{2}=[2 2]; %local outputs of Bob (party 2)
%nOutputs{3} = [2 2]; %local outputs of Charlie (party 3)

WChoi = [2]; %Where to apply Choi states
%that is, WChoi contains the label of initial parties one wants to apply a
%Choi map on (eg, WChoi = 2 means a Choi map will be apply to Bob, WChoi =
%[1 2] means a Choi map will be applied to Alice, and a Choi map will be
%applied to Bob as well, etc. )



%Vertices of all but last party
Vertices=cell(1,size(dims,2)-1);

%No-signaling vertices
%Vertices{1} = VNS22;
%Vertices{2} = Bobvertices



%Local vertices for the rest

if isempty(Vertices{1}) == 1
    Vertices{1}=local_polytope(length(nInputs{1}),nInputs{1},nOutputs{1}); %strategies for final parties coming from (initial) party 1
end

outA = prod(nOutputs{1}); %number of overall outputs

%produce a warning if local vertices are used where more than one local party
if length(dloc{1}) > 1 && size(Vertices{1},1) == prod(nOutputs{1}.^nInputs{1}), disp('Warning: there are several Alices but local vertices are currently used (does not correspond to a broadcast-local set); if needed, change Vertices{k} accordingly'), end

VerticesAllbutLast = Vertices{1};
for k=2:length(dims)-1
    
    if isempty(Vertices{k}) == 1
        Vertices{k}=local_polytope(length(nInputs{k}),nInputs{k},nOutputs{k}); %strategies for final parties coming from party k
    end
    
    %produce a warning if local vertices are used where more than one local party
    if length(dloc{k}) > 1 && size(Vertices{k},1) == prod(nOutputs{k}.^nInputs{k}), disp('Warning: there are several local parties in initial party '), disp(k), disp(' but local vertices are currently used (does not correspond to a broadcast-local set); if needed, change Vertices{k} accordingly'), end
    
    VerticesAllbutLast=Concatenate_Vertices(VerticesAllbutLast,Vertices{k},outA,prod(nOutputs{k})); %add vertices of party k to overall vertices
    
    outA = outA*prod(nOutputs{k}); %number of overall outputs
    
end


%produces a warning if a Choi state is used where there is only one final party
for k=1:length(WChoi)
    
    if length(dloc{WChoi(k)}) == 1
        
        disp('Warning: a Choi state is used on party'), disp(k), disp('while there is only one final correpsonding party (no need for Choi state)') 
        
    end
    
end
%%

%Run the heuristic search



HeuristicSearchBroadcastScript

