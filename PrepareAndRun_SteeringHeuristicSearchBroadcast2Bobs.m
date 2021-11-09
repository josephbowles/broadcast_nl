%Set the parameters for SteeringHeuristicSearchBroadcast with 2 Bobs

%This script is meant to set all the variables needed to run a heuristic
%search in the broadcast scenario

ob1=2; ny1=3;  ob2=2; ny2=2;   %number of inputs and outputs for each Bob
scenario=[ob1 ob2 ny1 ny2];

%initial state
d=2;
rho=IsotropicState(d,1);
sig=eye(d^2)/d^2;
dA=d; dB0=size(rho,1)/dA; dB1=dB0; dB2=dB0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IaIb=[ny1 ny2];
if IaIb == [2 2]
    load VNS22;
    DB=VNS22;  %Variable which contains the vertices of the scenario
elseif IaIb == [3 2]
    load VNS32;
    DB=VNS32;  %Variable which contains the vertices of the scenario
elseif IaIb == [3 3]
    load VNS33;
    DB=VNS33;  %Variable which contains the vertices of the scenario
elseif IaIb == [4 4]
    load VNS44;
    DB=VNS44;  %Variable which contains the vertices of the scenario
else
    'Please set manually the vertices of the NS polytope'
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run2Bobs