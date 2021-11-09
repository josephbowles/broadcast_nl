%Set the parameters for SteeringHeuristicSearchBroadcast with 2 Bobs

%This script is meant to set all the variables needed to run a heuristic
%search in the broadcast scenario

%scenario


ny=2; nz=2; nw=2;
scenario = [ny nz nw];
%state

d=2;
rho = IsotropicState(d,1);
rhoN = eye(4)/4;
%local initial and final dimensions
dA=2;
dB0=2;
dB=2;
dC=2;
dD=2;


IaIbIc=[ny nz nw];
if IaIbIc == [2 2 2]
    load VNS3BOBS;
    Vertices = VNS3BOBS(:,:,:,1:end);   %Variable which contains the vertices of the scenario
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    'Please set manually the vertices of the NS polytope'
    return
end


run3Bobs