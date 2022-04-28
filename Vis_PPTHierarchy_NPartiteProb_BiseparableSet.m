
function [vis,BellIneq,SDPOutput] = Vis_PPTHierarchy_NPartiteProb_BiseparableSet(Num_Outputs_Inpust_EachPartie,p_T,p_W,level,sys,sdpoptions)
% Suppose that the Bell scenario is (Oa,Ob,Oc,Sa,Sb,Sc)
% Num_Outputs_Input_EachPartie: [Oa,Ob,Oc;Sa,Sb,Sc]
% p_T,p_W: are in the following form  
% [P(a=1,b=1,c=1|x=1,y=1,z=1), P(a=2,b=1,c=1|x=1,y=1,z=1), ..., P(a=1,b=Oa-1,c=1|x=1,y=1,z=1),P(a=1,b=1,c=1|x=2,y=1,z=1), ... , P(a=Oa-1,b=1,c=1|x=Sa,y=1,z=1), P(a=1,b=1,c=1|x=Sa+1,y=1,z=1), P(a=1,b=1,c=1|x=1,y=2,z=1),
%  ... ,P(a=Oa-1,b=1,c=1|x=Sa,y=2,z=1), P(a=1,b=1,c=1|x=Sa+1,y=2,z=1),P(a=1,b=1,c=1|x=1,y=3,z=1), ... , P(a=1,b=1,c=1|x=Sa+1,y=Sb+1,z=Sz+1)]
% when the input index is Sa+1/ Sb+1/ Sc+1 then it stands for the
% identity operator of the corresponding party. The last term is then the constant term
% of the given Bell inequality
% level: is the level of the Moroder hierarchy wished to be implemented. 
% sys: is the index of the parties grouped as one group. For example, in a
% 4-partite scenario, sys = [1,3] means to solve the problem in the
% situation where parties are separated into [1,3] || [2,4]
% sdpoptions: is the sdpsetting used for yalmip
% requires: yalmip, vec, mat that comes with SeDuMi, Tx from Toby Cubitt's website,
% Generate_Moroder_Hierarchy_of_MomentMatrix.m and files that depend on these

%% 

% Extract the number of input and outputs from the probabilities given

N = size(Num_Outputs_Inpust_EachPartie,2);
TmpSys = N:-1:1;
sys = TmpSys(sys);

% Create the Moroder moment matrix
[MomentMatrix,Operator_gives_PhysicalMoment,TotalNumOperatorNeeded_EachPartie] = Generate_Moroder_Hierarchy_of_MomentMatrix(level,Num_Outputs_Inpust_EachPartie);



% This gives the number of starting operators used to construct the
% moment matrix; it is also the size of the moment matrix
% D0=length(Operator_gives_PhysicalMoment);
D0 = size(MomentMatrix,1);
dLocal=round(D0^(1/N));

% This gives the number of distinct variables required to run the SDP
% relaxation

F=[];

% Declare the decision variables
MomentVar = sdpvar(1,max(MomentMatrix(:)));

% Declare the decision variables
%MomentVar = sdpvar(1,max(MomentMatrix(:)));

% This ensures that the moment matrix has its (1,1) equals to 1, which
% enforces the normalization of the state giving rise to those moments
MomentVar(1) = 1; 
% Construct the moment matrix constraints
MomentConstraint = zeros(size(MomentMatrix));
for i = 1:max(MomentMatrix(:))
    MomentConstraint = MomentConstraint + MomentVar(i)*(MomentMatrix==i);
end
% Perform partial transposition of the matrix
MomentConstraint_Pt = Tx(MomentConstraint,sys,dLocal*ones(1,N));

 %MomentConstraint_Pt = Tx(MomentConstraint,sys,TotalNumOperatorNeeded_EachPartie);
 
% Now sets the PSD constraint of the moment matrix, and the moment matrices
% that result from
F = F + [MomentVar(Operator_gives_PhysicalMoment(:,end)) >=0];
F=F+[MomentConstraint >= 0];
F = F + [MomentConstraint_Pt>=0];
vis = sdpvar(1);
F = F + [vis >=0]; 
PTarget = vis*p_T + (1-vis)*p_W;
BellIneq_Constraint = [];
for i =1:size(Operator_gives_PhysicalMoment,1)-1;
    BellIneq_Constraint= BellIneq_Constraint+ [MomentVar(Operator_gives_PhysicalMoment(i,end)) == PTarget(i)];
end
F = F + [(BellIneq_Constraint):'BellIneq'];
SDPOutput=solvesdp(F,-vis,sdpoptions);
for i =1:size(Operator_gives_PhysicalMoment,1)-1;
    BellIneq(i)= dual(BellIneq_Constraint(i));
end
BellIneq(i+1)=0;
vis=double(vis);
yalmip('clear')

