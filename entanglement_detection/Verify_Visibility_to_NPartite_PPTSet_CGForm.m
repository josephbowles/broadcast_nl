function [Vis,Diag,BellIneq]=  Verify_Visibility_to_NPartite_PPTSet_CGForm (Scenario,P_CGV,Level,Grouping,options,P_Noise_CGV,MomentData)
% This code solves the visibility of a given correlation to the PPT set, characterized by the Moroder et. al. hierarchy with PPT constraints. 
% The argument 'Scenario' is a 2-dimensional matrix and is in the form
% [Oa|x=1, Oa|x=2, Oa|x=3, ?, Oa|x=N_x,
%  Ob|y=1, Ob|y=2, Ob|y=3, ?, Ob|y=N_y,
%  Oc|z=1, Oc|z=2, Oc|z=3  ...]
% 'P_CGV' is the CG vector of a correlation
% 'Level' is to specify the maximum order of the chosen local operators
% of each party to construct the moment matrix. It can be an arbitrary
% positive number.
% 
% 'Grouping' is a cell matrix in the follwoing format '{{{ },{ }, { } }, { { },{ },...,{ }}, { { },{ },..,{ } }, ..., {{ }, { }, ...,{ } } }'
% The length of the first layer, length(Grouping), specifies the number of
% DIFFERENT moment matrices is used in the following computation. 
% The length of the second layer, length(Grouping{k}), specifies the number
% of DIFFERENT PPT constraints imposed on the kth moment matrix. 
% For example, in a three party case, the following Grouping characterized
% the biseparable set of correlations { { {1} },{ {2} },{ {3} } }' 
%  Similaraly, with 'Grouping' being as {{{1},{2},{3}}}, an outer
%  approximation to the local polytope is constructed. 
% 
% 'options' is the sdpsettings for yalmip
% 'P_Noise_CGV' is the CG vector of the noisy distribution. By default it
% is the white-noise but a different distribution can be assigned. 
% 
% 'MomentData' is the data used to construct the moment matrix. If it is
% given, then one does not need to generater/load the corresponding data
% while sovling a huge number of points again and again. 
% It can be generated first by 
% [MomentData.FinalMomentMatrix,MomentData.Operator_gives_PhysicalMoment_abcxyz,~,~,MomentData.TotalNumOperatorNeeded_EachPartie] = Generate_Moroder_Hierarchy_of_MomentMatrix(Scenario,Level);

Normalization = P_CGV(end);
FinalMomentMatrix = MomentData.FinalMomentMatrix;
Operator_gives_PhysicalMoment_abcxyz = MomentData.Operator_gives_PhysicalMoment_abcxyz;
TotalNumOperatorNeeded_EachPartie = MomentData.TotalNumOperatorNeeded_EachPartie;
if find(Operator_gives_PhysicalMoment_abcxyz(:,end)==0)
    error('Not every physical terms are used in current approximation to quantum set. Try using higer order or switch to Moroder hierarchy')
end

NumMoment = max(vec(FinalMomentMatrix));
NumDifferntGrouping = length(Grouping);
N = size(Scenario,1);
TmpSys = N:-1:1;
F = [];
DistinctOpt = cell(1,NumDifferntGrouping);
for i =1:NumDifferntGrouping
    DistinctOpt{i} = sdpvar(max(vec(FinalMomentMatrix)),1);
end
NormalizationConstraint = 0;
for i =1:NumDifferntGrouping
    NormalizationConstraint = NormalizationConstraint + DistinctOpt{i}(1);
end
F = F+[NormalizationConstraint == Normalization];
PDecomp = cell(1,NumDifferntGrouping);
MomentConstraint = cell(1,NumDifferntGrouping);
MomentConstraint_PartialTranspose = cell(1,NumDifferntGrouping);
for i =1:NumDifferntGrouping
    MomentConstraint{i} =  zeros(size(FinalMomentMatrix));
    for k=1:NumMoment
        MomentConstraint{i} = MomentConstraint{i} + DistinctOpt{i}(k)*(FinalMomentMatrix==k);
    end
    F = F +[MomentConstraint{i}>=0];
        
    PDecomp{i} = DistinctOpt{i}(Operator_gives_PhysicalMoment_abcxyz(:,end));
    MomentConstraint_PartialTranspose{i} = cell(1,length(Grouping{i}));
    for k =1:length(Grouping{i})
        MomentConstraint_PartialTranspose{i}{k} = Tx(MomentConstraint{i},TmpSys(cell2mat(Grouping{i}{k})),TotalNumOperatorNeeded_EachPartie);
       F = F +[MomentConstraint_PartialTranspose{i}{k} >=0];
    end
    
end

Vis = sdpvar(1);
F = F +[Vis>=0];
P_V = Vis*P_CGV(:) + (1-Vis)*P_Noise_CGV;
PVar = 0 ;
for i =1:NumDifferntGrouping
    PVar = PVar + PDecomp{i};
end
BellIneq_Constraint = [ ];
for i =1:size(Operator_gives_PhysicalMoment_abcxyz,1)-1
    BellIneq_Constraint = BellIneq_Constraint + [PVar(i) == P_V(i)];
end
F = F + [(BellIneq_Constraint):'BellCoeff'];
% end

Diag=solvesdp(F,-Vis,options);
BellIneq = zeros(size(Operator_gives_PhysicalMoment_abcxyz,1)-1,1);
for i =1:size(Operator_gives_PhysicalMoment_abcxyz,1)-1
    BellIneq(i) = dual(BellIneq_Constraint(i));
end
BellIneq = [BellIneq;0];
Vis = double(Vis);
yalmip('clear')
end