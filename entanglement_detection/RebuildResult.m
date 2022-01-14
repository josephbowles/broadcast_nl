load PPTFinaldata.mat 
rho=Wern(1);
sig = eye(4)/4;
Scenario = ones(4,3)*2;
nx1 = 3; oa1 = 2; 
nx2 = 3; oa2 = 2; 
ny1 = 3; ob1 = 2;
ny2 = 3; ob2 = 2; 
dA0 = 2; 

prho = Rafael_pointAABB(rho,dA0,TA,A1,A2,TB,B1,B2);
prho = reshape(prho,[ob2,ob1,oa2,oa1,ny2,ny1,nx2,nx1]);
prho = permute(prho,[4,3,2,1,8,7,6,5]);
psigN= Rafael_pointAABB(sig,dA0,TA,A1,A2,TB,B1,B2);
psigN = reshape(psigN,[ob2,ob1,oa2,oa1,ny2,ny1,nx2,nx1]);
psigN= permute(psigN,[4,3,2,1,8,7,6,5]);
Grouping = {{{[1,2]}}};
%transform to Collins-Gisin
prhoT = Transform_NPartite_CGVector_between_MultiDimensionalMatrix('CGM',Scenario,Transform_NPartite_ProbabilityDistributions('full',Scenario,prho));
psigN = Transform_NPartite_CGVector_between_MultiDimensionalMatrix('CGM',Scenario,Transform_NPartite_ProbabilityDistributions('full',Scenario,psigN));
options = sdpsettings; 
options.solver = 'mosek';
% Primial feasibility for Conic and Linear problems, default 1e-8
options.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1e-6;
options.mosek.MSK_DPAR_INTPNT_TOL_PFEAS = 1e-6;
% Dual feasibility for Conic and Linear problems, default 1e-8
options.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS = 1e-6;
options.mosek.MSK_DPAR_INTPNT_TOL_DFEAS = 1e-6;
% Dualiy and complementarity gap for Conic and Linear problems, default 1e-8,1e-8,1e-16,1e-8,
options.mosek.MSK_DPAR_INTPNT_CO_TOL_MU_RED = 1e-6;
options.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-6;
%     options.mosek.MSK_DPAR_INTPNT_TOL_MU_RED = eps;
options.mosek.MSK_DPAR_INTPNT_TOL_REL_GAP=1e-8;
% Infeasibility for Conic and Linear problems, default 1e-10
% Bigeps when above precisions are not reached
options.mosek.MSK_DPAR_INTPNT_CO_TOL_NEAR_REL=1e3;
% Maximal number of interactions
options.mosek.MSK_IPAR_INTPNT_MAX_ITERATIONS = 3000;

% Verify visibility to PPT set with level 1 of the Moroder hierarchy 

MomentMatrixData = load('MomentmatrixData_Scenario_222_222_222_222_LocalLevel_1.mat');
[vis_level1,SDPOutput_level1,Dual_level1] =  Verify_Visibility_to_NPartite_PPTSet_CGForm(Scenario,vec(prhoT),1,Grouping,options,vec(psigN),MomentMatrixData);

% Verify visibility to PPT set with level 1.2 of the Moroder hierarchy to
% obtain a better visibiliy value 
MomentMatrixData = load('MomentmatrixData_Scenario_222_222_222_222_LocalLevel_1.2.mat');
[vis,SDPOutput,Dual] =  Verify_Visibility_to_NPartite_PPTSet_CGForm(Scenario,vec(prhoT),1.2,Grouping,options,vec(psigN),MomentMatrixData);
