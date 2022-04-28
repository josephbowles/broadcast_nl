% First test the Tsirelson point in the CHSH scenario with white noise
%% Create the Tsirelson point and white noise in full probability distribution form, i.e., P(a,b|x,y) 
rho=[0 0 0 0; 0 0.5 -0.5 0 ; 0 -0.5 0.5 0 ; 0 0 0 0];
Ma=cell(2,2);
Mb=cell(2,2);
pauliZ=[1 0 ; 0 -1];
pauliX=[0 1 ; 1 0 ];
pauliY=[0 -1i;1i 0 ];
Ma{1,1}=0.5*(eye(2)+pauliZ);
Ma{2,1}=eye(2)-Ma{1,1};
Ma{1,2}=0.5*(eye(2)+pauliX);
Ma{2,2}=eye(2)-Ma{1,2};
Mb{1,1}=0.5*(eye(2)-(pauliZ+pauliX)/2^0.5);
Mb{2,1}=eye(2)-Mb{1,1};
Mb{1,2}=0.5*(eye(2)+(pauliX-pauliZ)/2^0.5);
Mb{2,2}=eye(2)-Mb{1,2};

p_T=zeros(2,2,2,2);
for i=1:16
    [a,b,x,y]=ind2sub([2,2,2,2],i);
    p_T(a,b,x,y) =trace(rho*kron(Ma{a,x},Mb{b,y}));
end
old_p_T = p_T;
p_W = ones(2,2,2,2)/4;

%% Trasnform full P(a,b|x,y) into the CG form
% POutput = Transform_NPartite_ProbabilityDistributions(type,Scenario,PTarget,Normalization)
p_T = Transform_NPartite_ProbabilityDistributions('full',[2,2;2,2;],p_T);
p_W = Transform_NPartite_ProbabilityDistributions('full',[2,2;2,2;],p_W);
%% Solve the visibility to the biseparable set 
% [vis,BellIneq,SDPOutput] = Vis_PPTHierarchy_NPartiteProb_BiseparableSet(Num_Outputs_Inpust_EachPartie,p_T,p_W,level,sys,sdpoptions)
options=sdpsettings;
[Vis,BellIneq] = Vis_PPTHierarchy_NPartiteProb_BiseparableSet([2,2;2,2],p_T(:),p_W(:),1,[1],options);
assert(abs(Vis-1/sqrt(2))<=1e-6)
CHSHIneq = zeros(2,2,2,2);
CHSHIneq(:,:,1,1) = [1,-1;-1,1];
CHSHIneq(:,:,2,1) = [1,-1;-1,1]; CHSHIneq(:,:,1,2) = [1,-1;-1,1]; CHSHIneq(:,:,2,2) = -[1,-1;-1,1];
CHSHIneq_CGForm = Transform_NPartite_BellInequaliy('full',[2,2;2,2],CHSHIneq);
BellIneq = BellIneq/sqrt(2)*-4;
[BellIneq(:),CHSHIneq_CGForm(:)];

%% Lift the Tsirelson point to a 3 party scenario
p_T = zeros(2,2,2,2,2,2); 
p_T(:,:,1,:,:,1) = old_p_T; 
p_T(:,:,1,:,:,2) = old_p_T;
p_W = ones(2,2,2,2,2,2)/8;
p_W = Transform_NPartite_ProbabilityDistributions('full',[2,2;2,2;2,2],p_W);
p_T = Transform_NPartite_ProbabilityDistributions('full',[2,2;2,2;2,2],p_T);
[Vis1,BellIneq1] = Vis_PPTHierarchy_NPartiteProb_BiseparableSet([2,2,2;2,2,2],p_T,p_W,1,[3],options);
assert(abs(Vis1-1)<=1e-5);
[Vis2,BellIneq2] = Vis_PPTHierarchy_NPartiteProb_BiseparableSet([2,2,2;2,2,2],p_T,p_W,1,[2],options);

Vis1=Vis1
Vis2=Vis2


