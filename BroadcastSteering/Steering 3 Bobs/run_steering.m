
%scenario

ny=2; nz=2; nw=2;


rho = IsotropicState(2,1);

rhoN = eye(4)/4;

dA=2; 
dB0=2;
dB=2;
dC=2;
dD=2;


T = RandomSuperoperator([2 8]); 

ketV=zeros(1,8*2);
ketV(1)=1;
ketV(8*2)=1;
T=ketV'*ketV;

scenario = [ny nz nw];


B=rand_povms(dB,2,ny,1);

C=rand_povms(dC,2,nz,1);

D=rand_povms(dD,2,nw,1);


sigrho = Broadcast_TripartiteAssemblage_FC(rho,dA,T,B,C,D);

sigN = Broadcast_TripartiteAssemblage_FC(rhoN,dA,T,B,C,D);


Vertices=facets4fla(:,:,:,1:300); 

tic
fP = TripartiteLHS_FC(sigrho,sigN,scenario,Vertices); 

time_primal=toc

tic
fD = TripartiteLHS_FC_DUAL(sigrho,sigN,scenario,Vertices) 

time_dual=toc