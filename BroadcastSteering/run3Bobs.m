%Finds the smallest visibility 'v' for which v*rho+(1-v)*sig is steerable
%with 3 Bobs

%all_info contains a cell with all relevant information
%all_info{1,1} = all visibility values
%all_info{1,2} = sigma
%all_info{1,3} = scenario

%all_info{2,1} = eta
%all_info{2,2} = time4_evaluating_eta
%all_info{2,3} = B1
%all_info{2,4} = B2
%all_info{2,5} = T
%all_info{2,6} = Witness

rng('shuffle')

C=clock;
clock_string=strcat(int2str(C(1)),'_',int2str(C(2)),'_',int2str(C(3)),'_',int2str(C(4)),'_',int2str(C(5)),'_',int2str(round(C(6))))
info='3Bobs2InputsQubitIstropic_';
file_name=strcat(info,clock_string);

%Uncomment if you desire to save the useful .mat variables at the computer
% if ~exist('3Bobs2Inputs', 'dir')
%     mkdir 3Bobs2Inputs
% end


%creates the maximally mixed assemblage:
sigId=zeros(dA,dA,ny+1,nz+1,nw+1);
sigId(:,:,ny+1,nz+1,nw+1)=eye(dA)/dA;


all_info{2,1}=rho;
all_info{3,1}=rhoN;
all_info{4,1}=scenario;    

%random Choi and local POVMs
ketV=zeros(1,8*2);
ketV(1)=1;
ketV(8*2)=1;
T=ketV'*ketV;
U=RandomUnitary(2^3);
T=kron(eye(2),U)*T*kron(eye(2),U');

B=rand_povms(dB,2,ny,1);
C=rand_povms(dC,2,nz,1);
D=rand_povms(dD,2,nw,1);

%construct the two assemblages coming from rho and rhoN
sigrho = Broadcast_TripartiteAssemblage_FC(rho,dA,T,B,C,D);
sigN = Broadcast_TripartiteAssemblage_FC(rhoN,dA,T,B,C,D);

%find the critical visibility (and corresponding inequality)
tic;
qId=10^-6; 
f = TripartiteLHS_FC(sigrho,(1-qId)*sigN+qId*sigId,scenario,Vertices); 
time_SDP = toc;
vis=f{1};


%extract inequality and compute its bound

In=f{2}; In=make_hermitian(In);

q=vis; pq = q*sigrho + (1-q)*sigN;

S=0;
for y=1:ny+1
    for z=1:nz+1
        for w=1:nw+1
                S = S + In(:,:,y,z,w)*pq(:,:,y,z,w);
        end
    end
end

sc=trace(S);



all_info{1,2}=vis;
all_info{2,2}=time_SDP;
all_info{3,2}=B;
all_info{4,2}=C;
all_info{5,2}=D;
all_info{6,2}=T;
all_info{7,2}=In;
%Uncomment if you desire to save the useful .mat variables at the computer
% cd 3Bobs2Inputs
% save(file_name,'all_info')
% cd ..


%current visibility
visc=vis;


%store the successive visibilities in vector 'rep'
rep=0;
rep(1)=vis;
%stop when two successive visibilities are closer than 'prec'
prec=10^-5;

gap=1; k=2;
while gap > prec
    
    
    %current unsteerable state
    rhoc=visc*rho+(1-visc)*rhoN;
    
    %current settings
    init={T,C,D};
    
    %maximizes the inequality with a seesaw. tries 'limcount' times before
    %giving up 
    limcount=7;
    
    viol=sc-1; count=0; 
    while viol <= sc
        
        %for count==3 uses the previous settings as initial point of the
        %seesaw
        if count == 3
            f=SeesawTripartiteSteering_inequFC(In,rhoc,dA,scenario,dB,dC,dD,1,init);
        else
            f=SeesawTripartiteSteering_inequFC(In,rhoc,dA,scenario,dB,dC,dD,1);
        end
      
        yalmip clear
        
        viol=f{1};
        
        count=count+1;
        if count > limcount
            break
        end 
    end
    %extracts answer of the seesaw (and force SDP + normalization for all
    %variables)
    
    T=get_valid_channel(f{2},2,8); B=make_povms_valid(f{3}); C=make_povms_valid(f{4}); D=make_povms_valid(f{5});
    
%construct the two assemblages coming from rho and rhoN

sigrho = Broadcast_TripartiteAssemblage_FC(rho,dA,T,B,C,D);
sigN = Broadcast_TripartiteAssemblage_FC(rhoN,dA,T,B,C,D);

%find the critical visibility (and corresponding inequality)

qId=10^-6; 
f = TripartiteLHS_FC(sigrho,(1-qId)*sigN+qId*sigId,scenario,Vertices);

%new visibility
vis=f{1}
    rep(k)=vis;
    k=k+1;
    %new inequality
    In=f{2}; In=make_hermitian(In);
    
    %compute local bound of new inequality
S=0;

for y=1:ny+1
    for z=1:nz+1
        for w=1:nw+1
            
                
                S = S + In(:,:,y,z,w)*pq(:,:,y,z,w);

            
        end
    end
end


sc=trace(S);

all_info{1,1}=rep;
all_info{1,k}=vis;
all_info{2,k}=time_SDP;
all_info{3,k}=B;
all_info{4,k}=C;
all_info{5,k}=D;
all_info{6,k}=T;
all_info{7,k}=In;

%Uncomment if you desire to save the useful .mat variables at the computer
% cd 3Bobs2Inputs
% save(file_name,'all_info')
% cd ..

%difference between two succesive values
    gap=visc-vis;
    
    visc=vis;
    
end
