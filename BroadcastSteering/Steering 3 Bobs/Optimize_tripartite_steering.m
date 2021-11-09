%Finds the smallest visibility 'v' for which v*rho+(1-v)*rhoN is steerable
%in the broadcast scenario, using three Bobs 



%scenario

ny=2; nz=2; nw=2;
scenario = [ny nz nw];


%vertices

Vertices=facets4fla(:,:,:,1:300); 



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


%creates the maximally mixed assemblage:

sigId=zeros(dA,dA,ny+1,nz+1,nw+1);
sigId(:,:,ny+1,nz+1,nw+1)=eye(dA)/dA;




%first finds a point outside:

vis=2;
while vis > 1
    

%random Choi and local POVMs

T = RandomChoi(d,d^3); 



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
time = toc


vis=f{1};

end



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



%start the heuristic optimization


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
    
    T=ForceSDP(f{2}); B=make_povms_valid(f{3}); C=make_povms_valid(f{4}); D=make_povms_valid(f{5});
    
   
    

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


%difference between two succesive values
    gap=visc-vis;
    
    visc=vis;
    
end
