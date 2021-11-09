function sig=Broadcast_TripartiteAssemblage_FC(rho,dA,T,B,C,D)

%In the broadcast scenario with 2 outputs, computes the overall "full-correlator assemblage" sig(b c d |y z w) coming from state 'rho' (of local dimension dA for Alice), Choi state 'T', and local measurements, 'B' and 'C' and 'D', for Bob, Charlie and David

%Inputs: 

%'rho' is a (dA*dB0)x(dA*dB0) SDP matrix, corresponding to a bipartite state
%'dA' is a scalar, corresponding to Alice's local dimension
%'T' is a (dB0*dB*dC*dD)x(dB0*dB*dC*dD) SDP matrix, corresponding to a Choi state
%'B' is a dB x dB x 2 x ny matrix, corresponding to Bob's local POVMs
%'C' is a dB x dB x 2 x nz matrix, corresponding to Charlie's local POVMs
%'D' is a dB x dB x 2 x nw matrix, corresponding to David's local POVMs

%Output:

%'sig' is a dA x dA x (ny+1) x (nz+1) x (nw+1) matrix, corresponding to a tripartite
%assemblage, given in the full-correlation convention, where the last input
%of each party corresponds to the identity measurement (that is,
%sig(:,:,y,nz+1,nw+1) gives Bob's marginal, etc. )



%extracts scenario

ny=size(B,4);

nz=size(C,4);

nw=size(D,4);


dB0=length(rho)/dA;
dB=size(B,1);
dC=size(C,2);
dD=size(D,2);



%Define observables

OB = zeros(dB,dB,ny+1);
for y=1:ny
    
    OB(:,:,y) = B(:,:,1,y) - B(:,:,2,y);
    
end
OB(:,:,ny+1) =eye(dB);



OC = zeros(dC,dC,nz+1);
for z=1:nz
    
    OC(:,:,z) = C(:,:,1,z) - C(:,:,2,z);
    
end
OC(:,:,nz+1) =eye(dC);




OD = zeros(dD,dD,nw+1);
for w=1:nw
    
    OD(:,:,w) = D(:,:,1,w) - D(:,:,2,w);
    
end
OD(:,:,nw+1) =eye(dD);


%applies isometry

TOT=kron(PartialTranspose(rho,2,[dA dB0]),eye(dB*dC*dD))*kron(eye(dA),T);

rhot=PartialTrace(TOT,2,[dA dB0 dB dC dD]);



%computes the assemblage

sig=zeros(dA,dA,ny,nz,nw);



for y=1:ny+1
    for z=1:nz+1
        for w=1:nw+1

            
            sig(:,:,y,z,w)=PartialTrace(Tensor(eye(dA),OB(:,:,y),OC(:,:,z),OD(:,:,w))*rhot,[2 3 4],[dA dB dC dD]);
            
        end
        
    end
end







end