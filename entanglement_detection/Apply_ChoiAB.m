function f=Apply_ChoiAB(rho,TA,TB,dA,UA12,UB12)

%Applies the map corresponding to Choi state 'TA' for Alice, 'TB' for Bob, on state 'rho' (of local dimension 'dA', rho assumed symmetric if not specified) 
%(optional: Unitary UA12 on the output space of Alice, Unitary UB12 on the output space of Bob)

%'rho' is a dA*dB x dA*dB matrix, representing a bipartite quantum state
%'TA' is a dA*dAf x dA*dAf matrix, representing the Choi state of Alice 
%'TB' is a dB*dBf x dB*dBf matrix, representing the Choi state of Bob 

%optional variable:
%'dA' is Alice's local dimension


if nargin < 4
    
    dA0=sqrt(size(rho,1));
    dB0=size(rho,1)/dA0;  dAf=size(TA,1)/dA0; dBf=size(TB,1)/dB0;
    UA12=eye(dAf);
    UB12=eye(dBf);
else
    
    dA0=dA;
    
end

 dB0=size(rho,1)/dA0;  dAf=size(TA,1)/dA0; dBf=size(TB,1)/dB0;

  
if nargin < 6
    
    UA12=eye(dAf);
    UB12=eye(dBf);
    
end



%applies isometry

%Bob
TOT=kron(PartialTranspose(rho,2,[dA0 dB0]),eye(dBf))*kron(eye(dA0),TB); 

rhotb=kron(eye(dA0),UB12)*PartialTrace(TOT,2,[dA0 dB0 dBf])*kron(eye(dA0),UB12');

%Alice
dBt=dBf; 
rhos=Swap(rhotb,[1,2],[dA0 dBt]);

TOT=kron(PartialTranspose(rhos,2,[dBt dA0]),eye(dAf))*kron(eye(dBt),TA); 

rhots=PartialTrace(TOT,2,[dBt dA0 dAf]);

rhot=kron(UA12,eye(dBf))*Swap(rhots,[1,2],[dBt dAf])*kron(UA12',eye(dBf)); 


f=rhot;

end
