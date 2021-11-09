function f=SeesawTripartiteSteering_inequFC(ineq,rho,dA,scenario,dB,dC,dD,choiin,Initial) 

%In the broadcast scenario with three Bobs (with Bobs' dimensions 'dB', 'dC', 'dD') optimizes steering inequality 'ineq' for state 'rho' (of local dim 'dA'), in scenario 'scenario' (starting from optional points 'Initial')

%'ineq' is a dA x dA x (ny+1) x (nz+1) x (nw+1) matrix, interpreted as the steering inequality matrices I_{y,z,w} in full-correlator notation
%'rho' is a bipartite state
%'scenario' is of the form [ny nz nw] 
%'choiin' is an integer, wich tells whether the initial Choi state is taken as a projector (if choiin==1), or a generic Choi state (other values)

%optional variable:
%'Initial' is a 1x2 cell, containing initials values of T,C,D. That is, Initial = {T,C,D}


%Extracts scenario 

ny=scenario(1); nz=scenario(2); nw=scenario(3); 

dB0=size(rho,1)/dA;


It=ineq;

%initial point


if nargin < 9
    
    if choiin == 1
        
        Tc=RandomChoi(dB0,dB*dC*dD);
        
    else
        
        Tc=RandomChoi(dB0,dB*dC*dD);
        
    end
    
    
    
    Cc=rand_povms(dC,2,nz,1);
    Dc=rand_povms(dD,2,nw,1);
else
    
    Tc=Initial{1};  Cc=Initial{2}; Dc=Initial{3};
   
    
end





prec=10^-5;

former=1i;
c=1;
while (c > prec)


rr=optineqB(Tc,Cc,Dc);

Bc=rr{2};


rr=optineqC(Tc,Bc,Dc);

Cc=rr{2};


rr=optineqD(Tc,Bc,Cc);

Dc=rr{2};


rr=optineqT(Bc,Cc,Dc);

Tc=rr{2};


c = abs(rr{1} - former);

former=rr{1};

end



f={rr{1},Tc,Bc,Cc,Dc};



%%%sub-functions: 




function yy=optineqB(T,C,D)

    


%sdp variables

B=sdpvar(dB,dB,2,ny,'hermitian','complex');

R=[];

for x=1:ny
    for a=1:2
        
        R = R + (B(:,:,a,x) >= 0);
        
    end
    
    R = R + (sum(B(:,:,:,x),3)==eye(dB) ) ;
    
end



%Define observables


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

   
    %compute current score with new inequ
S=0;

for y=1:ny+1
    for z=1:nz+1
        for w=1:nw+1
            
            if y==ny+1
                OBB = eye(dB);
            else
                 OBB = B(:,:,1,y) - B(:,:,2,y);
            end
                
                S = S + trace(It(:,:,y,z,w)*PartialTrace(Tensor(eye(dA),OBB,OC(:,:,z),OD(:,:,w))*rhot,[2 3 4],[dA dB dC dD]));

            
        end
    end
end


O = real(S); 




 ops=sdpsettings('verbose',0,'warning',0,'solver','mosek');
 optimize(R,-O,ops);
 
 
yy={double(O),double(B)};

end

function yy=optineqC(T,B,D)

    


%sdp variables

C=sdpvar(dC,dC,2,nz,'hermitian','complex');

R=[];

for x=1:nz
    for a=1:2
        
        R = R + (C(:,:,a,x) >= 0);
        
    end
    
    R = R + (sum(C(:,:,:,x),3)==eye(dC) ) ;
    
end

%Define observables


OB = zeros(dC,dC,ny+1);
for y=1:ny
    
    OB(:,:,y) = B(:,:,1,y) - B(:,:,2,y);
    
end
OB(:,:,ny+1) =eye(dB);


OD = zeros(dD,dD,nw+1);
for w=1:nw
    
    OD(:,:,w) = D(:,:,1,w) - D(:,:,2,w);
    
end
OD(:,:,nw+1) =eye(dD);


%applies isometry

TOT=kron(PartialTranspose(rho,2,[dA dB0]),eye(dB*dC*dD))*kron(eye(dA),T);

rhot=PartialTrace(TOT,2,[dA dB0 dB dC dD]);

   
    %compute current score with new inequ
S=0;

for y=1:ny+1
    for z=1:nz+1
        for w=1:nw+1
            
            if z==nz+1
                OCC = eye(dC);
            else
                 OCC = C(:,:,1,z) - C(:,:,2,z);
            end
                
                S = S + trace(It(:,:,y,z,w)*PartialTrace(Tensor(eye(dA),OB(:,:,y),OCC,OD(:,:,w))*rhot,[2 3 4],[dA dB dC dD]));

            
        end
    end
end


O = real(S); 




 ops=sdpsettings('verbose',0,'warning',0,'solver','mosek');
 optimize(R,-O,ops);
 
 
yy={double(O),double(C)};

end

function yy=optineqD(T,B,C)

    


%sdp variables

D=sdpvar(dD,dD,2,nw,'hermitian','complex');

R=[];

for x=1:nw
    for a=1:2
        
        R = R + (D(:,:,a,x) >= 0);
        
    end
    
    R = R + (sum(D(:,:,:,x),3)==eye(dD) ) ;
    
end


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



%applies isometry

TOT=kron(PartialTranspose(rho,2,[dA dB0]),eye(dB*dC*dD))*kron(eye(dA),T);

rhot=PartialTrace(TOT,2,[dA dB0 dB dC dD]);



   
    %compute current score with new inequ
S=0;

for y=1:ny+1
    for z=1:nz+1
        for w=1:nw+1
            
            if w==nw+1
                ODD = eye(dD);
            else
                 ODD = D(:,:,1,w) - D(:,:,2,w);
            end
                
                S = S + trace(It(:,:,y,z,w)*PartialTrace(Tensor(eye(dA),OB(:,:,y),OC(:,:,z),ODD)*rhot,[2 3 4],[dA dB dC dD]));

            
        end
    end
end


O = real(S); 



 ops=sdpsettings('verbose',0,'warning',0,'solver','mosek');
 optimize(R,-O,ops);
 
 
yy={double(O),double(D)};

end


function yy=optineqT(B,C,D) 


%isometry

T=sdpvar(dB0*dB*dC*dD,dB0*dB*dC*dD,'hermitian','complex');


R = (T >= 0) + (PartialTrace(T,2,[dB0 dB*dC*dD]) == eye(dB0));



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



   
    %compute current score with new inequ
S=0;

for y=1:ny+1
    for z=1:nz+1
        for w=1:nw+1
            

                S = S + trace(It(:,:,y,z,w)*PartialTrace(Tensor(eye(dA),OB(:,:,y),OC(:,:,z),OD(:,:,w))*rhot,[2 3 4],[dA dB dC dD]));

            
        end
    end
end


O = real(S); 



 ops=sdpsettings('verbose',0,'warning',0,'solver','mosek');
 optimize(R,-O,ops);
 
 
yy={double(O),double(T)};

end


end

