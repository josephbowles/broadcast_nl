function f = MaxIneqBroadcastMeas(Ineq,rho,dims,Party,Choi,Meas)

%Optimizes broadcast inequality 'Ineq' (for state 'rho' of local dimensions 'dims') over measurements of party 'Party' (rest is set from 'Choi' and 'Meas')

%Inputs:

%'Ineq' is a 1 x oa*ob*oc*...*nx*ny*nz*... vector, interpreted as the coefficients c(abc...|xyz...) = [c(000..|000..) c(00..1|000) ... ]) of the inequality

%'rho' is a (dA*dB*dC..)x(dA*dB*dC..) PSD and trace-1 matrix, representing a multipartite quantum state

%'dims' is a 1Xn vector, containing the local dimensions: dims = [dA dB dC ...]

%'Party' is an integer, representing the party's labeling one wants to optimize over (that is, will maximize the value of 'Ineq' over measurements Meas{Party} )

%'Choi' is a 1xn cell, each entry containing a Choi state (that is, a PSD matrix T such that PartialTrace(T,2) = eye(d1) )

%'Meas' is a 1xN cell, each entry containing local POVMs in table form: Meas = { A(dA,dA,oa,nx), B(dB,dB,oa,ny), ... }
%That is, A(dA,dA,oa,nx) is a dA x dA x oa x nx matrix, such that
%A(:,:,a,x)>=0 and sum_n(A(:,:,:,x),3)=eye(dA) for all x=1,..,nx


%Outputs:

%'f' is a 1x2 cell, f{1} containing the maximal value of the inputed
%inequality, f{2} the corresponding measurements






%final state
dimsc=dims;
for k=1:length(dims)
    rho=ApplyChoiMap(rho,Choi{k},dimsc,k);
    dimsc(k)=dims(k)*size(rho,1)/prod(dimsc); 
end



%Create the (partial) Bell operator


N=length(Meas);

AllbutParty = 1:N;
AllbutParty(Party) = [];



Ac=Meas{Party};

dA=size(Ac,1); oac=size(Ac,3); nxc=size(Ac,4);


dvec=zeros(1,N); oa=zeros(1,N); nx=zeros(1,N);


for k=1:N
    
    dvec(k) = size(Meas{k},1);
    oa(k)=size(Meas{k},3);
    nx(k)=size(Meas{k},4);
    
end



It=convP2Table(Ineq,[nx oa]);


Meas{Party} = eye(dA);
nx(Party)=1;
oa(Party)=1;



FM=Meas{1};
LM=Meas{N};



%construct the partial measurements operator

vx=ones(1,N); vx(end)=0;

F=zeros(dA,dA,prod(nx)*prod(oa));

for LX=1:prod(nx)
    
    vx = update_vector(vx,nx);
    
    
    va=ones(1,N); va(end)=0;
    
    S0=FM(:,:,1,vx(1));
    for k=2:N-1
        M=Meas{k};
        S0=superkron(S0,M(:,:,1,vx(k)));
    end
    
    for LA=1:prod(oa)
        
        
        va=update_vector(va,oa);
        
        
        if va(N) == 1
            
            S0=FM(:,:,va(1),vx(1));
            for k=2:N-1
                M=Meas{k};
                S0=superkron(S0,M(:,:,va(k),vx(k)));
            end
            
            
        end
        
        
        posF = vec2ind([nx oa],[vx va]);
        
        F(:,:,posF) = PartialTrace(superkron(S0,LM(:,:,va(N),vx(N)))*rho,AllbutParty,dvec);
        
    end
end


%constructing the partial Bell operator 

BO=zeros(dA,dA,oac,nxc);

for xc=1:nxc
    for ac=1:oac
        
        vx=ones(1,N); vx(end)=0;
        
        for LX=1:prod(nx)
            
            vx = update_vector(vx,nx);
            
            
            va=ones(1,N); va(end)=0;          
            
            for LA=1:prod(oa)
                
                
                va=update_vector(va,oa);
                
                
                vx_all = vx;
                vx_all(Party) = xc;
                
                va_all = va;
                va_all(Party) = ac;
                
                pos =  vec2ind(size(It),[vx_all va_all]);
                
                posF =  vec2ind([nx oa],[vx va]);
               
                
                BO(:,:,ac,xc) = BO(:,:,ac,xc) + It(pos)*F(:,:,posF);
                
            end
            
        end
        
    end
    
end




%SDP


%sdp variables

A=sdpvar(dA,dA,oac,nxc,'hermitian','complex');

R=[];

S=0;
for x=1:nxc
    for a=1:oac
        
        
        R = R + (A(:,:,a,x) >= 0);
        
        S = S + real(trace(A(:,:,a,x)*BO(:,:,a,x)));
        
    end
    
    R = R + (sum(A(:,:,:,x),3)==eye(dA) ) ;
    
end


%Objective function

O = S;


ops=sdpsettings('verbose',0,'warning',0,'solver','mosek');
optimize(R,-O,ops);


f={double(O),double(A)};

end