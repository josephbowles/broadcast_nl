function f = MultipartiteBehaviour(rho,Meas)

%Computes the multipartite behaviour p(abc...|xyz...) coming from state 'rho' and local measurements 'Meas'

%'rho' is a (dA*dB*dC..)x(dA*dB*dC..) psd matrix, representing a multipartite quantum state
%'Meas' is a 1xN cell, each entry containing local POVMs in table form: Meas = { A(dA,dA,oa,nx), B(dB,dB,oa,ny), ... }
%That is, A(dA,dA,oa,nx) is a dA x dA x oa x nx matrix, such that
%A(:,:,a,x)>=0 and sum_n(A(:,:,:,x),3)=eye(dA) for all x=1,..,nx



%local dims and #inputs&outputs

N=length(Meas);

if N > 1

dvec=zeros(1,N); oa=zeros(1,N); nx=zeros(1,N);
ind=1;

for k=1:N
    
    dvec(k) = size(Meas{k},1);
    oa(k)=size(Meas{k},3);
    nx(k)=size(Meas{k},4);
    
end


p=zeros(1,prod(oa)*prod(nx));

FM=Meas{1};
LM=Meas{N};
        
vx=ones(1,N); vx(end)=0;
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
        
        p(ind)=real(trace(superkron(S0,LM(:,:,va(N),vx(N)))*rho));
        
        ind=ind+1;
        
    end
    
end




elseif N==1
    
    M=Meas{1};
    
    k=1;  p=zeros(1,size(M,4)*size(M,3));
    for x=1:size(M,4)
        for a=1:size(M,3)
            
            p(k)=trace(M(:,:,a,x)*rho);
            k=k+1;
            
        end
    end
    
    
end




f=p;


end

