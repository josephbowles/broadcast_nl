function f=make_povms_valid(A)

%Makes the set of POVMs 'A' valid 

%'A' is a d x d x oa x nx matrix, representing a set of 'nx' oa-outcome POVMs of dimension d


L=size(A); 

d=L(1); oa=L(3);

if length(L) < 4
    
    nx=1;
    
else
    
    nx=L(4);
    
end


for x=1:nx
    
    Av=A(:,:,:,x);
    
    for a=1:oa-1
        
        Av(:,:,a)=ForceSDP(Av(:,:,a));
        
    end

    
    Av(:,:,oa)=eye(d)-sum(Av(:,:,1:oa-1),3);
    
    
    if min(eig(Av(:,:,oa))) < 0
        
        c = max(eig(sum(Av(:,:,1:oa-1),3)))^-1; 
        
        Av = c*Av; 
        
        Av(:,:,oa)=eye(d)-sum(Av(:,:,1:oa-1),3);
        
    end
        


    A(:,:,:,x)=Av;
    
end


f=A;

end
    
    
    