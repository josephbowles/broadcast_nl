function f=rand_povms(d,na,nx,proj)

%Creates a random set of d-dimensional POVMs labeled by 'nx' (each POVM containing 'na' outcomes), (all projectors if 'proj' is set )

%'d' is an integer 
%'na' is an integer 
%'nx' is a vector of integers
%'proj' is a 0/1 switch



nxeff=prod(nx);


A=zeros(d,d,na,nxeff);


for x=1:nxeff
    
    
    if nargin < 4
        
    A(:,:,:,x)=RandPOVM(d,na);
    
    
    else
        
        A(:,:,:,x)=rand_projd(d,na);
        
    end
    
    
end




vr=[d d na nx ]; 

%L=length(nx);
% for k=1:L
%     vr(k+3)=nx(k);
% end


f=reshape(A,vr);

end
    

    
