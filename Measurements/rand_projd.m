function f=rand_projd(d,oa)

%creates a random 'oa'-outcome projective measurement of dimension 'd'
%If not specified, oa=d and all elements are rank-1


A=zeros(d,d,d);

%computational basis

for k=1:d
    
    A(k,k,k)=1;
    
end

%applies a random unitary on each element

U=RandomUnitary(d);

for k=1:d
    
    A(:,:,k)=U*A(:,:,k)*U';
    
end


if nargin == 2
    
    r=mod(d,oa);
    
    n=(d-r)/oa;
    
    
    B=zeros(d,d,oa);
    
 

    for k=1:r
 
        B(:,:,k)=sum(A(:,:,(k-1)*(n+1)+1:k*(n+1)),3);
        
    end
    
    if r==0
        fin=0;
    else
    fin=k*(n+1);
    end
    
    for k=r+1:oa
        in=fin+1;
        fin=in+n-1;

        B(:,:,k)=sum(A(:,:,in:fin),3);
        
    end
    
    A=B;
end


f=A;