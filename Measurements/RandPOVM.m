function f=RandPOVM(d,oa) 

%Generates a random 'oa'-outcome POVM in dimension 'd' (default: oa = d) 


if nargin < 2 
    oa=d;
end



A=RandomPOVM(d,oa); 

S=zeros(d,d,oa);

for k=1:oa 
    
    S(:,:,k)=A{k};
    
end


f=S;

end