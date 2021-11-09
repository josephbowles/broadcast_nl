

%Generates the no-signalling vertices from the no-signalling facets 

%The ouptut 'VNS' is understood as a  n x (oa*ob*nx*ny) matrix, each line of which is interpreted as a bipartite distribution p = p(ab|xy) = [p(00|00) p(01|00) .. p(0 ob|00) p(10|00) ... p(oa ob|00) p(00|01) ... ... p(oa ob|nx ny) ]



%scenario:

oa=2; ob=2; nx=2; ny=2;

L=oa*ob*nx*ny;




%Generates no-signalling matrix constraints: 


%Alice's no-signalling:

F=summation_matrix([nx ny oa ob],3);

block=ob*ny;
SX=zeros(block,L,nx-1);

FNSAB=zeros(1,L);
for x=1:nx-1
    
    SX(:,:,x) = F((x-1)*block+1:x*block,:) - F(x*block+1:(x+1)*block,:);
    
%     SX(:,:,x)*pQ

FNSAB=[FNSAB;SX(:,:,x)];
    
end

FNSAB=FNSAB(2:end,:);





%Bob's no-signalling:

F=summation_matrix([nx ny oa ob],4);

Fp=zeros(size(F));
for k=1:size(F,2)
    
    Fp(:,k)=change_ordering_p(F(:,k),[nx ny oa],[2 1 3]);
    
end



block=oa*nx;
SY=zeros(block,L,ny-1);

FNSBA=zeros(1,L);
for y=1:ny-1

    
    M = Fp((y-1)*block+1:y*block,:) - Fp(y*block+1:(y+1)*block,:);
   

    
    R=zeros(size(M));
for k=1:size(F,2)
    
    R(:,k)=change_ordering_p(M(:,k),[1 nx oa],[2 1 3]);
    
end


SY(:,:,y) = R;

FNSBA=[FNSBA; R];

%     SY(:,:,y)*pQ
    
end

FNSBA=FNSBA(2:end,:);

    


%Normalization matrix    


FN=zeros(nx*ny,L); 

k=1;
for x=1:nx
    for y=1:ny
        
        in=(x-1)*ny*oa*ob+(y-1)*oa*ob+1;
        fin=in+oa*ob-1;
        v=zeros(1,L); v(in:fin)=1;
        
        FN(k,:)=v;
        
        k=k+1;
        
    end
end
    


%Positivity matrix 

FP=eye(L,L);



%Overall equality constraints:

Feq=[FN;FNSAB;FNSBA];

beq=zeros(size(Feq,1),1); beq(1:nx*ny)=1;

    

%Inequality constraints:

Fineq=FP;

bineq=zeros(size(FP,1),1);



%Gives the vertices

V=lcon2vert(-Fineq,bineq,Feq,beq); 

VNS=1/2*round(2*V);

    
    
