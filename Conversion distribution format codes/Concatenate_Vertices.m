function f=Concatenate_Vertices(DA,DB,oa,ob)

%Concatenates distributions 'DA' and 'DB' into a single distribution, where 'DA' is a distribution with 'oa' outcomes, 'DB' with 'ob' outcomes

%Input format: 'DA' and 'DB' are  n x (oa*ob*nx*ny) matrices, each line of which is interpreted as a bipartite distribution p = p(ab|xy) = [p(00|00) p(01|00) .. p(0 ob|00) p(10|00) ... p(oa ob|00) p(00|01) ... ... p(oa ob|nx ny) ]


 

LA=size(DA);
LB=size(DB);




V=zeros(LA(1)*LB(1),LA(2)*LB(2));
n=1;
for k=1:LA(1)
    for l=1:LB(1)
        
%         k, l 
%         
%         size(DA(k,:)), size(DB(l,:)), 
        
        V(n,:)=concatenate_prob(DA(k,:),DB(l,:),oa,ob);
        
%         f=V(1,:); return
%         
        n=n+1;
        
    end
end

f=V;

end
        