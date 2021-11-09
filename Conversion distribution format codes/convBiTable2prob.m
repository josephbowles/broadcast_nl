function f=convBiTable2prob(T)

%Converts a collection of bipartite distributions 'T' (in the table form T(a,b,x,y,k)) to the vector form v=[p(00|00) p(01|00) ... ] 

%Input format: 'T' is a oa x ob x nx x ny x n matrix, interpreted as collection of distribution T(a,b,x,y,k) (for k=1:n)

%Output format: 'p' is a  n x (oa*ob*nx*ny) matrix, each line of which is interpreted as a bipartite distribution p(k,:) = p(ab|xy) = [p(00|00) p(01|00) .. p(0 ob|00) p(10|00) ... p(oa ob|00) p(00|01) ... ... p(oa ob|nx ny) ]


L=size(T);

R=zeros(L(5),prod(L(1:4)));


for k=1:L(5)

vt=permute(T(:,:,:,:,k),[3 4 1 2]);

R(k,:)=convTable2P(vt,[L(3) L(4) L(1) L(2)]);

end


f=R;

end
