function f=convBiprob2Table(p,oa,ob,nx,ny)

%Converts a collection of bipartite distributions 'p' (in scenario [oa,ob,nx,ny]) from the vector form v=[p(00|00) p(01|00) ... ] to a table T(a,b,x,y)

%Input format: 'p' is a  n x (oa*ob*nx*ny) matrix, each line of which is interpreted as a bipartite distribution p(k,:) = p(ab|xy) = [p(00|00) p(01|00) .. p(0 ob|00) p(10|00) ... p(oa ob|00) p(00|01) ... ... p(oa ob|nx ny) ]

%Output format: 'f' is a oa x ob x nx x ny x n matrix, interpreted as collection of distribution T(a,b,x,y,k) (for k=1:n)

R=zeros(oa,ob,nx,ny,size(p,1));

for k=1:size(p,1)

    
vt=convP2Table(p,[nx,ny,oa,ob]);

R(:,:,:,:,k)=permute(vt,[3 4 1 2]);

end


f=R;

end
