function f=conv_probmat2vec(M,oa,ob) 

%Converts a matrix of bipartite distributions 'M' to a single vector (Alice has 'oa' outputs, Bob 'ob')

%Input format: 'M' is a (oa*nx) x (ob*ny) matrix, where each oa x ob block contains the probability distribution associated to the inputs x and y (that is, p(ab|xy) is located on row (x-1)*a+a and column (y-1)*b+b) 

%Output format: 'p' is a 1x(oa*ob*nx*ny) vector, interpreted as a bipartite distribution p = p(ab|xy) = [p(00|00) p(01|00) .. p(0 ob|00) p(10|00) ... p(oa ob|00) p(00|01) ... ... p(oa ob|nx ny) ]


D=size(M);

nx=D(1)/oa;
ny=D(2)/ob;

p=zeros(1,oa*ob*nx*ny);

%loop over blocks: 

for x=1:nx
    for y=1:ny
        
        B=M((x-1)*oa+1:x*oa,(y-1)*ob+1:y*ob);
        
        v=reshape(B',1,oa*ob); 
        
        in=oa*ob*( (x-1)*ny + (y-1))+1; fin=in+oa*ob-1;
        
        p(in:fin)=v; 
        
    end
end

f=p;

end

        