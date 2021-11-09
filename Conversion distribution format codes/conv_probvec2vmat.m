function f=conv_probvec2vmat(p,oa,ob,nx) 

%Converts a probability vector to a matrix of bipartite distributions 'M' (Alice has 'oa' outputs and 'nx' inputs, Bob 'ob' outputs)

%Input format: 'p' is a 1x(oa*ob*nx*ny) vector, interpreted as a bipartite distribution p = p(ab|xy) = [p(00|00) p(01|00) .. p(0 ob|00) p(10|00) ... p(oa ob|00) p(00|01) ... ... p(oa ob|nx ny) ]

%Output format: 'M' is a (oa*nx) x (ob*ny) matrix, where each oa x ob block contains the probability distribution associated to the inputs x and y (that is, p(ab|xy) is located on row (x-1)*a+a and column (y-1)*b+b) 



ny=length(p)/(oa*ob*nx);

M=zeros(oa*nx,ob*ny);

%loop over blocks: 

for x=1:nx
    for y=1:ny
        
        
        in=oa*ob*( (x-1)*ny + (y-1))+1; fin=in+oa*ob-1;
        
        v=p(in:fin);
        
        v=reshape(v,ob,oa);
        
        M((x-1)*oa+1:x*oa,(y-1)*ob+1:y*ob)=v';
        
        
    end
end

f=M;

end

        