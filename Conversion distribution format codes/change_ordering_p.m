function f=change_ordering_p(p,dims,order)

%Change the ordering of a each line of 'p', where a line is of the form p(abc...) (that is, written as [p(00..0) p(00..1) ... ]) following the permutation 'order'


%'p' is a nxd matrix

%'dims' is a 1xN vector of integers, containing the dimensions of a,b,c,etc. More precisely, dims=[da db dc ... ]. One needs prod(dims)=d

%'order' is a 1xN vector of permutation (ie, a vector containing all integers from 1 to n, in an arbitrary order)

pf=zeros(size(p,1),prod(dims));
for k=1:size(p,1)

    v=p(k,:);
    
D=convP2Table(v,dims);

D=permute(D,order); 

perdims=dims(order);

pf(k,:)=convTable2P(D,perdims);

end

f=pf;

end