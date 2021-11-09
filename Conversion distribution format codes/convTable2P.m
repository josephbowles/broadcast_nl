function f=convTable2P(p,dims)

%Converts a multi-dimensional matrix p(a,b,c,...) to a vector p(abc...) of the form [p(000..0) p(000..1) ... ] 

%'p' is a na x nb x nc x .. arbitrary matrix

%'dims' is a 1xn vector of integers, containing the dimensions of a,b,c,etc. More precisely, dims=[da db dc ... ]. One needs prod(dims)=d


p=permute(p,flip(1:length(dims)));


f=reshape(p,1,numel(p));

end