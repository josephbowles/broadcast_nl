function f=convP2Table(p,dims) 

%Converts a vector p(abc...) of the form [p(000..0) p(00...1) ... ] to a multi-dimensional matrix p(a,b,c,...)

%'p' is a 1xd vector 

%'dims' is a 1xn vector of integers, containing the dimensions of a,b,c,etc. More precisely, dims=[da db dc ... ]. One needs prod(dims)=d


T=reshape(p,flip(dims));

f=permute(T,flip(1:length(dims)));


end