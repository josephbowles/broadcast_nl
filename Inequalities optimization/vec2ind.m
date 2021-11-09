function LinInd = vec2ind(siz,vec)

%Converts vector 'vec' into correspoding linear index, for a matrix of size 'siz'

%That is, if LinInd = vec2ind(size(A),vec), then A(LinInd) = A(vec(1),vec(2),...)

%'siz' is a 1xn vector of integers
%'vec' is a 1xn vectors of integers


vecC=vec2cell(vec); 

LinInd=sub2ind(siz,vecC{:});

end