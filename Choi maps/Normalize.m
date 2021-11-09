function v=Normalize(x)

%Normalizes each line of the matrix 'x' (or the vector if a column-vector is inputed)

%'x' is a n x m matrix, each line of which is interpreted as a vector in R^m

v=x;

if size(x,2) ~= 1
    
for l=1:size(x,1)
    
v(l,:)=x(l,:)/norm(x(l,:));

end

else
    
    v=x/norm(x);
end

end