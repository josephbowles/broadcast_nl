function f=make_hermitian(A)

%Makes each operator in 'A' hermitian (outputs 1/2*(A(:,:,x)+A(:,:,x)') for all x )

%'A' is a d x d x n matrix, where 'n' can be multi-dimensional 



L=size(A);

%if 'A' contains a single matrix, artificially interpret it as a set of matrices

if length(L) < 3
    L=[size(A) 1];
end

%if 'A' is labelled with a multi-dimensional label, gather everything into a one-dimensional one
if length(L) > 3
    
    A=reshape(A,[L(1) L(2) prod(L(3:end))]); 
    
end


for x=1:size(A,3)
    
    A(:,:,x) = (A(:,:,x)+A(:,:,x)')/2;
    
end


A=reshape(A,L); 

f=A;

end