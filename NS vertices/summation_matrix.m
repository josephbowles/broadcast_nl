function f=summation_matrix(dims,sum_index)

%Generates the matrix A such that A*p' corresponds to a summation on p over the variable 'sum_index' (p is interpreted as p(abc...) with the usual convention p=[p(00..00) p(00..01) ... p(100...0) ... ])


%'dims' is a 1xd vector of integers, containing the dimensions of a,b,c,etc. More precisely, dims=[da db dc ... ]. One needs prod(dims)=d
%'sum_index' is an integer




n=sum_index;


T=convP2Table(zeros(1,prod(dims)),dims);


%dimension of the matrix:

L=prod(dims)/dims(n);

M=zeros(L,prod(dims));


dimsum=dims; dimsum(n)=[];

vecloop=loop_vector(dimsum); 



for k=1:length(vecloop)
    
    v=vecloop(k,:);
    
    Cin=num2cell(v(1:n-1));
    Cfin=num2cell(v(n:end));
    
    T(:)=0;
    T(Cin{:},:,Cfin{:})=1;
    
    M(k,:)=convTable2P(T,dims);
    
end


f=M;

% 
% %check the result:
% T=convP2Table(p,dims); 
% norm(unique(reshape(sum_n(T,n),1,L))-unique(f*p')')

end
    

