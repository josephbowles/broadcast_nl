function f=RandomChoi_proj(d1,d2)

%Creates a random projective Choi state of dimension 'd1'*'d2' 

%More precisely, creates a d1*d2 x d1*d2 positive semi-definite matrix, with PartialTrace(T,2,[d1 d2])=eye(d1) 



%sdp

T=sdpvar(d1*d2,d1*d2,'hermitian','complex'); 

R = (T >= 0) + (PartialTrace(T,2,[d1 d2]) == eye(d1)); 


%random objective

C=Normalize(normrnd(0,1,d1*d2,1)+1i*normrnd(0,1,d1*d2,1)); C=C*C'; 

O =real( trace(C*T));



 ops=sdpsettings('verbose',0,'warning',0,'solver','mosek');
 optimize(R,-O,ops);
 
 
 f=double(T);
 
end

