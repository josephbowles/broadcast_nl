function f=concatenate_prob(p1,p2,na,nb)

%Merge two distributions 'p1'=pA(a|x) and 'p2'=pB(b|y) into a the global distribution p = [ p(00|00) p(01|00) ... p(na nb |000) p(00|01) ... p(na nb | ma mb ) ], where 'na', 'nb' are the respective numbers of outputs 

% 'p1' and 'p2' are assumed to be vectors of the form p = [ p(0|0) p(1|0) ... p(n|0) p(0|1) ... p(n|m) ]


ma=length(p1)/na; 

mb=length(p2)/nb;

M1=reshape(p1,na,ma);

M2=reshape(p2,nb,mb);

f=vec(kron(M1,M2));

