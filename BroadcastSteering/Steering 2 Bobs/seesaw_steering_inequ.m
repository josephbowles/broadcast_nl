function f=seesaw_steering_inequ(ineq,rho,dA,scenario,dB1,dB2,choiin,Initial) 

%In Rafael scenario (with Bobs' dimensions 'dB1', 'dB2') optimizes steering inequality 'ineq' for state 'rho' (of local dim 'dA'), in scenario 'scenario' (starting from optional points 'Initial')

%'ineq' is a dA x dA x ob1 x ob2 x ny1 x ny2 matrix, interpreted as the steering inequality matrices I_{ob1,ob2,ny1,ny2}
%'rho' is a bipartite state
%'scenario' is of the form [ob1 ob2 ny1 ny2] 
%'choiin' is an integer, wich tells whether the initial Choi state is taken as a projector (if choiin==1), ChoiJoe style (if choiin==2), or a generic Choi state (other values)

%optional variable:
%'Initial' is a 1x2 cell, containing initials values of T, B2. That is, Initial = {T,B2}


%Rafael scenario

ob1=scenario(1); ny1=scenario(3);  ob2=scenario(2); ny2=scenario(4);


dB0=size(rho,1)/dA;


It=ineq;

%initial point


if nargin < 8
    
    if choiin == 1
        
        Tc=RandomChoi_proj(dB0,dB1*dB2);
        
    elseif choiin == 2
        
        U=RandomUnitary(dB1*dB2);
        Tc=Choijoe(U,dB1);
        
    else
        
        Tc=RandomChoi(dB0,dB1*dB2);
        
    end
    
    
    
    B2c=rand_povms(dB2,ob2,ny2,1);
    
else
    
    Tc=Initial{1};  B2c=Initial{2};
    
end





prec=10^-5;

former=1i;
c=1;
while (c > prec)


rr=optineqB1(Tc,B2c);

B1c=rr{2};

rr=optineqB2(Tc,B1c);

B2c=rr{2};

rr=optineqT(B1c,B2c);

Tc=rr{2};

c = abs(rr{1} - former);

former=rr{1};

end



f={rr{1},Tc,B1c,B2c};



%%%sub-functions: 




function y=optineqB1(T,B2) 

    

%applies isometry

TOT=kron(PartialTranspose(rho,2,[dA dB0]),eye(dB1*dB2))*kron(eye(dA),T); 

rhot=PartialTrace(TOT,2,[dA dB0 dB1 dB2]);

%sdp variables

B1=sdpvar(dB1,dB1,ob1,ny1,'hermitian','complex');

R=[];

for x=1:ny1
    for a=1:ob1
        
        R = R + (B1(:,:,a,x) >= 0);
        
    end
    
    R = R + (sum(B1(:,:,:,x),3)==eye(dB1) ) ;
    
end



S=0;


for y1=1:ny1
    for y2=1:ny2
        
        
        for b1=1:ob1
            for b2=1:ob2
                
                
                
                S = S + trace(It(:,:,b1,b2,y1,y2) * PartialTrace(Tensor(eye(dA),B1(:,:,b1,y1),B2(:,:,b2,y2))*rhot,[2 3],[dA dB1 dB2]));
                
                
            end
        end
        
    end
end



O=real(S);



 ops=sdpsettings('verbose',0,'warning',0,'solver','mosek');
 optimize(R,-O,ops);
 
 
y={double(O),double(B1)};

end


function y=optineqB2(T,B1) 

    

%applies isometry

TOT=kron(PartialTranspose(rho,2,[dA dB0]),eye(dB1*dB2))*kron(eye(dA),T); 

rhot=PartialTrace(TOT,2,[dA dB0 dB1 dB2]);


%sdp variables

B2=sdpvar(dB2,dB2,ob2,ny2,'hermitian','complex');

R=[];

for x=1:ny2
    for a=1:ob2
        
        R = R + (B2(:,:,a,x) >= 0);
        
    end
    
    R = R + (sum(B2(:,:,:,x),3)==eye(dB2) ) ;
    
end


S=0; 

for y1=1:ny1
    for y2=1:ny2
        
        
        for b1=1:ob1
            for b2=1:ob2
                
                
                
                S = S + trace(It(:,:,b1,b2,y1,y2) * PartialTrace(Tensor(eye(dA),B1(:,:,b1,y1),B2(:,:,b2,y2))*rhot,[2 3],[dA dB1 dB2]));
                
                
            end
        end
        
    end
end



O=real(S);



 ops=sdpsettings('verbose',0,'warning',0,'solver','mosek');
 optimize(R,-O,ops);
 
 
y={double(O),double(B2)};

end


function y=optineqT(B1,B2) 


%isometry

T=sdpvar(dB0*dB1*dB2,dB0*dB1*dB2,'hermitian','complex');


R = (T >= 0) + (PartialTrace(T,2,[dB0 dB1*dB2]) == eye(dB0));



%applies isometry

TOT=kron(PartialTranspose(rho,2,[dA dB0]),eye(dB1*dB2))*kron(eye(dA),T); 

rhot=PartialTrace(TOT,2,[dA dB0 dB1 dB2]);

S=0; 

for y1=1:ny1
    for y2=1:ny2
        
        
        for b1=1:ob1
            for b2=1:ob2
                
                
                
                S = S + trace(It(:,:,b1,b2,y1,y2) * PartialTrace(Tensor(eye(dA),B1(:,:,b1,y1),B2(:,:,b2,y2))*rhot,[2 3],[dA dB1 dB2]));
                
                
            end
        end
        
    end
end



O=real(S);



 ops=sdpsettings('verbose',0,'warning',0,'solver','mosek');
 optimize(R,-O,ops);
 
 
y={double(O),double(T)};

end


end

