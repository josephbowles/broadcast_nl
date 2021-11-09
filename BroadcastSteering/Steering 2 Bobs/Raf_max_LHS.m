function f=Raf_max_LHS(sigrho,sigN,scenario,DB)

%Finds the largest q such that q*'prho'+(1-q)*'psig' is unsteerable in the Rafael scenario, characterized by 'scenario', and using vertices 'DB' for the  Bobs


%'sigrho' is a dA x dA x ob1 x ob2 x ny1 x ny2 matrix, interpreted as an assemblage in Rafael scenario
%'sigN' is a dA x dA x ob1 x ob2 x ny1 x ny2 matrix, interpreted as an assemblage in Rafael scenario
%'scenario' is of the form [ob1 ob2 ny1 ny2]
%'DB' is a L x (ob1*ob2*ny1*ny2) matrix, each line of which is interpreted as a distribution DB(k,:) = p(b1b2|y1y2) = [p(00|00) p(01|00) ... p(00|01) .. ]


%Rafael scenario

ob1=scenario(1); ny1=scenario(3);  ob2=scenario(2); ny2=scenario(4);

dA=size(sigrho,1); 

%convert vertices to a table
L=size(DB,1);

Df=zeros(ny1,ny2,ob1,ob2,L);

for k=1:L
    
    Df(:,:,:,:,k)=convP2Table(DB(k,:),[ny1 ny2 ob1 ob2]);
    
end



%%%%SDP

cvx_begin quiet

%variables:

variable q nonnegative
variable X(dA,dA,L) hermitian semidefinite;
dual variable F

% for k=1:L
%     
%     sig(:,:,k) == semidefinite(dA,dA);
%     
% end




q <= 1;


%LHS decomposition

  S=cvx(zeros(dA,dA,ob1,ob2,ny1,ny2));
for y1=1:ny1
    for y2=1:ny2
        
        for b1=1:ob1
            for b2=1:ob2
                
  
                for k=1:L
               S(:,:,b1,b2,y1,y2) = S(:,:,b1,b2,y1,y2) + X(:,:,k)*Df(y1,y2,b1,b2,k);
  %                  S = S + X(:,:,k)*Df(y1,y2,b1,b2,k);
                end
                
                
%             F(:,:,b1,b2,y1,y2):   S ==  q*sigrho(:,:,b1,b2,y1,y2) + (1-q)*sigN(:,:,b1,b2,y1,y2);
%             F(:,:,b1,b2,y1,y2): S == q*sigrho(:,:,b1,b2,y1,y2) + (1-q)*sigN(:,:,b1,b2,y1,y2) ;

            end
        end
        
        
    end
end



% 
F: S == q*sigrho + (1-q)*sigN  ;



% %explicit local decomposition
%
% prhoM=convP2Table(prho,[nx ny1 ny2 oa ob1 ob2]);
% psigM=convP2Table(psig,[nx ny1 ny2 oa ob1 ob2]);
%
%
% for x=1:nx
%     for y1=1:ny1
%         for y2=1:ny2
%
%
%             for a=1:oa
%                 for b1=1:ob1
%                     for b2=1:ob2
%
%                         S=0;
%                         for l=1:Lam
%
%                             DAm=convP2Table(DA(l,:),[nx oa]);
%                             S = S + DAm(x,a)*pNS(b1,b2,y1,y2,l) ;
%
%                         end
%
%                         R = R + (S == q*prhoM(x,y1,y2,a,b1,b2) + (1-q)*psigM(x,y1,y2,a,b1,b2));
%
%                     end
%                 end
%             end
%
%
%         end
%     end
% end



%objective function

O = q;



maximize(O);

cvx_end



f={O,F};


end