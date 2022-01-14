function f=Rafael_pointAABB(rho,dA0,TA,A1,A2,TB,B1,B2,UA12,UB12)

%In Rafael scenario with 2 Alices and 2 Bobs, computes the overall distribution p(a1a2b1b2|x1x2y1y2) coming from state 'rho' (of local Alice dimension 'dA'), isometries 'TA' and 'TB', and measurements 'A1', 'A2', 'B1' and 'B2'
%(optional: Unitary UA12 on the output space of Alice, Unitary UB12 on the output space of Bob)


%The output is given in the usual convention p(a1a2b1b2|x1x2y1y2) = [ p(0000|0000) p(0001|0000) ... ]



%extracts scenario

oa1=size(A1,3); nx1=size(A1,4);

oa2=size(A2,3); nx2=size(A2,4);

ob1=size(B1,3); ny1=size(B1,4);

ob2=size(B2,3); ny2=size(B2,4);


dA1=size(A1,1);
dA2=size(A2,1);
dB0=length(rho)/dA0;
dB1=size(B1,1);
dB2=size(B2,2);


if nargin < 9
    
    UA12=eye(dA1*dA2);
    UB12=eye(dB1*dB2);
    
end

%applies isometry

rhot=Apply_ChoiAB(rho,TA,TB,dA0,UA12,UB12);

% dA1=size(A1,1);
% dA2=size(A2,1);
% dB0=length(rho)/dA0;
% dB1=size(B1,1);
% dB2=size(B2,2);
% %Bob
% TOT=kron(PartialTranspose(rho,2),eye(dB1*dB2))*kron(eye(dA0),TB);
%
% rhotb=PartialTrace(TOT,2,[dA0 dB0 dB1 dB2]);
%
%
% %Alice
% dBt=dB1*dB2; dAfin=dA1*dA2;
% rhos=Swap(rhotb,[1,2],[dA0 dBt]);
%
% TOT=kron(PartialTranspose(rhos,2),eye(dAfin))*kron(eye(dBt),TA);
%
% rhots=PartialTrace(TOT,2,[dBt dA0 dAfin]);
%
% rhot=Swap(rhots,[1,2],[dBt dAfin]);


%computes the probability distribution

p=zeros(nx1,nx2,ny1,ny2,oa1,oa2,ob1,ob2);

for x1=1:nx1
    for x2=1:nx2
        for y1=1:ny1
            for y2=1:ny2
                
                for a1=1:oa1
                    for a2=1:oa2
                        for b1=1:ob1
                            for b2=1:ob2
                                
                                
                                
                                p(x1,x2,y1,y2,a1,a2,b1,b2)=trace(superkron(A1(:,:,a1,x1),A2(:,:,a2,x2),B1(:,:,b1,y1),B2(:,:,b2,y2))*rhot);
                                
                                
                                if (p(x1,x2,y1,y2,a1,a2,b1,b2) ) < -10^-8
                                    
                                    disp('bug, negative probability')
                                    
                                    return
                                    
                                end
                                
                                if imag(p(x1,x2,y1,y2,a1,a2,b1,b2) ) > 10^-8
                                    
                                    disp('bug, imaginary probability')
                                    
                                    return
                                    
                                end
                                
                            end
                        end
                        
                    end
                    
                end
            end
            
        end
        
    end
    
end




f=real(convTable2P(real(p),[nx1 nx2 ny1 ny2 oa1 oa2 ob1 ob2]));

end