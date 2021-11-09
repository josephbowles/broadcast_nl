function f=Rafael_assemblage(rho,dA,T,B1,B2)

%In Rafael scenario, computes the overall assemblage sig(b1b2|y1y2) coming from state 'rho' (of local dimension dA), Choi state 'T', and measurements, 'B1' and 'B2'


%extracts scenario

ob1=size(B1,3); ny1=size(B1,4);

ob2=size(B2,3); ny2=size(B2,4);


dB0=length(rho)/dA;
dB1=size(B1,1);
dB2=size(B2,2);





%applies isometry

TOT=kron(PartialTranspose(rho,2,[dA dB0]),eye(dB1*dB2))*kron(eye(dA),T);

rhot=PartialTrace(TOT,2,[dA dB0 dB1 dB2]);



%computes the assemblage

sig=zeros(dA,dA,ob1,ob2,ny1,ny2);



for y1=1:ny1
    for y2=1:ny2        
        
        for b1=1:ob1
            for b2=1:ob2

                
                sig(:,:,b1,b2,y1,y2)=PartialTrace(Tensor(eye(dA),B1(:,:,b1,y1),B2(:,:,b2,y2))*rhot,[2 3],[dA dB1 dB2]);
                
                
            end
        end
            
    end
end






f=sig;

end