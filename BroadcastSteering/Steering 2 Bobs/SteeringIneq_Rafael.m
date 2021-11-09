function f=SteeringIneq_Rafael(In,sig)

%Evaluates the steering inequality 'In' on Rafael assemblage 'sig' 

%'ineq' is a dA x dA x ob1 x ob2 x ny1 x ny2 matrix, interpreted as the steering inequality matrices I_{ob1,ob2,ny1,ny2}
%'sig' is a dA x dA x ob1 x ob2 x ny1 x ny2 matrix, interpreted as an assemblage in Rafael scenario


L=size(In); ob1=L(3); ob2=L(4); ny1=L(5); ny2=L(6); 

S=0;
for y1=1:ny1
    for y2=1:ny2
        
        for b1=1:ob1
            for b2=1:ob2
                
  S = S + In(:,:,b1,b2,y1,y2)*sig(:,:,b1,b2,y1,y2); 
                
            end
        end
        
        
    end
end


f=trace(S); 

end