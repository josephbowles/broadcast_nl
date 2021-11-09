function f=update_vector(v,Maxval,LeftOrRight)

%Updates vector 'v' to its next allowed value, where each component maximum value is given by 'Maxval'. By default changes last entry first, flips this conventions if LeftOrRight is set to 2

if nargin < 3
    
    LeftOrRight=1;
    
end

N=length(v);


if LeftOrRight==1
    
    
    for k=N:-1:1
        
        
        if v(k) < Maxval(k)
    
            v(k) = v(k)+1;
            v(k+1:end)=1;
            break
        end
        
        
    end
    
    
elseif LeftOrRight==2
    
    for k=1:N
        
        if v(k) < Maxval(k)
            
            v(k) = v(k)+1;
            v(1:k-1)=1;
            break
        end
        
        
    end
    
end



f=v;


end