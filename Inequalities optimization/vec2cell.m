function f = vec2cell(v)

%Converts a vector 'v' into a cell 


C=cell(size(v));

for k=1:length(v)
    
    C{k}=v(k);
    
end


f=C;

end