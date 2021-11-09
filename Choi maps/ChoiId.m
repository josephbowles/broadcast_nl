function f = ChoiId(d)

%Generates the Choi state of the idendity map of dimension 'd'

PHIP=MaxEntangled(d)*MaxEntangled(d)'; 

f = d*PHIP; 

end