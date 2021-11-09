function f = Matgen_Dax(Oa,Ia)

%Generate the D_\lambda(a|x)



nA=Ia;
mA=Oa;
Ndet = mA^nA; 
Idm = eye(mA);
SingleParty = zeros(nA*mA,Ndet);
multf = linspace(nA-1,0,nA); %used for 'binary' filling of arrays
multf = mA.^multf;


%% SINCE WE WANT VARIABLE NUMBER OF MEASUREMENTS USE ODOMETER

odo_n = nA;             % number of for loops
odo_k = mA*ones(1,nA);  % number of iterations inside each for loop in order
odo_inc = ones(1,nA);   % increment value for each for loop
odo_sp = zeros(1,nA);   %starting position of odometer

odo = odo_sp;
exitflat = 0;
while exitflat < 1
    
    for i = odo_n:-1:1
        if odo(i)+odo_inc(i) <= odo_sp+(odo_k(i)-1)*odo_inc(i)
            odo(i) = odo(i) + odo_inc(i);
            break
        else
            odo(i) = odo_sp(i);
        end
        
        if i == 1
            exitflat = 1;
        end
    end
    
    % ---- inside of nested loops starts here
    
    temp = [];
    for j = 1:nA
        temp = [temp; Idm(:,odo(j)+1)];
    end
    SingleParty(:,1+multf*odo') = temp;
 
    
    % ---- inside of nested loops ends here
end

f=reshape(SingleParty,Oa,Ia,Ndet);




end