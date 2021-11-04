function corr = ToCorrelatorNotationWithInOut(settings,outputs)
    % TODO make it so I don't define x,y,z manually but variable
    % depending on length(settings)
    x=settings(1);
    y=settings(2);
    z=settings(3);
    
    Ax = join(['A',string(x),''],''); 
    Ax = sym(char(Ax));
    By = join(['B',string(y),''],''); 
    By = sym(char(By));
    Cz = join(['C',string(z),''],''); 
    Cz = sym(char(Cz));
    
    a=outputs(1)-1; % so we have {0,1} instead of {1,2}
    b=outputs(2)-1;
    c=outputs(3)-1;

    %if a>3 || b>3 || c>3
    %   printf("Doesn't make sense to use the correlator notation for more than 2 outputs in this case.") 
    %end
    
    corr = 1/8 * (1 + ... 
        ((-1)^a)*Ax        + ((-1)^b)*By        + ((-1)^c)*Cz + ...
        ((-1)^(a+b))*Ax*By + ((-1)^(a+c))*Ax*Cz + ((-1)^(b+c))*By*Cz + ...
        ((-1)^(a+b+c))*Ax*By*Cz );
end