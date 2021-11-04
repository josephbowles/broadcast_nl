function channel = giveChannelRAND(dimx,dimy)
    id = eye(dimx);
    
    e = cell(dimx);
    for i=1:dimx
       e{i} = id(:,i); 
    end
    
    o = cell(dimy);
    for i = 1:dimy
       o{i}=zeros(dimy,1);
    end
    
    o{1} = (-1+2*rand(dimy,1))+1i*(-1+2*rand(dimy,1));
    o{1} = o{1}/norm(o{1});
    
    channel = o{1} * e{1}';
    for i=2:dimx
       randvec = (-1+2*rand(dimy,1))+1i*(-1+2*rand(dimy,1));
       proyector = 0;
       for j=1:dimy
          proyector = proyector + o{j}*o{j}';
       end
       o{i} = randvec - proyector * randvec;
       o{i} = o{i} / norm(o{i});
       
       channel = channel + o{i} * e{i}';
    end

    % check it's an isometry
    diff = channel' * channel - eye(dimx);
    assert(norm(diff)<1e-8, 'Channel is not an isometry within the specified precision.');
    
end

