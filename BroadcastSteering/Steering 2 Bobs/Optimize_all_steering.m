

%Finds the smallest visibility 'v' for which v*rho+(1-v)*sig is steerable in Rafael scenario


%Rafael scenario

ob1=2; ny1=2;  ob2=2; ny2=2;   DB=VNS33;  %Insert all the vertices

% %3 outputs:
% oa=3; nx=3;   ob1=3; ny1=3;  ob2=3; ny2=3;

scenario=[ob1 ob2 ny1 ny2];




%initial state


% rho=MMM(1,pi/8);
% sig=eye(4)/4;
% dA=2; dB0=size(rho,1)/dA; dB1=dB0; dB2=dB0;

% 
% rho=HorodeckiState(1/2)
% sig=HorodeckiState(0);
% dA=3; dB0=size(rho,1)/dA; dB1=dB0; dB2=dB0;


d=2;
rho=IsotropicState(d,1);
sig=eye(4)/4;
dA=d; dB0=size(rho,1)/dA; dB1=dB0; dB2=dB0;

% 
% d=3;
% rho=WernerState(d,1);
% sig=WernerState(d,0);
% dA=d; dB0=size(rho,1)/dA; dB1=dB0; dB2=dB0;


% rho=MaxEntangled(d)*MaxEntangled(d)';
% sig=eye(d^2)/d^2;
% dA=d; dB0=size(rho,1)/dA; dB1=dB0; dB2=dB0;


% rho=rand_stateNL;
% dA=2; dB0=size(rho,1)/dA; dB1=dB0; dB2=dB0;
% sig=eye(dA*dB0)/(dA*dB0)^2;

% rho=RandomDensityMatrix(6);
% dA=3; dB0=size(rho,1)/dA; dB1=dB0; dB2=dB0;
% sig=eye(dA*dB0)/(dA*dB0)^2;


%creates the maximally mixed assemblage:

sigId=zeros(dA,dA,ob1,ob2,ny1,ny2);
for y1=1:ny1
    for y2=1:ny2
        
        
        for b1=1:ob1
            for b2=1:ob2
                
                
                
                sigId(:,:,b1,b2,y1,y2)=eye(dA)/(dA*ob1*ob2);
                
            end
        end
        
    end
end



%first finds a point outside:

vis=2;
while vis >= 1/sqrt(3)
    
    %
    %     B1=rand_povms(dB1,ob1,ny1);
    %
    %     B2=rand_povms(dB2,ob2,ny2);
    
    
  U=RandomUnitary(dB1*dB2);
    T=Choijoe(U,dB1);

%      T=RandomChoi_proj(dB0,dB1*dB2); yalmip clear
     %T=Tjoe;
      %T=Tfla;
   
      
      T = GHZChoi, 
    
    B1=rand_povms(dB1,ob1,ny1,1);
    %B1=proj_from_bloch(eye(3));
    B1=proj_from_bloch([1 0 0;0 0 1;0 1 0]);
    
    B2=rand_povms(dB2,ob2,ny2,1);
B2=proj_from_bloch([1 0 0;0 0 1;0 1 0]);  

U = RandomUnitary(2); for y2=1:3, for b2=1:2, B2(:,:,b2,y2) = U*B2(:,:,b2,y2)*U' ; end, end,

    sigrho=Rafael_assemblage(rho,dA,T,B1,B2);
    

    sigN=Rafael_assemblage(sig,dA,T,B1,B2);
    

   q=1-10^-8; q=1; f=Raf_max_LHS(sigrho,q*sigN+(1-q)*sigId,scenario,DB)
    cvx_clear
    

   
%      cvx_clear
    
    vis=f{1};
    

    
    if isnan(vis)

       
        vis=2;
        
    end
    
end

In=f{2}; In=make_hermitian(In);

q=vis; pq = q*sigrho + (1-q)*sigN;

S=0;
for y1=1:ny1
    for y2=1:ny2
        
        for b1=1:ob1
            for b2=1:ob2
                
  S = S + In(:,:,b1,b2,y1,y2)*pq(:,:,b1,b2,y1,y2); 
                
            end
        end
        
        
    end
end

sc=trace(S),



%  pause


visc=vis;

rep=0;
rep(1)=vis;

prec=10^-5;

gap=1; k=2;
while gap > prec
    
    rhoc=visc*rho+(1-visc)*sig;
    
    init={T,B2};
    
    viol=sc-1; count=0; limcount=7;
    while viol <= sc
        
        if count == 3
            f=seesaw_steering_inequ(In,rhoc,dA,scenario,dB1,dB2,2,init);
        else
            f=seesaw_steering_inequ(In,rhoc,dA,scenario,dB1,dB2,2);
        end
        % f=seesaw_inequ_without_T(In,rhoc,scenario,T);
        
        yalmip clear
        
        viol=f{1}
        
        count=count+1;
        if count > limcount
            break
        end
        
    end
    count
    

    
    
    T=ForceSDP(f{2}); B1=make_povms_valid(f{3}); B2=make_povms_valid(f{4});
    
   
    
    

    sigrho=Rafael_assemblage(rho,dA,T,B1,B2);
    

    sigN=Rafael_assemblage(sig,dA,T,B1,B2);

    
    q=1-10^-8; q=1; f=Raf_max_LHS(sigrho,q*sigN+(1-q)*sigId,scenario,DB)
    cvx_clear
    
    vis=f{1},
    
    
    % if vis < 1/sqrt(2)
    %
    %     return
    %
    % end
    
    rep(k)=vis;
    k=k+1;
    
    %new inequality
    In=f{2}; In=make_hermitian(In);
    
    %compute current score with new inequ
    
    q=vis; pq = q*sigrho + (1-q)*sigN;

S=0;
for y1=1:ny1
    for y2=1:ny2
        
        for b1=1:ob1
            for b2=1:ob2
                
  S = S + In(:,:,b1,b2,y1,y2)*pq(:,:,b1,b2,y1,y2); 
                
            end
        end
        
        
    end
end

sc=trace(S),
    
    gap=visc-vis;
    
    visc=vis;
    
end


