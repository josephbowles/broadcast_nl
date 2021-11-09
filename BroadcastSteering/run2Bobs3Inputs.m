%Finds the smallest visibility 'v' for which v*rho+(1-v)*sig is steerable
%with 2 Bobs

%all_info contains a cell with all relevant information
%all_info{1,1} = all visibility values
%all_info{1,2} = sigma
%all_info{1,3} = scenario

%all_info{2,1} = eta
%all_info{2,2} = time4_evaluating_eta
%all_info{2,3} = B1
%all_info{2,4} = B2
%all_info{2,5} = T
%all_info{2,6} = Witness

clear
rng('shuffle')
C=clock;
clock_string=strcat(int2str(C(1)),'_',int2str(C(2)),'_',int2str(C(3)),'_',int2str(C(4)),'_',int2str(C(5)),'_',int2str(round(C(6))))
info='2Bobs3InputsQubitIstropic_';
file_name=strcat(info,clock_string);
if ~exist('2Bobs3Inputs', 'dir')
    mkdir 2Bobs3Inputs
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load VNS33;
DB=VNS33;  %Variable which contains the vertices of the scenario
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ob1=2; ny1=3;  ob2=2; ny2=3;   %number of inputs and outputs for each Bob
scenario=[ob1 ob2 ny1 ny2];

%initial state
d=2;
rho=IsotropicState(d,1);
sig=eye(d^2)/d^2;
dA=d; dB0=size(rho,1)/dA; dB1=dB0; dB2=dB0;
% END

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


all_info{2,1}=rho;
all_info{3,1}=sig;
all_info{4,1}=scenario;


%Prepare the initial assamblage

B1=rand_povms(dB1,ob1,ny1,1);  %Measurements for Bob1
B2=rand_povms(dB2,ob2,ny2,1);  %Measurments for Bob2
T=Choijoe(RandomUnitary(d^2),d);   %Choi state of the broadcast map (Alice to Bob1 and Bob2)
sigrho=Rafael_assemblage(rho,dA,T,B1,B2);  %Construct the ass 
sigN=Rafael_assemblage(sig,dA,T,B1,B2); %Construct the noise ass

%Evaluate the Robstness of sig_rho
q=1-10^-8; q=1; 
tic;
f=Raf_max_LHS(sigrho,q*sigN+(1-q)*sigId,scenario,DB);
time_SDP=toc;



vis=f{1};
visibility4initialpoint=vis

if vis>=0.9999  %If condition to ensure that a violation was found with the initial point
    disp('NO VIOLATION WAS FOUND!!! DANGER')
    pause
end
ineq=f{2};
ineq=make_hermitian(ineq);
In=ineq;

all_info{1,2}=vis;
all_info{2,2}=time_SDP;
all_info{3,2}=B1;
all_info{4,2}=B2;
all_info{5,2}=T;
all_info{6,2}=ineq;
cd 2Bobs3Inputs
save(file_name,'all_info')
cd ..


%Evaluate the unsteerable bound using a point in the boundary
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
sc=trace(S);

visc=vis;
rep=0;
rep(1)=vis;
prec=10^-5;
all_info{1,1}=rep;

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
        viol=f{1};
        count=count+1;
        if count > limcount
            break
        end
    end
    count
    T=get_valid_channel(f{2},2,4); B1=make_povms_valid(f{3}); B2=make_povms_valid(f{4});
    sigrho=Rafael_assemblage(rho,dA,T,B1,B2);
    sigN=Rafael_assemblage(sig,dA,T,B1,B2);
    q=1-10^-8; q=1; 
    tic
    f=Raf_max_LHS(sigrho,q*sigN+(1-q)*sigId,scenario,DB);
    time_SDP = toc;
    cvx_clear
    
    vis=f{1},
    rep(k)=vis;
    k=k+1;
    
    %new inequality
    In=f{2}; In=make_hermitian(In);

all_info{1,1}=rep;
all_info{1,k}=vis;
all_info{2,k}=time_SDP;
all_info{3,k}=B1;
all_info{4,k}=B2;
all_info{5,k}=T;
all_info{6,k}=ineq;
cd 2Bobs3Inputs
save(file_name,'all_info')
cd ..

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
sc=trace(S);
    gap=visc-vis;
    visc=vis;
end
