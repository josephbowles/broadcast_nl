%Heuristic search in the broadcast scenario

rng('shuffle'); %Ensures that the Heuristic starts with different random seeds

%This code tries to find the smallest visibility 'v' such that v*rho + (1-v)*rhoNoise is
%broadcast nonlocal, where 'rho' is typically nonlocal while 'rhoNoise' is
%separable (see https://arxiv.org/abs/2007.16034)

%The variables can be set using script "PrepareAndRun_HeuristicSearchBroadcast.m 


%Variables needed :

%'rho' is an arbitrary multipartite state, typically entangled/nonlocal
%'rhoNoise' is an arbitrary multipartite state of the same dimension as 'rho', 'rhoN' is typically separable and represents the noisy part of the state
%'WChois' sets the (initial) parties where one wants to apply a Choi map
%'dloc' sets the local dimensions of the final parties
%'nInputs' sets the number of inputs of each final party
%'nOutputs' sets the number of outputs of each final party
%'VerticesAllbutLast' sets the strategies of all parties, but the last (where no-signaling is already built-in)



%Less important degrees of freedom (already set in the script):
%'prec' sets the gap between two succesive visibilites at which the search stops (default: 10^-5)
%'limcount' sets the maximum number of tries to find a violation of the current best inequality (default: 5)

%One might also enter manually some of the initial parameters of the search
%(variables 'Meas' and 'Choi' in the first loop), or can set which Choi states/ which measurements to optimize the
%inequalities on (variables 'WhichChoi' and 'WhichMeas')





Nparties = length(dims); %number of initial parties



%Output dimensions of Choi states

dOut=zeros(1,Nparties);
for k=1:length(dloc)
    
    dOut(k) = prod(dloc{k});
    
end


%Polytope last party 
VerticesB=local_polytope(length(nInputs{Nparties}),nInputs{Nparties},nOutputs{Nparties}); %change only if other strategies than no-signaling polytope are required for last party
%local vertices are generated here, but one then take an arbitrary linear combination of them -instead of convex.
%This is equivalent to considering the set of all no-signaling strategies ((see e.g. Theorem 1 of https://arxiv.org/pdf/0804.4859.pdf)


%numbers of overall outputs for last party
outB=prod(nOutputs{Nparties});



%flat scenario
dlocvec = cell2mat(dloc);
nInputsflat = cell2mat(nInputs);
nOutputsflat = cell2mat(nOutputs);




%Search a point outside

%stores succesive values of visibility in a vector
Visibility_vector=0;

vis = 2;
while vis > 1
    
    
    %Random (projective) Choi states
    Choi=cell(1,length(dims));
    for k=1:length(dims)
        if ismember(k,WChoi)
            Choi{k}=RandomChoi_proj(dims(k),dOut(k));
        else
            Choi{k}=ChoiId(dims(k));
        end
    end
    
    
    %Random (projective) measurements
    Meas=cell(1,length(nInputsflat));
    for k=1:length(nInputsflat)
        Meas{k}= rand_povms(dlocvec(k),nOutputsflat(k),nInputsflat(k),1);
    end
    
    
    %Enter here optional fixed initial Choi states/measurements
    %     Meas{1} = blah ;
    %
    %     Choi{2} = blouh ;
    
    
    
    %Final states
    
    rhoFinal = ApplyManyChois(rho,Choi,dims);
    
    rhoNoiseFinal = ApplyManyChois(rhoNoise,Choi,dims);
    
    
    %Quantum behaviours
    
    BroadcastBehaviourRhoFinal = MultipartiteBehaviour(rhoFinal,Meas);
    
    BroadcastBehaviourRhoNoiseFinal = MultipartiteBehaviour(rhoNoiseFinal,Meas);
    
    %add some uniform noise, for stability (not clear it's useful with
    %InteriorPoint algorithm)
     p=10^-5;
     BroadcastBehaviourRhoNoiseFinal = (1-p)*BroadcastBehaviourRhoNoiseFinal + p*MaximallyMixedBehaviour(Nparties,nInputsflat,nOutputsflat);
%     
    
    %Visibility
    
    [vis,~,~,~,lambda] = p_maxPConvexAffineabxy(BroadcastBehaviourRhoFinal',BroadcastBehaviourRhoNoiseFinal',VerticesAllbutLast,VerticesB,[outA,outB]);
    
    vis = vis(1)
    
%     %If bug in the linprog stops the loop
%     if length(vis) > 1
%         disp 'Error in linprog, flag is'
%         disp(vis(2))
%         return
%     end
    
    
end


%Extract inequality
Inequality=-lambda.eqlin;

%Compute the local bound of the inequality
violc = Inequality'*(vis*BroadcastBehaviourRhoFinal'+(1-vis)*BroadcastBehaviourRhoNoiseFinal');


%Store the results
all_info_HeuristicIsotropicState = cell(1,5);

all_info_HeuristicIsotropicState{1,1}=vis;  all_info_HeuristicIsotropicState{1,2}=Choi;  all_info_HeuristicIsotropicState{1,3}=Meas;

all_info_HeuristicIsotropicState{1,4}=Inequality;  all_info_HeuristicIsotropicState{1,5}=violc;


%Start the heuristic optimization (optimizes inequality, then computes
%visibility of resulting point and new inequality, and so on)


CurrentVis=vis;

%Store the succesives visibilites into a vector
Visibility_vector(1)=vis(1);

%Stop when two succesives visibilities are closer than 'prec'
prec=10^-5;


%Start the loop
gap=1; iteration = 2;
while gap > prec
    
    
    %Current best state
    rhoc=CurrentVis*rho+(1-CurrentVis)*rhoNoise;
    
    
    %Seesaw of the current inequality over measurements and Chois
    
    WhichParties = 1:length(dlocvec); %Which parties to optimize over
    WhichChois = WChoi; %Which Choi states to optimize over
    
    %Tries 'limcount' times to violate the current inequality, stops when succeeds
    limcount=5;
    
    viol=violc-1;
    count=0;
    while viol < violc+10^-5
        
        %Tries using previous settings as starting point
        if count == 2
            
            f=SeeSawBroadcastIneq(Inequality,rhoc,dims,WhichParties,WhichChois,Choi,Meas);
            g=f;
            
            
        else
            
            %Random initial points:
            
            
            %Random (projective) Choi states
            ChoiR=cell(1,length(dims));
            for k=1:length(dims)
                if ismember(k,WChoi)
                    ChoiR{k}=RandomChoi_proj(dims(k),dOut(k));
                else
                    ChoiR{k}=ChoiId(dims(k));
                end
            end
            
            
            %Random (projective) measurements
            MeasR=cell(1,length(nInputsflat));
            for k=1:length(nInputsflat)
                MeasR{k}= rand_povms(dlocvec(k),nOutputsflat(k),nInputsflat(k),1);
            end
            
            
            f=SeeSawBroadcastIneq(Inequality,rhoc,dims,WhichParties,WhichChois,ChoiR,MeasR);
            
        end
        
        yalmip clear
        
        
        %stores the current violation
        viol=f{1};
        
        
        count=count+1;
        
        %stops after 'limcount' attempts and keep the run using previous
        %settings
        if count > limcount && viol <= violc+10^-5
            
            f=g;
            break
            
        end
        
    end
    
    
    
    %Extract measurements and Choi states
    
    Meas = f{2};
    
    Choi = f{3};
    
    
    
    %Make all measurements and Chois valid
    
    for k=1:length(Meas)
        Meas{k}=make_povms_valid(Meas{k});
    end
    
    for k=1:length(Choi)
        Choi{k}=make_choi_valid(Choi{k},dims(k),dOut(k));
    end
    
    
    %Final states
    
    rhoFinal = ApplyManyChois(rho,Choi,dims);
    
    rhoNoiseFinal = ApplyManyChois(rhoNoise,Choi,dims);
    
    
    %Quantum behaviours
    
    BroadcastBehaviourRhoFinal = MultipartiteBehaviour(rhoFinal,Meas);
    
    BroadcastBehaviourRhoNoiseFinal = MultipartiteBehaviour(rhoNoiseFinal,Meas);
    
    
    %Visibility

    %add some uniform noise, for stability (not clear it's useful with
    %InteriorPoint algorithm)
     p=10^-5;
     BroadcastBehaviourRhoNoiseFinal = (1-p)*BroadcastBehaviourRhoNoiseFinal + p*MaximallyMixedBehaviour(Nparties,nInputsflat,nOutputsflat);
%     
    
    
    [vis,~,~,~,lambda] = p_maxPConvexAffineabxy(BroadcastBehaviourRhoFinal',BroadcastBehaviourRhoNoiseFinal',VerticesAllbutLast,VerticesB,[outA,outB]);
    
    vis = vis(1);
%     %If bug in the linprog stops the loop
%     if length(vis) > 1
%         disp 'Error in linprog, flag is'
%         disp(vis(2))
%         return
%     end
    
    
    
    
    Visibility_vector(iteration)=vis(1),
    
    
    %New inequality and bound
    
    Inequality=-lambda.eqlin;
    
    violc = Inequality'*(vis*BroadcastBehaviourRhoFinal'+(1-vis)*BroadcastBehaviourRhoNoiseFinal');
    
    
    %Store the results
    
    all_info_HeuristicIsotropicState{iteration,1}=vis;  all_info_HeuristicIsotropicState{iteration,2}=Choi;  all_info_HeuristicIsotropicState{iteration,3}=Meas;
    
    all_info_HeuristicIsotropicState{iteration,4}=Inequality;  all_info_HeuristicIsotropicState{iteration,5}=violc;
    
    
    
    gap=CurrentVis-vis;
    
    CurrentVis=vis;
    
    iteration = iteration + 1;
    
end








