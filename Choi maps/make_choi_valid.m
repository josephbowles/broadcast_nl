function [Cout,distance_from_Cin_to_Cout,eigenvalues_of_Cout,PartialTrace_on_output_space] = make_choi_valid(T, din, dout)

% The input should be the Choi operator which is "close" to a quantum channel.
%That is, a bipartite operator T \in H_{in}\otimes H_{out} where T is close
%to be positive, and PartialTrace(T,2,[d_in d_out]) is close to identity
%din and dout are the input space and output space dimension

%The output of this function is a Choi operator which satisfies the Channel
%conditions "matlab exactly". (Where matlab means that they are constructed
%in a way that the constraints are exact, but due to matlab float variable
%approximations, it will not be exact


Cin = T;

Cout=(Cin+Cin')/2; %Ensures that Cout is self adjoint

Cout = Cout - kron(PartialTrace(Cout,2,[din dout]),eye(dout)/dout) + trace(Cout)*eye(din*dout)/(din*dout);
%Ensures that the Channel causal order TP condition is satisfied

%Uncoment this line if you want to impose that the Choi operator
%corresponds to a unital map;
    %Cout = Cout - kron(eye(din)/din, PartialTrace(Cout,1,[din dout])) + trace(Cout)*eye(din*dout)/(din*dout);
%Ensures that the Channel causal order Unital condition is satisfied

Cout=(Cout+Cout')/2; %Ensures that Cout is self adjoint
%(In principle we don't need that, but we add it to reduce matlabs
%imprecisions)

min_eig=min(real(eig(Cout)));
     if min_eig<=0
        Cout= Cout - min_eig*eye(din*dout) + 10^(-9)*eye(din*dout);
    end
%Ensures that Cout is positive definite;

Cout=Cout/trace(Cout)*din;
%Impose that Cout has the correct trace


 distance_from_Cin_to_Cout=norm(Cin-Cout);
 
 
if distance_from_Cin_to_Cout>0.01
    disp('Warning, the output is very far from the input!!')
%     "Your input Choi operator is very far from being a quantum channel"
%     "The Norm distance is:"
%     distance_from_Cin_to_Cout=distance_from_Cin_to_Cout;
end


% Final double check:
eigenvalues_of_Cout=eig(Cout);
PartialTrace_on_output_space = PartialTrace(Cout,2,[din dout]);


end