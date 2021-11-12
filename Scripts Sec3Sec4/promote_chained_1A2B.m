

%% Add relevant directories to the path

mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'Aux functions Sec3Sec4',filesep);
newdir3 = strcat(newdir,filesep,'QETLAB-0.9',filesep);
addpath(newdir2);
addpath(genpath(newdir3));

%%

C1 = sym('C1');
C2 = sym('C2');
B1 = sym('B1');
B2 = sym('B2');
B3 = sym('B3');
A1 = sym('A1');
A2 = sym('A2');
A3 = sym('A3');
A4 = sym('A4');

ins = [3,3];
outs = [2,2];
% expression taken from eq 24 in https://arxiv.org/pdf/1607.08182.pdf
ChainedCorrelator = A1*B1 - A1*B3 + A2*B1 + A2*B2 + A3*B2 + A3*B3;
FLAG_Use01obsInsteadOfCorrelator = false;
ChainedProjector = ToProjectorNotation(ChainedCorrelator,ins,outs,FLAG_Use01obsInsteadOfCorrelator);
ChainedProbability = FromProjectorNotationToProbability(ChainedProjector, ins, outs);
ChainedBellcoeffs = GetBellCoeffsFromProbSymIneq(ChainedProbability,ins,outs);
ChainedLocalBound = ClassicalOptInequalityNoBroadcast(ChainedBellcoeffs) % should be 4 for 3322 chained

ins = [4,3,2];
outs = [2,2,2];
ModifiedChained = ChainedCorrelator * (C1 + C2) + ChainedLocalBound * A4 *(C1 - C2);
ModifiedChained = expand(ModifiedChained);

ModifiedChainedProjector = ToProjectorNotation(ModifiedChained,ins,outs,FLAG_Use01obsInsteadOfCorrelator);
ModifiedChainedProbability = FromProjectorNotationToProbability(ModifiedChainedProjector, ins, outs);
ModifiedChainedBellcoeffs = GetBellCoeffsFromProbSymIneq(ModifiedChainedProbability,ins,outs);
ModifiedChainedLocalBound = ClassicalOptInequalityNoBroadcast(ModifiedChainedBellcoeffs)

bellcoeffs = ModifiedChainedBellcoeffs;

%% !! bellcoeffs called as b(x,y,z,a,b,c) and NOT b(a,b,c,x,y,z)
save('bellcoeffsChained3322broadcast.mat','bellcoeffs','ins','outs');

