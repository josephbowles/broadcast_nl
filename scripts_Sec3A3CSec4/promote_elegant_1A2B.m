%% Add relevant directories to the path

mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'aux_funcSec3A3CSec4',filesep);
newdir3 = strcat(newdir,filesep,'QETLAB-0.9',filesep);
addpath(newdir2);
addpath(genpath(newdir3));

%%

C1 = sym('C1');
C2 = sym('C2');
B1 = sym('B1');
B2 = sym('B2');
B3 = sym('B3');
B4 = sym('B4');
A1 = sym('A1');
A2 = sym('A2');
A3 = sym('A3');
A4 = sym('A4');


ins = [3,4];
outs = [2,2];
% expression taken from eq 24 in https://arxiv.org/pdf/1607.08182.pdf
EBICorrelator = A1*B1 + A1*B2 - A1*B3 - A1*B4 + A2*B1 - A2*B2 + ...
                A2*B3 - A2*B4 + A3*B1 - A3*B2 - A3*B3 + A3*B4;
FLAG_Use01obsInsteadOfCorrelator = false;
EBIProjector = ToProjectorNotation(EBICorrelator,ins,outs,FLAG_Use01obsInsteadOfCorrelator);
EBIProbability = FromProjectorNotationToProbability(EBIProjector, ins, outs);
EBIBellcoeffs = GetBellCoeffsFromProbSymIneq(EBIProbability,ins,outs);
EBILocalBound =  ClassicalOptInequalityNoBroadcast(EBIBellcoeffs) % should be 4 for 3322 chained

ins = [4,4,2];
outs = [2,2,2];
ModifiedEBI = EBICorrelator * (C1 + C2) + EBILocalBound * A4 *(C1 - C2);
ModifiedEBI = expand(ModifiedEBI);

ModifiedEBIProjector = ToProjectorNotation(ModifiedEBI,ins,outs,FLAG_Use01obsInsteadOfCorrelator);
ModifiedEBIProbability = FromProjectorNotationToProbability(ModifiedEBIProjector, ins, outs);
ModifiedEBIBellcoeffs = GetBellCoeffsFromProbSymIneq(ModifiedEBIProbability,ins,outs);
ModifiedEBILocalBound = ClassicalOptInequalityNoBroadcast(ModifiedEBIBellcoeffs)

bellcoeffs = ModifiedEBIBellcoeffs;

%% !! bellcoeffs called as bellcoeffs(x,y,z,a,b,c) and NOT bellcoeffs(a,b,c,x,y,z)
save('bellcoeffsEBIbroadcast.mat','bellcoeffs','ins','outs');

