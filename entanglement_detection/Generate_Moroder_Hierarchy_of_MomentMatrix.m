function [FinalMomentMatrix,Operator_gives_PhysicalMoment,DetailofEachMomentIndex,DetailofEachOperator] = Generate_Moroder_Hierarchy_of_MomentMatrix(LocalLevel,Num_Outputs_Inpust_EachPartie)
% Num_Outputs_Inpust_EachPartie =  [Oa,Ob,Oc;Sa,Sb,Sc]
% Order of PhysicalMoment =
% [P(a=1,b=1,c=1|x=1,y=1,z=1), P(a=2,b=1,c=1|x=1,y=1,z=1), ..., P(a=1,b=Oa-1,c=1|x=1,y=1,z=1),P(a=1,b=1,c=1|x=2,y=1,z=1), ... , P(a=Oa-1,b=1,c=1|x=Sa,y=1,z=1), P(a=1,b=1,c=1|x=Sa+1,y=1,z=1), P(a=1,b=1,c=1|x=1,y=2,z=1),
%  ... ,P(a=Oa-1,b=1,c=1|x=Sa,y=2,z=1), P(a=1,b=1,c=1|x=Sa+1,y=2,z=1),P(a=1,b=1,c=1|x=1,y=3,z=1), ... ]
Ceil_Level = ceil(LocalLevel);
Floor_Level = floor(LocalLevel);
Num_of_MarginalOperatorinCGForm = [Num_Outputs_Inpust_EachPartie(1,:)-1;Num_Outputs_Inpust_EachPartie(2,:)];
Num_Parties = length(Num_Outputs_Inpust_EachPartie);
Num_Operator_DifferentOrder = zeros(2*Ceil_Level,Num_Parties);
Num_of_BoxesToFilled = 2*Ceil_Level; % Each box represent choosing one element from measurements.
OperatorDetail_Index = cell(1,Num_Parties);


SaveDataFileName ='';
for i =1:2
    SaveDataFileName = strcat(SaveDataFileName,'_');
    for t = 1:Num_Parties
        SaveDataFileName = strcat(SaveDataFileName,num2str(Num_Outputs_Inpust_EachPartie(i,t)));
    end
end

SaveDataFileName = strcat('MomentmatrixData_LocalLevel_',num2str(LocalLevel),SaveDataFileName,'.mat');


try
    load(SaveDataFileName);
catch
    'Constructing momentmatrix'

% Generate possible operator for each parties
for k_partie = 1:Num_Parties
    
    TmpNumMarginalOpe = prod(Num_of_MarginalOperatorinCGForm(:,k_partie));
    Num_of_PossibilieCombinationsofMeasurements = TmpNumMarginalOpe^(Num_of_BoxesToFilled);
    TmpOpe = zeros(Num_of_PossibilieCombinationsofMeasurements,Num_of_BoxesToFilled+3);%+1 to record index of each possibile combinations of measurements, +1 to record the index after "flip" the measurements,
    %+1 to record the properties of projector.
    IndexOfBeingOrthogonal = Num_of_BoxesToFilled+3; IndexOfOperator_BeforeFlip = Num_of_BoxesToFilled+1; IndexOfOperator_AfterFlip = Num_of_BoxesToFilled+2;
    TmpDetail_MeasurementsFilled = cell(1,Num_of_BoxesToFilled+1); % Extral column is used to compare different elements.
    for  i = 1:Num_of_PossibilieCombinationsofMeasurements % All possibile combinations of measurements
        [TmpDetail_MeasurementsFilled{1:Num_of_BoxesToFilled}] = ind2sub(TmpNumMarginalOpe*ones(1,Num_of_BoxesToFilled),i);
        FLAG_of_Being_Orthogonal = 0;
        for j = 1:Num_of_BoxesToFilled
            if TmpDetail_MeasurementsFilled{j} == TmpDetail_MeasurementsFilled{j+1},
                TmpDetail_MeasurementsFilled{j} = 0;
            elseif ceil(TmpDetail_MeasurementsFilled{j}/(Num_Outputs_Inpust_EachPartie(1,k_partie)-1)) == ceil(TmpDetail_MeasurementsFilled{j+1}/(Num_Outputs_Inpust_EachPartie(1,k_partie)-1))
                FLAG_of_Being_Orthogonal = 1;
            end
        end
        TTT = cell2mat(TmpDetail_MeasurementsFilled);
        TmpOpe(i,Num_of_BoxesToFilled-length(find(TTT~=0))+1:Num_of_BoxesToFilled) = TTT(find(TTT~=0));
        TmpOpe(i,IndexOfBeingOrthogonal) = FLAG_of_Being_Orthogonal;
    end
    for i =2:Num_of_PossibilieCombinationsofMeasurements
        TmpOpe = sortrows(TmpOpe);
        if isequal(TmpOpe(i,:),TmpOpe(i-1,:))
            TmpOpe(i-1,:) = 0;
            TmpOpe = sortrows(TmpOpe);
        end
    end
    TmpOpe = TmpOpe(find(TmpOpe(:,Num_of_BoxesToFilled)~=0,1):Num_of_PossibilieCombinationsofMeasurements,:);
    TmpOpe(:,IndexOfOperator_BeforeFlip) = 1:size(TmpOpe,1);
    for i =1:size(TmpOpe,1) % Find the flip index
        TmpTarget = TmpOpe(i,1:Num_of_BoxesToFilled);
        TmpTarget(find(TmpTarget~=0)) = flip(TmpTarget(find(TmpTarget~=0)));
        TmpOpe(find(ismember(TmpOpe(:,1:Num_of_BoxesToFilled),repmat(TmpTarget,[size(TmpOpe,1),1]),'rows'),2,'last'),IndexOfOperator_AfterFlip) = i;
    end
    OperatorDetail_Index{k_partie} = TmpOpe; % All information of operators
    for L = Num_of_BoxesToFilled:-1:1
        if L ~=1
            Num_Operator_DifferentOrder(Num_of_BoxesToFilled-L+1,k_partie) = find(TmpOpe(:,L-1)~=0,1) - find(TmpOpe(:,L)~=0,1);
        else
            Num_Operator_DifferentOrder(Num_of_BoxesToFilled-L+1,k_partie) = size(TmpOpe,1) - find(TmpOpe(:,1)~=0,1)+1;
        end
    end
end

% Now create the moment matrix, by first creating matrix for operators then
% decide the correponding moment index.

if Ceil_Level == Floor_Level
    if Ceil_Level == 1
        TotalNumOperatorNeeded_EachPartie =  Num_Operator_DifferentOrder(1:Ceil_Level,:)+1;
    else
        TotalNumOperatorNeeded_EachPartie =  sum(Num_Operator_DifferentOrder(1:Ceil_Level,:))+1;
    end
elseif Floor_Level == 1
    TotalNumOperatorNeeded_EachPartie = Num_Operator_DifferentOrder(1:Floor_Level,:)+1;
    TotalNumOperatorNeeded_EachPartie = round(TotalNumOperatorNeeded_EachPartie + Num_Operator_DifferentOrder(Ceil_Level,:)*(LocalLevel-Floor_Level));
else
    TotalNumOperatorNeeded_EachPartie =  sum(Num_Operator_DifferentOrder(1:Floor_Level,:))+1;
    TotalNumOperatorNeeded_EachPartie = round(TotalNumOperatorNeeded_EachPartie + Num_Operator_DifferentOrder(Ceil_Level,:)*(LocalLevel-Floor_Level));
end
Operator_In_First_Row = zeros(prod(TotalNumOperatorNeeded_EachPartie),Num_Parties);
Operator_In_First_Column = zeros(prod(TotalNumOperatorNeeded_EachPartie),Num_Parties);

TmpDetailOperator_EachElement = cell(1,Num_Parties);
for i = 1:prod(TotalNumOperatorNeeded_EachPartie)
    [TmpDetailOperator_EachElement{:}] = ind2sub(TotalNumOperatorNeeded_EachPartie,i);
    Operator_In_First_Row(i,:) = cell2mat(TmpDetailOperator_EachElement);
end

for i =1:prod(TotalNumOperatorNeeded_EachPartie)
    for k_partie =1:Num_Parties
        if Operator_In_First_Row(i,k_partie) == 1
            Operator_In_First_Column(i,k_partie) = 1;
        elseif OperatorDetail_Index{k_partie}(Operator_In_First_Row(i,k_partie)-1,IndexOfBeingOrthogonal) ==1
            Operator_In_First_Row(i,k_partie) = 0;
            Operator_In_First_Column(i,k_partie) = 0;
        else
            Operator_In_First_Column(i,k_partie) = OperatorDetail_Index{k_partie}(Operator_In_First_Row(i,k_partie)-1,IndexOfOperator_AfterFlip)+1;
            
        end
        %             Operator_In_First_Column(i,k_partie) = OperatorDetail_Index{k_partie}(Operator_In_First_Row(i,k_partie)-1,IndexOfOperator_AfterFlip)+1;
        %             if OperatorDetail_Index{k_partie}(Operator_In_First_Row(i,k_partie)-1,IndexOfBeingOrthogonal) ==1
        %                 Operator_In_First_Row(i,k_partie) = 0;
        %             end
        %             if OperatorDetail_Index{k_partie}(Operator_In_First_Column(i,k_partie)-1,IndexOfBeingOrthogonal) ==1
        %                 Operator_In_First_Column(i,k_partie) = 0;
        %             end
    end
end


CombinationOperator_EachPartie = zeros(Num_Parties,2); % [IndexOperator_From_Row,IndexOperator_From_Col; IndexOperator_From_Row,IndexOperator_From_Col]
Reduced_Operator_EachPartie = zeros(prod(TotalNumOperatorNeeded_EachPartie)^2,Num_Parties); TmpRow_Operator = zeros(1,Num_of_BoxesToFilled); TmpCol_Operator = zeros(1,Num_of_BoxesToFilled);
FinalMomentMatrix = zeros(prod(TotalNumOperatorNeeded_EachPartie)); Tmp_MomentMatrix = zeros(prod(TotalNumOperatorNeeded_EachPartie)^2,Num_Parties+2);% Here k_partie columns + 1 for IndexOfMomentMatrix + 1 for SortedMomentMatrixIndex

for IndexofMomentMatrix = 1:prod(TotalNumOperatorNeeded_EachPartie)^2
    [Index_Row,Index_Col] = ind2sub(ones(1,k_partie)*prod(TotalNumOperatorNeeded_EachPartie),IndexofMomentMatrix);
    % First find the combination of operators for  each partie
    if IndexofMomentMatrix == 1
        FinalMomentMatrix(Index_Row,Index_Col) =1;
        Reduced_Operator_EachPartie(IndexofMomentMatrix,:) = 1;
    else
        for k_partie = 1:Num_Parties
            CombinationOperator_EachPartie(k_partie,:) = [Operator_In_First_Column(Index_Row,k_partie),Operator_In_First_Row(Index_Col,k_partie)];
            if find(CombinationOperator_EachPartie(k_partie,:)==0,1)
                Reduced_Operator_EachPartie(IndexofMomentMatrix,k_partie) = 0;
            elseif length(find(CombinationOperator_EachPartie(k_partie,:)==1)) == 2,
                Reduced_Operator_EachPartie(IndexofMomentMatrix,k_partie) = 1;
            elseif length(find(CombinationOperator_EachPartie(k_partie,:)==1)) & length(find(CombinationOperator_EachPartie(k_partie,:)~=1))
                Reduced_Operator_EachPartie(IndexofMomentMatrix,k_partie) = CombinationOperator_EachPartie(k_partie,find(CombinationOperator_EachPartie(k_partie,:)~=1));
            else
                TmpDetailOperator_RowCol = zeros(1,Num_of_BoxesToFilled);
                TmpRow_Operator = OperatorDetail_Index{k_partie}(CombinationOperator_EachPartie(k_partie,1)-1,1:Num_of_BoxesToFilled);
                TmpCol_Operator = OperatorDetail_Index{k_partie}(CombinationOperator_EachPartie(k_partie,2)-1,1:Num_of_BoxesToFilled);
                TmpDetailOperator_RowCol(Num_of_BoxesToFilled-length([TmpRow_Operator(find(TmpRow_Operator~=0)),TmpCol_Operator(find(TmpCol_Operator~=0))])+1:Num_of_BoxesToFilled)...
                    = [TmpRow_Operator(find(TmpRow_Operator~=0)),TmpCol_Operator(find(TmpCol_Operator~=0))];  % Add an additional column so that I can compare the whole elements
                TmpOperator = [TmpDetailOperator_RowCol,0];
                for t =1:Num_of_BoxesToFilled
                    if TmpOperator(t) == TmpOperator(t+1)
                        TmpOperator(t) = 0;
                    end
                end
                TmpDetailOperator_RowCol = zeros(1,Num_of_BoxesToFilled);
                TmpDetailOperator_RowCol(Num_of_BoxesToFilled-length(find(TmpOperator~=0))+1:Num_of_BoxesToFilled) = TmpOperator(find(TmpOperator~=0));
                TTT =  find(ismember(OperatorDetail_Index{k_partie}(:,1:Num_of_BoxesToFilled),repmat(TmpDetailOperator_RowCol,[size(OperatorDetail_Index{k_partie},1),1]),'rows'));
                if OperatorDetail_Index{k_partie}(TTT,IndexOfBeingOrthogonal)
                    Reduced_Operator_EachPartie(IndexofMomentMatrix,k_partie)= 0 ;
                else
                    Reduced_Operator_EachPartie(IndexofMomentMatrix,k_partie)=TTT+1;
                end
                
            end
        end
    end
    Tmp_MomentMatrix(IndexofMomentMatrix,1:Num_Parties) = Reduced_Operator_EachPartie(IndexofMomentMatrix,:);
end

% Sort moment matrix
Tmp_MomentMatrix(:,Num_Parties+1) = 1:prod(TotalNumOperatorNeeded_EachPartie)^2;
SortOrder = 1;
TmpTarget_Flip = zeros(1,Num_Parties);
for i =1:prod(TotalNumOperatorNeeded_EachPartie)^2
    if Tmp_MomentMatrix(i,Num_Parties+2)==0;
        TmpTarget_NoFlip = Tmp_MomentMatrix(i,1:Num_Parties);
        if find(TmpTarget_NoFlip ==0)
            TmpTarget_Flip = TmpTarget_NoFlip;
        else
            for k_partie =1:Num_Parties
                if Tmp_MomentMatrix(i,k_partie) == 1,
                    TmpTarget_Flip(k_partie) = 1;
                else
                    TmpTarget_Flip(k_partie) = OperatorDetail_Index{k_partie}(Tmp_MomentMatrix(i,k_partie)-1,IndexOfOperator_AfterFlip)+1;
                end
            end
        end
        if isequal(TmpTarget_NoFlip,TmpTarget_Flip)
            TTT = find(ismember(Tmp_MomentMatrix(:,1:Num_Parties),repmat(TmpTarget_NoFlip,[prod(TotalNumOperatorNeeded_EachPartie)^2,1]),'rows'));
            Tmp_MomentMatrix(TTT,Num_Parties+2) = SortOrder;
            SortOrder = SortOrder + 1;
        else
            TTT = find(ismember(Tmp_MomentMatrix(:,1:Num_Parties),repmat(TmpTarget_NoFlip,[prod(TotalNumOperatorNeeded_EachPartie)^2,1]),'rows'));
            Tmp_MomentMatrix(TTT,Num_Parties+2) = SortOrder;
            TTT = find(ismember(Tmp_MomentMatrix(:,1:Num_Parties),repmat(TmpTarget_Flip,[prod(TotalNumOperatorNeeded_EachPartie)^2,1]),'rows'));
            Tmp_MomentMatrix(TTT,Num_Parties+2) = SortOrder;
            SortOrder = SortOrder + 1;
        end
    end
end
for i =1:Num_Parties
    Tmp_MomentMatrix(find(Tmp_MomentMatrix(:,i)==0),Num_Parties+2)=0;
end
TTT = unique(Tmp_MomentMatrix(:,Num_Parties+2));
if TTT(1) == 0,
    for i =1:length(TTT)-1
        Tmp_MomentMatrix(find(Tmp_MomentMatrix(:,Num_Parties+2)==TTT(i+1)),Num_Parties+2)=i;
    end
else
    for i =1:length(TTT)
        Tmp_MomentMatrix(find(Tmp_MomentMatrix(:,Num_Parties+2)==TTT(i)),Num_Parties+2)=i;
    end
end
for IndexofMomentMatrix = 1:prod(TotalNumOperatorNeeded_EachPartie)^2
    [Index_Row,Index_Col] = ind2sub(ones(1,k_partie)*prod(TotalNumOperatorNeeded_EachPartie),IndexofMomentMatrix);
    FinalMomentMatrix(Index_Row,Index_Col) = Tmp_MomentMatrix(IndexofMomentMatrix,Num_Parties+2);
end

Operator_gives_PhysicalMoment = zeros(prod(Num_Operator_DifferentOrder(1,:)+1),Num_Parties+1);
TTT = cell(1,Num_Parties);
for i =1:prod(Num_Operator_DifferentOrder(1,:)+1)
    [TTT{:}] = ind2sub(Num_Operator_DifferentOrder(1,:)+1,i);
    Operator_gives_PhysicalMoment(i,1:Num_Parties) = cell2mat(TTT);
    Operator_gives_PhysicalMoment(i,Num_Parties+1)= unique(Tmp_MomentMatrix(find(ismember(Tmp_MomentMatrix(:,1:Num_Parties),repmat(Operator_gives_PhysicalMoment(i,1:Num_Parties),[size(Tmp_MomentMatrix,1),1]),'rows')),Num_Parties+2));
end
Operator_gives_PhysicalMoment(:,1:Num_Parties) = Operator_gives_PhysicalMoment(:,1:Num_Parties)-1;
for i =1:Num_Parties,
    Tmp = Operator_gives_PhysicalMoment(:,i); 
    Tmp(find(Tmp==0)) = max(Tmp)+1;
    Operator_gives_PhysicalMoment(:,i) = Tmp;
end
TmpOperator_gives_PhysicalMoment = zeros(size(Operator_gives_PhysicalMoment,1),2*Num_Parties+1);
T = cell(1,2); 
Num_Outputs_Inpust_EachPartie(1,:) = Num_Outputs_Inpust_EachPartie(1,:)-1;
for i =1:size(Operator_gives_PhysicalMoment,1)
    for k = 1:Num_Parties
        if Operator_gives_PhysicalMoment(i,k) <= prod(Num_Outputs_Inpust_EachPartie(:,k))
            [T{:}] = ind2sub(Num_Outputs_Inpust_EachPartie(:,k)',Operator_gives_PhysicalMoment(i,k) );
        else
            [T{:}] = ind2sub([1,Num_Outputs_Inpust_EachPartie(2,k)+1],prod([1,Num_Outputs_Inpust_EachPartie(2,k)+1]));
        end
        TmpOperator_gives_PhysicalMoment(i,1+(k-1)*2:k*2) = [T{:}];
    end    
end
TmpOperator_gives_PhysicalMoment(:,end) = Operator_gives_PhysicalMoment(:,end);
TmpOperator_gives_PhysicalMoment = sortrows(TmpOperator_gives_PhysicalMoment,[2*(Num_Parties:-1:1),2*(Num_Parties:-1:1)-1]);
Operator_gives_PhysicalMoment = TmpOperator_gives_PhysicalMoment; 
DetailofEachOperator = struct('OperatorDetail_Index',{OperatorDetail_Index},...
    'Num_Operator_DifferentOrder',Num_Operator_DifferentOrder);
DetailofEachMomentIndex = Tmp_MomentMatrix;

save(SaveDataFileName,'FinalMomentMatrix','Operator_gives_PhysicalMoment','DetailofEachMomentIndex','DetailofEachOperator');
end
end

