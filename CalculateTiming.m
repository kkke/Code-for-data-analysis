%% This code is used to extract the duration of central sampling and delay 
animalID = {'RVKC205','RVKC235','RVKC236','RVKC297','RVKC298','RVKC364'};
path = 'C:\Users\Ke Chen\Dropbox\TasteDiscriminationAnalysis_withKe\';

for j = 1:length(animalID)
    cd([path,animalID{j}])
    load('Summary.mat')
    for i = 1:length(summaryData)
        cen_L = summaryData(i).LickTs.Central;
        [~,cen_duration{i} cen_boutC{i}]= lickbout(cen_L);
    end
    
    s{j} = cell2mat(cen_duration);
    c{j} = cell2mat(cen_boutC);
end
%%
duration = cell2mat(s);
count    = cell2mat(c);
%% Calculate the delay
for j = 1:length(animalID)
    cd([path,animalID{j}])
    load('Summary.mat')
    for k = 1: length(summaryData)
        for i = 1:size(summaryData(k).raw_decision,1)
            delay(k,i) = summaryData(k).raw_decision{i,7}-summaryData(k).raw_decision{i,6}(end);
        end
        d{j} = delay(:)';
    end
end
delay = cell2mat(d);
delay(find(delay ==0))=[];
mean(delay)
std(delay)
%% Calculate the reaction time
for j = 1:length(animalID)
    cd([path,animalID{j}])
    load('Summary.mat')
    for k = 1: length(summaryData)
        for i = 1:size(summaryData(k).raw_decision,1)
            % calculate the left reaction time
            if  summaryData(k).raw_decision{i,2} ==1 & summaryData(k).raw_decision{i,3}==1
                reaction(k,i) =   summaryData(k).raw_decision{i,8}(1)-summaryData(k).raw_decision{i,7};
            elseif summaryData(k).raw_decision{i,2} ==2 & summaryData(k).raw_decision{i,3}==0
                reaction(k,i) =   summaryData(k).raw_decision{i,8}(1)-summaryData(k).raw_decision{i,7};
            elseif summaryData(k).raw_decision{i,2} ==1 & summaryData(k).raw_decision{i,3}==0
                reaction(k,i) =   summaryData(k).raw_decision{i,9}(1)-summaryData(k).raw_decision{i,7};
            else  summaryData(k).raw_decision{i,2} ==2 & summaryData(k).raw_decision{i,3}==1
                reaction(k,i) =   summaryData(k).raw_decision{i,9}(1)-summaryData(k).raw_decision{i,7};
            end
        end
        r{j,k} = reaction(:)';
        clear reaction
        
    end
end

delay = cell2mat(d);
delay(find(delay ==0))=[];
mean(delay)
std(delay)

