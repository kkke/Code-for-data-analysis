clear
addpath('C:\Users\ke-roberto\Documents\MATLAB\Code_Ke\')
%% load raw data and extract the stimulation signals
thr = 2;
file = dir('*.rhd');
dataRaw = read_Intan(file.name);
[data.OptSt_on,data.OptSt_off] = Timing_onset_offset(dataRaw.analog(4,:), dataRaw.ts, thr,30,1);

%% load the summary file for the animal
cd('D:\Data_backUp\SummaryMiceDiscrimination\')
load('Summary.mat')
%% Let's make sure the every analysis is based on the most recent data
trial = summaryData(end);
%% For planning phase
for i = 1: size(trial.raw_decision,1)
    if ~isempty(trial.raw_decision{i,6})
        central_last(i) = trial.raw_decision{i,6}(end) + trial.raw_decision{i,5};
        central_first(i) = trial.raw_decision{i,6}(1) + trial.raw_decision{i,5};
        Left =  trial.raw_decision{i,8};
        Right = trial.raw_decision{i,9};
        lateral(i) = min([Left,Right])+ trial.raw_decision{i,5};
    end
end
sampling_duration = central_last - central_first;
delay_duration    = lateral-central_last;
k=1;
for i = 1: size(trial.raw_decision,1)
    for j = 1:length(data.OptSt_on)
        if data.OptSt_on(j)> central_last(i) && data.OptSt_on(j)< central_last(i)+2
            sti(k) = i;
            k = k+1;
        end
    end
end
%% Get the performance of the trials on stimulaiton
performance_sti(:,1) = cell2mat(trial.raw_decision(sti,3));
performance_sti(:,2) = cell2mat(trial.raw_decision(sti,2));
L_R = length(find(performance_sti(:,2)==1))/(length(find(performance_sti(:,2)==1))+length(find(performance_sti(:,2)==2)));
avg_per_sti = sum(performance_sti(:,1))/size(performance_sti,1);
%% get the performance of the trials without stimulation
perf = cell2mat(trial.raw_decision(:,3));
perf(sti)=[];
avg_per_non = sum(perf)/length(perf);
figure;
h1 = bar([avg_per_non,avg_per_sti],'FaceColor','w','EdgeColor','k','LineWidth',1);
ylim([0,1])
set(gca,'xticklabel',{'Non-Sti','Sti'})
set(h1,'BarWidth',0.5)

%%
clearvars -except delay_duration sampling_duration avg_per_non avg_per_sti
