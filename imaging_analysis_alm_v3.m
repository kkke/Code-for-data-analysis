% Imaging analysis for GC project
%% rigister your image
% segment your image with Suite2p
edit segmentation_Suite2P
%%
regtiff
segmentation_CNMF
%% load imaging data and ROI
load('data_CNMF.mat')
clearvars -except F0 F_dff cc Cn jsf
% clear
Y_r = F0.*F_dff+F0;
%% load the event
[analog,trial] = process_intan_v1('RVKC366_180808.rhd');
%%
%% generate a trial structure; I know each trial contains 60 frames; and licking data are aligned to the tone;
dF = reshape(full(Y_r), size(Y_r,1),60,[]);
for i = 1:length(trial)
%     trial(i).licks = info(i).lick - info(i).tone;
    trial(i).trace = dF(:,:,i);
end
save('data_CNMF_taste.mat','Y_r','trial','cc','Cn','jsf')
%% get the timestamps of each frame after 5 time averaging
for i = 1:length(trial)
    idx = 3:5:300;
    trial(i).framT =trial(i).Frame(idx);
end
%% Smooth the time-series with gaussian kerner
for i = 1:length(trial)
    for j = 1:size(trial(1).trace,1)
        trial(i).traceSmooth(j,:) = gaussmooth(trial(i).trace(j,:),5,1);
    end  
end
%% calculate the dF/F0; F0 is 1 s before the tone
for i = 1:length(trial)
    idx = find(trial(1).framT>-1 & trial(1).framT<0);
    baseline =mean(trial(i).traceSmooth(:,idx),2);
    trial(i).traceSmooth_dF = (trial(i).traceSmooth-repmat(baseline,1,size(trial(i).traceSmooth,2)))./repmat(baseline,1,size(trial(i).traceSmooth,2));
end
%%
t = 0.01;
rw = 1;
for j = 1:size(F0,1)
    idx = find(trial(1).framT>-1 & trial(1).framT <0);
    C_idx1 = find(trial(1).framT>0 & trial(1).framT <rw);
    for i = 1:length(trial)
        baseline(i)   = mean(trial(i).traceSmooth_dF(j,idx),2);
        C_1(i)        = mean(trial(i).traceSmooth_dF(j,C_idx1),2);
    end
    [c_p(1),c_h(1)] = ranksum(baseline,C_1);
    if mean(C_1)< mean(baseline) || c_p(1)>t || mean(C_1)<0 ;
        c_h(1) = 0;
    end
    resp(j).CueRes = c_h(1);
end
%% trying to extract the licking activity
% load('summary.mat') % summary.mat is saved be beha_2p
% plotLickIntan(summarydata)
% realign to the first lick after tone
for i = 1:length(trial)
    idx = find(trial(i).licks >1);
    if ~isempty(idx)
        trial(i).time_lick = trial(i).framT-trial(i).licks(idx(1));
    else
        trial(i).time_lick = [];   
    end
end

for i = 1:length(trial)
    idx = find(trial(i).time_lick>-4 & trial(i).time_lick<4);  % should be 50
    if length(idx)==50
        trial(i).TLick = trial(i).time_lick(idx);
        trial(i).LickTrace = trial(i).traceSmooth_dF(:,idx);
    elseif length(idx) ==49
        idx(end+1) = idx(end)+1;
        trial(i).TLick = trial(i).time_lick(idx);
        trial(i).LickTrace = trial(i).traceSmooth_dF(:,idx);       
    end 
end
%%
neuron = trial2neuron(trial);
%% statistical test here the baseline is the 1 s before the cue; test lick response
t = 0.01;
rw = 2;
clear C_1 baseline
for j = 1:size(F0,1)
    idx = find(trial(1).framT>-1 & trial(1).framT <0);
    C_idx1 = find(trial(1).TLick>0 & trial(1).TLick <rw);
    for i = 1:length(trial)
        baseline(i)   = mean(trial(i).traceSmooth_dF(j,idx),2);
        C_1(i)        = mean(trial(i).LickTrace(j,C_idx1),2);
    end
    [c_p(1),c_h(1)] = ranksum(baseline,C_1);
    avg =  mean(neuron(j).Licktrace_dF,1);
    if mean(C_1)< mean(baseline) || c_p(1)>t || mean(C_1)<0 || max(avg(C_idx1))<0.05
        c_h(1) = 0;
        
    end
    resp(j).LickRes = c_h(1);

end
%%
save('data_CNMF_taste.mat','Y_r','trial','cc','Cn','jsf','neuron','resp')

%%
plot_dF_lick(5,neuron)
plot_dF_cue(15,neuron)
%%
resp2 = pdResponse(neuron,trial,2);