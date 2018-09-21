% Imaging analysis for GC project
%% rigister your image
% segment your image with Suite2p
%% load imaging data and ROI
clear
load('F_RVKC378_180904_plane1_proc')
idx = [dat.stat.iscell];
Fcell = dat.Fcell{1,1}(find(idx==1),:);  % extract the neuron;
FcellNeu = dat.FcellNeu{1,1}(find(idx==1),:);
% clearvars -except Y_r cc Cn
% clear
% filename = 'reg2.tif'
% data=MulRoiANA(filename)
% load C_df data from CNMF result
% %%
% save('data.mat','data')
%% load the event
[analog,trial] = process_intan_v2('RVKC378_180904.rhd');

%% load the licking
% analog2 = SettingTwoPhoton('RVKC368_180731');
% for i = 1:length(analog2)
%     info(i) = analog2p_gc(analog2(i).data,1000);
% end
% %% check the signal
% figure
% for i = 1:3
%     subplot(1,3,i)
% %     plot(info(i).time, analog2(i).data(:,3));
%     plot(info(i).time, -analog2(i).data(:,2));
%     hold on
%     scatter(info(i).lick, 0.2*ones(size(info(i).lick)));
% end
%% generate a trial structure; I know each trial contains 60 frames; and licking data are aligned to the tone;
dF = reshape(full(Fcell), size(Fcell,1),60,[]);
for i = 1:length(trial)
    trial(i).trace = dF(:,:,i);
end
save('data.mat','Fcell','trial','dat','FcellNeu')
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
%% Plot Tone response
for j =1:6
n = j;
figure;
for i = 1: length(trial)
    plot(trial(i).framT, trial(i).traceSmooth_dF(n,:),'Color',[0,0,0,0.5])
    hold on
end
title(['Neuron # ', num2str(n)])
xlim([-1.5,7])

end

%% align to the tastant
for i = 1:length(trial)
    if ~isnan(trial(i).S)
        trial(i).Time_Taste = trial(i).framT-trial(i).S(1);
    elseif ~isnan(trial(i).N)
        trial(i).Time_Taste = trial(i).framT-trial(i).N(1);
    elseif ~isnan(trial(i).CA)
        trial(i).Time_Taste = trial(i).framT-trial(i).CA(1);
    elseif ~isnan(trial(i).Q)
        trial(i).Time_Taste = trial(i).framT-trial(i).Q(1);
    elseif ~isnan(trial(i).W)
        trial(i).Time_Taste = trial(i).framT-trial(i).W(1);
    end
end
for i = 1:length(trial)
    idx = find(trial(i).Time_Taste>-4 & trial(i).Time_Taste<4);  % should be 50
    if length(idx)==50
        trial(i).T = trial(i).Time_Taste(idx);
        trial(i).Taste = trial(i).traceSmooth_dF(:,idx);
    elseif length(idx) ==49 && idx(end)+1 <=60
        idx(end+1) = idx(end)+1;
        trial(i).T = trial(i).Time_Taste(idx);
        trial(i).Taste = trial(i).traceSmooth_dF(:,idx);
    else
        trial(i).T = trial(i).Time_Taste(idx);
        trial(i).Taste = trial(i).traceSmooth_dF(:,idx);
    end 
end
%%
figure;
for i = 1:length(trial)
    if length(trial(i).T) ==50
        T(i,:) = trial(i).T;
    else
        a = length(trial(i).T);
        for j = 1: (50-a)
            trial(i).T(a+j) = trial(i).T(a) + j * 0.161; % in case that there are less than  50 bin, I add more bins and pad the trace with 0;
            trial(i).Taste(:,a+j) = zeros(size(trial(1).Taste,1),1);
        end
    end
end
T = mean(T,1);
for i = 1:length(trial)
    trial(i).Tpro = T; % creat a proximate time for all tastant, as tastant may jitter a little bit.
end

%% reorganize the data
neuron = trial2neuron5tastant(trial);
%% stats for each tastant
%% statistical test
t = 0.05;
rw = 3;
for j = 1:length(neuron)
    idx = find(neuron(j).T>-1 & neuron(j).T <0);
    T_idx1 = find(neuron(j).T>0 & neuron(j).T <rw);
    
   % S respose
    S_baseline    = mean(neuron(j).S_Taste_dF(:,idx),2);
    S_Taste_1     = mean(neuron(j).S_Taste_dF(:,T_idx1),2);   
    [p(1),h(1)] = ranksum(S_baseline,S_Taste_1);
    if mean(S_Taste_1)< mean(S_baseline) || p(1)>t || mean(S_Taste_1)<0
        h(1) = 0;
    end

    M_baseline    = mean(neuron(j).M_Taste_dF(:,idx),2); % 2nd taste
    M_Taste_1     = mean(neuron(j).M_Taste_dF(:,T_idx1),2);
    [p(2),h(2)] = ranksum(M_baseline,M_Taste_1);
    if mean(M_Taste_1)< mean(M_baseline)|| p(2)>t || mean(M_Taste_1)<0 
        h(2) = 0;
    end
    
    CA_baseline    = mean(neuron(j).CA_Taste_dF(:,idx),2); % 3rd taste
    CA_Taste_1     = mean(neuron(j).CA_Taste_dF(:,T_idx1),2);
    [p(3),h(3)] = ranksum(CA_baseline,CA_Taste_1);
    if mean(CA_Taste_1)< mean(CA_baseline)|| p(3)>t
        h(3) = 0;
    end
    
    Q_baseline    = mean(neuron(j).Q_Taste_dF(:,idx),2); % 4th taste
    Q_Taste_1     = mean(neuron(j).Q_Taste_dF(:,T_idx1),2);
    [p(4),h(4)] = ranksum(Q_baseline,Q_Taste_1);
    if mean(Q_Taste_1)< mean(Q_baseline)|| p(4)>t
        h(4) = 0;
    end   
    
    W_baseline    = mean(neuron(j).W_Taste_dF(:,idx),2); % 4th taste
    W_Taste_1     = mean(neuron(j).W_Taste_dF(:,T_idx1),2);
    [p(5),h(5)] = ranksum(W_baseline,W_Taste_1);
    if mean(W_Taste_1)< mean(W_baseline)|| p(5)>t || mean(W_Taste_1)<0
        h(5) = 0;
    end    
    resp(j).Sres = h(1);
    resp(j).Mres = h(2);
    resp(j).CAres = h(3);
    resp(j).Qres = h(4);
    resp(j).Wres = h(5);
end
%%
% resp1 = tasteResponse(neuron);
%%
plot_dF(8,neuron)

%% statistical test here the baseline is the 1 s before the cue; test cue response
t = 0.01;
rw = 2;
for j = 1:length(neuron)
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
    clear baseline C_1
end
%% check the cue respone
plot_dF_cue(15,neuron)

%%
save('data.mat','trial','Fcell','neuron','resp','dat','FcellNeu')
%% trying to extract the licking activity
% load('summary.mat') % summary.mat is saved be beha_2p
% plotLickIntan(summarydata)
% realign to the first lick after tone
for i = 1:length(trial)
    idx = find(trial(i).licks >1.9);
    trial(i).time_lick = trial(i).framT-trial(i).licks(idx(1));
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
for i = 1:length(trial)
    if length(trial(i).TLick)==50
        TLick(i,:) = trial(i).TLick;
    else
        error('Something is wrong')
    end
end

for i = 1:length(neuron)
    neuron(i).TLick = mean(TLick,1);
    for j = 1:length(trial)
    neuron(i).Licktrace_dF(j,:) = trial(j).LickTrace(i,:);
    end    
end
%% statistical test here the baseline is the 1 s before the cue; test lick response
t = 0.05;
rw = 1;
for j = 1:length(neuron)
    idx = find(trial(1).framT>-1 & trial(1).framT <0);
    C_idx1 = find(trial(1).TLick>0 & trial(1).TLick <rw);
    for i = 1:length(trial)
        baseline(i)   = mean(trial(i).traceSmooth_dF(j,idx),2);
        C_1(i)        = mean(trial(i).LickTrace(j,C_idx1),2);
    end
    [c_p(1),c_h(1)] = ranksum(baseline,C_1);
    if mean(C_1)< mean(baseline) || c_p(1)>t || mean(C_1)<0 ;
        c_h(1) = 0;
    end

    resp(j).LickRes = c_h(1);

end
%%

%%
save('data.mat','trial','Fcell','neuron','resp','dat')
%%
resp3 = tasteResponse3(neuron,trial);
resp2 = tasteResponse2(neuron,trial);
save('data.mat','trial','Fcell','neuron','resp','dat','resp3')

%% trying to find the location of the active neuron
% plot the spatial map
ind = [dat.stat.iscell];
loc = dat.stat(find(ind ==1));
image = dat.ops.mimg1;
figure;imshow(image,[100,3000])
hold on
BW =zeros(512,512);
for i = 1:length(loc) 
    for j = 1:length(loc(i).xpix)
        BW(loc(i).ypix(j), loc(i).xpix(j))=1;
    end
    [B{i},L{i}] = bwboundaries(BW,'noholes');
    BW = zeros(512,512);
 end
xlim([1,512]);ylim([1,512])
for k = 1:length(B)
   boundary = B{k};
   plot(boundary{1}(:,2), boundary{1}(:,1), 'r', 'LineWidth', 1)
   text(mean(boundary{1}(:,2)),mean(boundary{1}(:,1)),num2str(k),'Color',[0,0,0],'FontSize',16)
end
%% find the proportion of responsive neuron
cue_resp = [resp3.CueRes];
idx_cue  = find(cue_resp==1);
lick_resp = [resp3.LickRes];
idx_lick = find(lick_resp==1);
idx_lick = setdiff(idx_lick,idx_cue);
%% export as pdf
addpath('E:\MATLAB\Imaging Analysis\ExportPdf')
for i = 1:length(idx_cue)
    plot_dF_cue(idx_cue(i),neuron)
    export_fig(sprintf('cueRes%d', i), '-pdf');
    list{i} = ['cueRes',num2str(i),'.pdf']
end
append_pdfs('summary.pdf',list{:})
%% summarize the data across animals: Session1
summaryImaging_v1



%%

% prop_taste = length(find(sum(resp,1)>0))/size(resp,2);