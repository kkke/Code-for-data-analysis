% Imaging analysis for GC project
%% rigister your image
 edit regtiff
%% load imaging data and ROI
load('data_CNMF.mat')
clearvars -except Y_r
% clear
% filename = 'reg2.tif'
% data=MulRoiANA(filename)
% load C_df data from CNMF result
% %%
% save('data.mat','data')
%% load the event
[analog,trial] = process_intan_v2('RVKC340_180612.rhd');
%% load the licking
analog2 = SettingTwoPhoton('RVKC340_180612');
for i = 1:length(analog2)
    info(i) = analog2p_gc(analog2(i).data,1000);
end
%% check the signal
figure
for i = 1:3
    subplot(1,3,i)
    plot(info(i).time, analog2(i).data(:,3));
    hold on 
    scatter(info(i).lick, 0.2*ones(size(info(i).lick)));
end
%% generate a trial structure; I know each trial contains 60 frames; and licking data are aligned to the tone;
dF = reshape(full(Y_r), size(Y_r,1),60,[]);
for i = 1:length(trial)
    trial(i).licks = info(i).lick - info(i).tone;
    trial(i).trace = dF(:,:,i);
end
save('data.mat','Y_r','trial')
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
for j =101:117
n = j;
figure;
for i = 1: length(trial)
    plot(trial(i).framT, trial(i).traceSmooth_dF(n,:),'Color',[0,0,0,0.5])
    hold on
end
title(['Neuron # ', num2str(n)])
xlim([-1.5,7])

end
%% Plot the Sucrose trial
% for j = 1:size(trial(1).trace,1)
    n = 105;
    figure;
    for i = 1: length(trial)
        if ~isnan(trial(i).S)
        plot(trial(i).framT, trial(i).traceSmooth_dF(n,:),'Color',[0,0,0,0.5])
        hold on
        end
    end
    title(['Neuron # ', num2str(n)])
    xlim([-1.5,7])
% end
%% Plot the Maltose trial
% for j = 1:size(trial(1).trace,1)
    n = 105;
    figure;
    for i = 1: length(trial)
        if ~isnan(trial(i).N)
        plot(trial(i).framT, trial(i).traceSmooth_dF(n,:),'Color',[0,0,0,0.5])
        hold on
        end
    end
    title(['Neuron # ', num2str(n)])
    xlim([-1.5,7])
% end
%% Plot the quinine trial
    n = 105;
    figure;
    for i = 1: length(trial)
        if ~isnan(trial(i).CA)
        plot(trial(i).framT, trial(i).traceSmooth_dF(n,:),'Color',[0,0,0,0.5])
        hold on
        end
    end
    title(['Neuron # ', num2str(n)])
    xlim([-1.5,7])
    
%% Plot the Cyclohexamide trial
    n = 105;
    figure;
    for i = 1: length(trial)
        if ~isnan(trial(i).Q)
        plot(trial(i).framT, trial(i).traceSmooth_dF(n,:),'Color',[0,0,0,0.5])
        hold on
        end
    end
    title(['Neuron # ', num2str(n)])
    xlim([-1.5,7])
    
%% Plot the water trial
    n = 105;
    figure;
    for i = 1: length(trial)
        if ~isnan(trial(i).W)
        plot(trial(i).framT, trial(i).traceSmooth_dF(n,:),'Color',[0,0,0,0.5])
        hold on
        end
    end
    title(['Neuron # ', num2str(n)])
    xlim([-1.5,7])
% end
%% Visualize the averaged response
% for n = 1:27
% n =12;
n =105;
it = 1;
for i = 1: length(trial)
    if ~isnan(trial(i).S)
        S_trace_dF(it,:) = trial(i).traceSmooth_dF(n,:);
        it =1+it;
    end
end

% n =12;
it = 1;
for i = 1: length(trial)
    if ~isnan(trial(i).N)
        M_trace_dF(it,:) = trial(i).traceSmooth_dF(n,:);
        it =1+it;
    end
end

it = 1;
for i = 1: length(trial)
    if ~isnan(trial(i).CA)
        Q_trace_dF(it,:) = trial(i).traceSmooth_dF(n,:);
        it =1+it;
    end
end

it = 1;
for i = 1: length(trial)
    if ~isnan(trial(i).Q)
        Cy_trace_dF(it,:) = trial(i).traceSmooth_dF(n,:);
        it =1+it;
    end
end

figure;
plot(trial(1).framT,mean(S_trace_dF,1))
hold on
plot(trial(1).framT,mean(M_trace_dF,1))
plot(trial(1).framT,mean(Q_trace_dF,1))
plot(trial(1).framT,mean(Cy_trace_dF,1))
xlim([-1.5,7])
legend({'S','M','Q','Cy'})
% end
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
    elseif length(idx) ==49
        idx(end+1) = idx(end)+1;
        trial(i).T = trial(i).Time_Taste(idx);
        trial(i).Taste = trial(i).traceSmooth_dF(:,idx);
    end 
end
%%
figure;
for i = 1:length(trial)
    T(i,:) = trial(i).T;
end
T = mean(T,1);
for i = 1:length(trial)
    trial(i).Tpro = T; % creat a proximate time for all tastant, as tastant may jitter a little bit.
end
%% plot the average response aligned to tastant
n =105;
it = 1;
for i = 1: length(trial)
    if ~isnan(trial(i).S)
        S_Taste_dF(it,:) = trial(i).Taste(n,:);
        it =1+it;
    end
end

% n =12;
it = 1;
for i = 1: length(trial)
    if ~isnan(trial(i).N)
        M_Taste_dF(it,:) = trial(i).Taste(n,:);
        it =1+it;
    end
end

it = 1;
for i = 1: length(trial)
    if ~isnan(trial(i).W)
        W_Taste_dF(it,:) = trial(i).Taste(n,:);
        it =1+it;
    end
end


figure;
plot(T,mean(S_Taste_dF,1))
hold on
plot(T,mean(M_Taste_dF,1))
plot(T,mean(W_Taste_dF,1))
ylabel('dF/F')
xlabel('Time (s)')
xlim([-4,4])
legend('Sucrose','Maltose','Water')
% end
%Plot the Sucrose trial
% for j = 1:size(trial(1).trace,1)
% n = 12;
figure;
for i = 1: size(S_Taste_dF,1)  
    plot(T, S_Taste_dF(i,:),'Color',[0,0,0,0.2])
    hold on
end
plot(T,mean(S_Taste_dF,1),'k')
title(['Neuron # ', num2str(n), ' Sucrose'])
xlim([-4,4])
ylabel('dF/F')
xlabel('Time (s)')   
% Plot the Maltose trial    
% n = 12;
figure;
for i = 1: size(M_Taste_dF,1)  
    plot(T, M_Taste_dF(i,:),'Color',[0,0,0,0.2])
    hold on
end
plot(T,mean(M_Taste_dF,1),'k')
title(['Neuron # ', num2str(n), ' Maltose'])
xlim([-4,4])
ylabel('dF/F')
xlabel('Time (s)')

figure;
for i = 1: size(W_Taste_dF,1)  
    plot(T, W_Taste_dF(i,:),'Color',[0,0,0,0.2])
    hold on
end
plot(T,mean(W_Taste_dF,1),'k')
title(['Neuron # ', num2str(n), ' Water'])
xlim([-4,4])
ylabel('dF/F')
xlabel('Time (s)')
%%
save('data.mat','Y_r','trial')
%% reorganize the data
neuron = trial2neuron5tastant(trial);
%% stats for each tastant
%% statistical test
for j = 1:length(neuron)
    % j = 15;
    idx = find(neuron(j).T>-1 & neuron(j).T <0);
    T_idx1 = find(neuron(j).T>0 & neuron(j).T <3);
%     T_idx2 = find(neuron(j).T>1 & neuron(j).T <2);
%     T_idx3 = find(neuron(j).T>2 & neuron(j).T <3);
    S_baseline    = mean(neuron(j).S_Taste_dF(:,idx),2);
    S_Taste_1     = mean(neuron(j).S_Taste_dF(:,T_idx1),2);
%     S_Taste_2   = mean(neuron(j).S_Taste_dF(:,T_idx2),2);
%     S_Taste_3    = mean(neuron(j).S_Taste_dF(:,T_idx3),2);
    
    [p(1),h(1)] = ranksum(S_baseline,S_Taste_1);
    if mean(S_Taste_1)< mean(S_baseline);
        h(1) = 0;
    end
%     [p(2),h(2)] = ranksum(S_baseline,S_Taste_2);
%     if mean(S_Taste_2)< mean(S_baseline);
%         h(2) = 0;
%     end
%     [p(3),h(3)] = ranksum(S_baseline,S_Taste_3);
%     if mean(S_Taste_3)< mean(S_baseline);
%         h(3) = 0;
%     end
    M_baseline    = mean(neuron(j).M_Taste_dF(:,idx),2); % 2nd taste
    M_Taste_1     = mean(neuron(j).M_Taste_dF(:,T_idx1),2);
    [p(2),h(2)] = ranksum(M_baseline,M_Taste_1);
    if mean(M_Taste_1)< mean(M_baseline);
        h(2) = 0;
    end
    
    CA_baseline    = mean(neuron(j).CA_Taste_dF(:,idx),2); % 3rd taste
    CA_Taste_1     = mean(neuron(j).CA_Taste_dF(:,T_idx1),2);
    [p(3),h(3)] = ranksum(CA_baseline,CA_Taste_1);
    if mean(CA_Taste_1)< mean(CA_baseline);
        h(3) = 0;
    end
    
    Q_baseline    = mean(neuron(j).Q_Taste_dF(:,idx),2); % 4th taste
    Q_Taste_1     = mean(neuron(j).Q_Taste_dF(:,T_idx1),2);
    [p(4),h(4)] = ranksum(Q_baseline,Q_Taste_1);
    if mean(Q_Taste_1)< mean(Q_baseline);
        h(4) = 0;
    end   
    
    W_baseline    = mean(neuron(j).W_Taste_dF(:,idx),2); % 4th taste
    W_Taste_1     = mean(neuron(j).W_Taste_dF(:,T_idx1),2);
    [p(5),h(5)] = ranksum(W_baseline,W_Taste_1);
    if mean(W_Taste_1)< mean(W_baseline);
        h(5) = 0;
    end   
    
    
    neuron(j).Sres = h(1);
    neuron(j).Mres = h(2);
    neuron(j).CAres = h(3);
    neuron(j).Qres = h(4);
    neuron(j).Wres = h(5);
end
%%
plot_dF(116,neuron)

%% statistical test here the baseline is the 1 s before the cue; all tastant are tested together
% for j = 1:length(neuron)
%     % j = 15;
%     idx = find(trial(1).framT>-1 & trial(1).framT <0);
%     T_idx1 = find(trial(1).T>0 & trial(1).T <1);
%     T_idx2 = find(trial(1).T>1 & trial(1).T <2);
%     T_idx3 = find(trial(1).T>2 & trial(1).T <3);
%     for i = 1:length(trial)
%         baseline(i) = mean(trial(i).traceSmooth_dF(j,idx),2);
%         Taste_1(i)    = mean(trial(i).Taste(j,T_idx1),2);
%         Taste_2(i)    = mean(trial(i).Taste(j,T_idx2),2);
%         Taste_3(i)    = mean(trial(i).Taste(j,T_idx3),2);
%     end
%     [p(1),h(1)] = ranksum(baseline,Taste_1);
%     if mean(Taste_1)< mean(baseline);
%         h(1) = 0;
%     end
%     [p(2),h(2)] = ranksum(baseline,Taste_2);
%     if mean(Taste_2)< mean(baseline);
%         h(2) = 0;
%     end
%     [p(3),h(3)] = ranksum(baseline,Taste_3);
%     if mean(Taste_3)< mean(baseline);
%         h(3) = 0;
%     end
%     neuron(j).TasteRes.p = p;
%     neuron(j).TasteRes.h = sum(h);
%     neuron(j).TasteResponse = sum(h);
% end
%%
save('data.mat','trial','Y_r','neuron')