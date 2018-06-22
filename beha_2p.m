%% behaviral analysis for 2p rig
%% Step 1: load data from Intan
file = dir('*.rhd');
filename = file.name
data = read_Intan(filename);
[Frame,~] = Timing_onset_offset(data.event(7,:), data.ts, 0.5,10,1);
[Tone,~]  = Timing_onset_offset(data.event(6,:), data.ts, 0.5,10,1);
[W,~]     = Timing_onset_offset(data.event(5,:), data.ts, 0.5,10,1);
[Q,~]     = Timing_onset_offset(data.event(4,:), data.ts, 0.5,10,0);
[CA,~]    = Timing_onset_offset(data.event(3,:), data.ts, 0.5,10,0);
[N,~]     = Timing_onset_offset(data.event(2,:), data.ts, 0.5,10,0);
[S,~]     = Timing_onset_offset(data.event(1,:), data.ts, 0.5,10,1);
W_ps      = W(find(diff(W)<1)); % there is an asummption here, the water event after tone shoul be very close, less than 1 s.
summarydata.filename = filename(1:end-4);
summarydata.Frame = Frame;
summarydata.Tone  = Tone;
summarydata.W     = W; % Line 5 is W
summarydata.Q     = Q; % Line 4 is Q
summarydata.CA    = CA; % Line 3 is CA
summarydata.N     = N; % line 2 is N
summarydata.S     = S; % Line1 is S;
summarydata.W_ps   = W_ps;
summarydata.analog = data.analog;
summarydata.ts     = data.ts;
clearvars -except summarydata;
%% check the licking signal
figure;plot(summarydata.ts(1:500000),-summarydata.analog(1:500000))
thr = -0.2
[licks,~] = Timing_onset_offset(-summarydata.analog, summarydata.ts, thr,100,1);
%% visualize the licking
figure;plot(summarydata.ts(1:500000),-summarydata.analog(1:500000))
hold on
scatter(licks)
%% 
summarydata.licks = licks;
clearvars -except summarydata;
%% Step 2: load data from MScan
% analog signal should be extracted by manually export data to txt file
% ANSI coding
% file = dir('*.txt');
% filename = file.name
% analog = importtext(filename);
% analog2.data = analog(:,2:4);
% info = analog2p_gc(analog2.data,1000);
% figure
% plot(info.time, analog2.data(:,3));
% hold on 
% scatter(info.lick, 0.2*ones(size(info.lick)));
% summarydata.info     = info;  % data acquired by scope
% clearvars -except summarydata
%% check the difference of the Tone recorded between Intan and Scope' there is slowly drift of tone; 30 min 20 ms drifts.
% comp      = summarydata.Tone - summarydata.info.tone;
% diff_toen = max(abs(diff(comp))) 
%% generate trial structure to align data
% summarydata.licks = summarydata.info.lick+mean(comp);
%% Plot the PSTH of licks aligned to the onset of tone
% Tone_lick = spike2eventRasteandPSTH_NP (summarydata.licks,summarydata.Tone, 100, -2000, 8000);
% figure;
% subplot(2,1,1)
% for i = 1:length(Tone_lick.spikeraster)
%     if ~isempty(Tone_lick.spikeraster(i))        
%             h1=scatter(Tone_lick.spikeraster(i).times,i*ones(size(Tone_lick.spikeraster(i).times)),6,'rv','filled');
%             hold on
%     end
% end
% rectangle('Position',[0,0,2,i],'FaceColor',[0.5 .5 .5,0.5],'EdgeColor',[0.5 .5 .5,0])
% xlim([-2,8])
% xlabel('Time (s)')
% ylabel('Trial #')
% set(gca,'TickDir','out');
% 
% subplot(2,1,2)
% plot(Tone_lick.timepoint,Tone_lick.FR_avg,'r')
% rectangle('Position',[0,0,2,max(Tone_lick.FR_avg)],'FaceColor',[0.5 .5 .5,0.5],'EdgeColor',[0.5 .5 .5,0])
% xlim([-2,8])
% xlabel('Time (s)')
% ylabel('Lick Rate (Hz)')
% set(gca,'TickDir','out')
% box off
%%
save('summary.mat','summarydata')


