%% data analysis by Intan
% add directory
addpath('C:\Users\ke-roberto\Documents\MATLAB\Code_Ke\')
%% load data
file = dir('*.rhd');
data = read_Intan(file.name);
%% extract the events
thr = 2;
[CentraLicks,~] = Timing_onset_offset(data.analog(1,:), data.ts, thr,30,0); % get the central licks

message1='Which line is used as the water?';
uiwait(msgbox(message1))
id = input('Which line is used as the water?')
[W,~] = Timing_onset_offset(data.event(id,:), data.ts, 0.5,30,0); % get the central licks


%% plot licking raster based on trial sequence: 2 s before the Tone; 8 s after the tone
Tone_lick = spike2eventRasteandPSTH_NP (CentraLicks,W, 100, 0, 8000);
figure;
subplot(2,1,1)
for i = 1:length(W)
    if ~isempty(Tone_lick.spikeraster(i))
        scatter(Tone_lick.spikeraster(i).times,i*ones(size(Tone_lick.spikeraster(i).times)),6,'rv','filled');
        hold on
    end  
end
% rectangle('Position',[0,0,2,i],'FaceColor',[0.5 .5 .5,0.5],'EdgeColor',[0.5 .5 .5,0])
xlim([0,8])
xlabel('Time (s)')
ylabel('Trial #')
set(gca,'TickDir','out');

subplot(2,1,2)
plot(Tone_lick.timepoint,Tone_lick.FR_avg,'r')
xlim([0,8])
xlabel('Time (s)')
ylabel('Lick Rate (Hz)')
set(gca,'TickDir','out')
%%
