%% analysis of licking parameters for habituation
%% Step 1: load the Intan data file
file = dir('*.rhd');
filename = file.name
data = read_Intan(filename);
%% Only used for habituation with water
[W,~]     = Timing_onset_offset(data.event(5,:), data.ts, 0.5,10,1);
[Tone,~]  = Timing_onset_offset(data.event(6,:), data.ts, 0.5,10,0);
summarydata.Tone  = Tone;
summarydata.W     = W; % Line 5 is W
summarydata.analog = data.analog;
summarydata.ts     = data.ts;
clearvars -except summarydata;
%% check the licking signal
figure;plot(summarydata.ts(1:500000),-summarydata.analog(1:500000)) % the licking signal is inverted, so i take the opposite value of it

thr = input('Set the threshold')
[licks,~] = Timing_onset_offset(-summarydata.analog, summarydata.ts, thr,100,1);
summarydata.licks = licks;
clearvars -except summarydata;
%% get the water event
summarydata.W_ps      = summarydata.W(find(diff(summarydata.W)<1)); % there is an asummption here, the water event after tone shoul be very close, less than 1 s.

%% plot licking raster based on trial sequence: 2 s before the Tone; 8 s after the tone
Tone_lick = spike2eventRasteandPSTH_NP (summarydata.licks,summarydata.Tone, 100, -2000, 8000);
figure;
subplot(2,1,1)
for i = 1:length(summarydata.Tone)
    if ~isempty(Tone_lick.spikeraster(i))
        scatter(Tone_lick.spikeraster(i).times,i*ones(size(Tone_lick.spikeraster(i).times)),6,'rv','filled');
        hold on
    end  
end
rectangle('Position',[0,0,2,i],'FaceColor',[0.5 .5 .5,0.5],'EdgeColor',[0.5 .5 .5,0])
xlim([-2,8])
xlabel('Time (s)')
ylabel('Trial #')
set(gca,'TickDir','out');

subplot(2,1,2)
plot(Tone_lick.timepoint,Tone_lick.FR_avg,'r')
rectangle('Position',[0,0,2,max(Tone_lick.FR_avg)],'FaceColor',[0.5 .5 .5,0.5],'EdgeColor',[0.5 .5 .5,0])
xlim([-2,8])
xlabel('Time (s)')
ylabel('Lick Rate (Hz)')
set(gca,'TickDir','out')
%%
for i = 1:length(summarydata.Tone)
    if isempty(Tone_lick.spikeraster(i))
        licks(i) = 0;
    else
        licks(i) = length(find(Tone_lick.spikeraster(i).times > 0)); % only licks afte the tone are counted.
    end
end
figure;
plot(licks,'-o','Color',[1,0,0])
xlabel('Trial #')
ylabel('Licks per trial') % total number of licks between 0 and 8 s after the tone