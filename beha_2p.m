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
summarydata.filename = filename(1:end-4);
summarydata.Frame = Frame;
summarydata.Tone  = Tone;
summarydata.W     = W; % Line 5 is W
summarydata.Q     = Q; % Line 4 is Q
summarydata.CA    = CA; % Line 3 is CA
summarydata.N     = N; % line 2 is N
summarydata.S     = S; % Line1 is S;
summarydata.analog = data.analog;
summarydata.ts     = data.ts;
clearvars -except summarydata;
%% check the licking signal
figure;plot(summarydata.ts(1:500000),-summarydata.analog(1:500000)) % the licking signal is inverted, so i take the opposite value of it

thr = input('Set the threshold')
[licks,~] = Timing_onset_offset(-summarydata.analog, summarydata.ts, thr,100,1);
summarydata.licks = licks;
clearvars -except summarydata;
%%
% S  --> Line1;
% N  --> Line2;
% CA --> Line3; 
% Q  --> Line4; 
% W  --> Line5:
summarydata.W_ps      = summarydata.W(find(diff(summarydata.W)<1)); % there is an asummption here, the water event after tone shoul be very close, less than 1 s.
summarydata.S_ps      = summarydata.S(1:2:end); % there is an asummption here, there are also two consecutive event, take the first event;
summarydata.N_ps      = summarydata.N(1:2:end);
summarydata.CA_ps      = summarydata.CA(1:2:end);
summarydata.Q_ps      = summarydata.Q(1:2:end);
%% trying to extract which trial deliver which tastant
L1 = [ones(size(summarydata.S_ps));summarydata.S_ps];
L2 = [2*ones(size(summarydata.N_ps));summarydata.N_ps];
L3 = [3*ones(size(summarydata.CA_ps));summarydata.CA_ps];
L4 = [4*ones(size(summarydata.Q_ps));summarydata.Q_ps];
L5 = [5*ones(size(summarydata.W_ps));summarydata.W_ps];
L = [L1, L2, L3, L4, L5];
[~,I] = sort(L(2,:));
L_trial = L(:,I);
save('summary.mat','summarydata')
%% plot licking raster based on trial sequence
Tone_lick = spike2eventRasteandPSTH_NP (summarydata.licks,summarydata.Tone, 100, -2000, 8000);
figure;
subplot(2,1,1)

for i = 1:length(Tone_lick.spikeraster)
    if ~isempty(Tone_lick.spikeraster(i))
        id = L_trial(1,find(L_trial(2,:)>summarydata.Tone(i) & L_trial(2,:)<summarydata.Tone(i)+5));
        if ~isempty(id)
            switch id
                case 1
                    scatter(Tone_lick.spikeraster(i).times,i*ones(size(Tone_lick.spikeraster(i).times)),6,'cv','filled');
                case 2
                    scatter(Tone_lick.spikeraster(i).times,i*ones(size(Tone_lick.spikeraster(i).times)),6,'mv','filled');
                case 3
                    scatter(Tone_lick.spikeraster(i).times,i*ones(size(Tone_lick.spikeraster(i).times)),6,'bv','filled')
                case 4
                    scatter(Tone_lick.spikeraster(i).times,i*ones(size(Tone_lick.spikeraster(i).times)),6,'rv','filled')
                case 5
                    scatter(Tone_lick.spikeraster(i).times,i*ones(size(Tone_lick.spikeraster(i).times)),6,'kv','filled')
            end
        end
        
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
box off

%% plot the licking raster based on different tastant
figure
for i = 1:length(summarydata.S_ps)
    id = find(summarydata.Tone<summarydata.S_ps(i) & summarydata.Tone>summarydata.S_ps(i)-5);
    if ~isempty(id)
        h1 = scatter(Tone_lick.spikeraster(id).times,i*ones(size(Tone_lick.spikeraster(id).times)),6,'cv','filled');
    end   
    hold on
end
for j = 1:length(summarydata.N_ps)
    id = find(summarydata.Tone<summarydata.N_ps(j) & summarydata.Tone>summarydata.N_ps(j)-5);
    i = i+1;
    if ~isempty(id)
        h2 = scatter(Tone_lick.spikeraster(id).times,i*ones(size(Tone_lick.spikeraster(id).times)),6,'mv','filled');
    end   
end
for j = 1:length(summarydata.CA_ps)
    id = find(summarydata.Tone<summarydata.CA_ps(j) & summarydata.Tone>summarydata.CA_ps(j)-5);
    i = i+1;
    if ~isempty(id)
        h3 = scatter(Tone_lick.spikeraster(id).times,i*ones(size(Tone_lick.spikeraster(id).times)),6,'bv','filled');
    end   
end
for j = 1:length(summarydata.Q_ps)
    id = find(summarydata.Tone<summarydata.Q_ps(j) & summarydata.Tone>summarydata.Q_ps(j)-5);
    i = i+1;
    if ~isempty(id)
        h4 = scatter(Tone_lick.spikeraster(id).times,i*ones(size(Tone_lick.spikeraster(id).times)),6,'rv','filled');
    end   
end
for j = 1:length(summarydata.W_ps)
    id = find(summarydata.Tone<summarydata.W_ps(j) & summarydata.Tone>summarydata.W_ps(j)-5);
    i = i+1;
    if ~isempty(id)
        h5 = scatter(Tone_lick.spikeraster(id).times,i*ones(size(Tone_lick.spikeraster(id).times)),6,'kv','filled');
    end   
end
legend([h1,h2,h3,h4,h5],{'S','M','Q','C','W'})
rectangle('Position',[0,0,2,i],'FaceColor',[0.5 .5 .5,0.5],'EdgeColor',[0.5 .5 .5,0])
xlim([-2,8])
xlabel('Time (s)')
ylabel('Trial #')
set(gca,'TickDir','out');
%%
save('summary.mat','summarydata')
