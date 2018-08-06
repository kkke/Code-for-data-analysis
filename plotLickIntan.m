function plotLickIntan(summarydata)
% Input: summarydata
%        summarydata is save by runing beha_2p.m
Tone_lick = spike2eventRasteandPSTH_NP (summarydata.licks,summarydata.Tone, 100, -2000, 8000);
figure
subplot(2,1,1)
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

subplot(2,1,2)
for j = 1:length(summarydata.S_ps)
    id = find(summarydata.Tone<summarydata.S_ps(j) & summarydata.Tone>summarydata.S_ps(j)-5);
    S_lick(j,:) = Tone_lick.scmatrix(id,:);  
end
S_lick_avg = sum(S_lick,1)./(j*Tone_lick.Binsize*0.001);
plot(Tone_lick.timepoint,S_lick_avg,'c')
hold on

for j = 1:length(summarydata.N_ps)
    id = find(summarydata.Tone<summarydata.N_ps(j) & summarydata.Tone>summarydata.N_ps(j)-5);
    N_lick(j,:) = Tone_lick.scmatrix(id,:);  
end
N_lick_avg = sum(N_lick,1)./(j*Tone_lick.Binsize*0.001);
plot(Tone_lick.timepoint,N_lick_avg,'m')

for j = 1:length(summarydata.CA_ps)
    id = find(summarydata.Tone<summarydata.CA_ps(j) & summarydata.Tone>summarydata.CA_ps(j)-5);
    CA_lick(j,:) = Tone_lick.scmatrix(id,:);  
end
CA_lick_avg = sum(CA_lick,1)./(j*Tone_lick.Binsize*0.001);
plot(Tone_lick.timepoint,CA_lick_avg,'b')

for j = 1:length(summarydata.Q_ps)
    id = find(summarydata.Tone<summarydata.Q_ps(j) & summarydata.Tone>summarydata.Q_ps(j)-5);
    Q_lick(j,:) = Tone_lick.scmatrix(id,:);  
end
Q_lick_avg = sum(Q_lick,1)./(j*Tone_lick.Binsize*0.001);
plot(Tone_lick.timepoint,Q_lick_avg,'r')

for j = 1:length(summarydata.W_ps)
    id = find(summarydata.Tone<summarydata.W_ps(j) & summarydata.Tone>summarydata.W_ps(j)-5);
    W_lick(j,:) = Tone_lick.scmatrix(id,:);  
end
W_lick_avg = sum(W_lick,1)./(j*Tone_lick.Binsize*0.001);
plot(Tone_lick.timepoint,W_lick_avg,'k')
ylabel('Licking rate')
xlabel('Time (s)')
legend({'S','N','CA','Q','W'})
