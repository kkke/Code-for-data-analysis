function plotLickIntan(summarydata)
% Input: summarydata
%        summarydata is save by runing beha_2p.m
Tone_lick = spike2eventRasteandPSTH_NP (summarydata.licks,summarydata.Tone, 100, -2000, 8000);
figure
subplot(2,1,1)
event = {'S_ps','N_ps','CA_ps','Q_ps','W_ps'};
color = {'cv','mv','bv','rv','kv'};
k= 1;
for j = 1:length(event)
    for i = 1:length(summarydata.(event{j}))
        id = find(summarydata.Tone<summarydata.(event{j})(i) & summarydata.Tone>summarydata.(event{j})(i)-5);
        if ~isempty(id)
            h(j) = scatter(Tone_lick.spikeraster(id).times,k*ones(size(Tone_lick.spikeraster(id).times)),6,color{j},'filled');
            k = k+1;
        end
        hold on
    end
end

legend(h,{'S','M','Q','CA','W'})
rectangle('Position',[0,0,2,k],'FaceColor',[0.5 .5 .5,0.5],'EdgeColor',[0.5 .5 .5,0])
xlim([-2,8])
xlabel('Time (s)')
ylabel('Trial #')
set(gca,'TickDir','out');


taste = {'S','N','CA','Q','W'};
color = {'c','m','b','r','k'};
subplot(2,1,2)
for i = 1:length(event)
    for j = 1:length(summarydata.(event{i}))
        id = find(summarydata.Tone<summarydata.(event{i})(j) & summarydata.Tone>summarydata.(event{i})(j)-5);
        lick.(taste{i})(j,:) = Tone_lick.scmatrix(id,:);
    end
    avg = sum(lick.(taste{i}),1)./(j*Tone_lick.Binsize*0.001);
    plot(Tone_lick.timepoint,avg,color{i})
    hold on
end
ylabel('Licking rate')
xlabel('Time (s)')
legend({'S','N','CA','Q','W'})

%%
clear lick
t =1;
taste = {'S','N','CA','Q','W'};
event = {'S_ps','N_ps','CA_ps','Q_ps','W_ps'};
figure
subplot(2,1,1)
for i = 1:length(taste)
    lick.(taste{i}) = spike2eventRasteandPSTH_NP(summarydata.licks,summarydata.(event{i}), 100, -3000, 5000);
    for j = 1:length(summarydata.(event{i}))
        if ~isempty(lick.(taste{i}).spikeraster(j).times)
            id = find(lick.(taste{i}).spikeraster(j).times>0);
            licktime.(taste{i})(j) = lick.(taste{i}).spikeraster(j).times(id(1))+summarydata.(event{i})(j);
        end
    end
end
%%
color = {'cv','mv','bv','rv','kv'};
k =1;
for i = 1:length(event)
    for j = 1:length(summarydata.(event{i}))
        if ~isempty(lick.(taste{i}).spikeraster(j).times)
            h(i) = scatter(lick.(taste{i}).spikeraster(j).times,k*ones(size(lick.(taste{i}).spikeraster(j).times)),6,color{i},'filled');
            hold on
            k = k+1;
        end
    end
end
xlabel('Time (s)')
ylabel('Trial #')
xlim([-3,5])

subplot(2,1,2)
color = {'c','m','b','r','k'};

for i = 1:length(taste)
    plot(lick.(taste{i}).timepoint,lick.(taste{i}).FR_avg,color{i})
    hold on
end
xlabel('Time (s)')
ylabel('Licking rate (Hz)')
xlim([-3,5])

%%
color = {'cv','mv','bv','rv','kv'};
figure;
subplot(2,1,1)
k = 1;
for i = 1:length(taste)
    lick.(taste{i}) = spike2eventRasteandPSTH_NP(summarydata.licks,licktime.(taste{i}), 100, -3000, 5000);
    for j = 1:length(summarydata.(event{i}))
        if ~isempty(lick.(taste{i}).spikeraster(j).times)
            h1(i) = scatter(lick.(taste{i}).spikeraster(j).times,k*ones(size(lick.(taste{i}).spikeraster(j).times)),6,color{i},'filled');
            hold on
            k = k+1;
        end
    end
end
xlabel('Time (s)')
ylabel('Trial #')
xlim([-3,5])
subplot(2,1,2)
color = {'c','m','b','r','k'};
for i = 1:length(taste)
    plot(lick.(taste{i}).timepoint,lick.(taste{i}).FR_avg,color{i})
    hold on
end
xlabel('Time (s)')
ylabel('Licking rate (Hz)')
xlim([-3,5])

for i = 1:length(taste)
    for j = 1:length(lick.(taste{i}).spikeraster)
        lickcount.(taste{i}).counts(j,1)= length(find(lick.(taste{i}).spikeraster(j).times>=0 & lick.(taste{i}).spikeraster(j).times<1));
        lickcount.(taste{i}).counts(j,2)= length(find(lick.(taste{i}).spikeraster(j).times>=0 & lick.(taste{i}).spikeraster(j).times<2));
        lickcount.(taste{i}).counts(j,3)= length(find(lick.(taste{i}).spikeraster(j).times>=0 & lick.(taste{i}).spikeraster(j).times<3));
    end
end
for i = 1:length(taste)
    m1(i,:) = mean(lickcount.(taste{i}).counts);
    sem(i,:) = std(lickcount.(taste{i}).counts)./sqrt(size(lickcount.(taste{i}).counts,1));    
end

figure
bar(m1')
ylim([0,25])
ylabel('Lick Count')

% hold on
% errorbar(sem')
% figure;
% plot()


% figure;
% [h1,x1] = ecdf(lickrate_s);
% hold on
% [h2,x2] = ecdf(lickrate_n);
% [h3,x3] = ecdf(lickrate_c);
% [h4,x4]= ecdf(lickrate_q);
% [h5,x5]= ecdf(lickrate_w);
% 
% plot(x1,h1,'c')
% hold on
% plot(x2,h2,'m')
% plot(x3,h3,'b')
% plot(x4,h4,'r')
% plot(x5,h5,'k')
