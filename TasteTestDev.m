function [stats,stats2] = TasteTestDev(all)
data = all;
tasteResponse = {'S_Taste_dF','M_Taste_dF','CA_Taste_dF','Q_Taste_dF','W_Taste_dF'};
trace = {'S_trace_dF','M_trace_dF','CA_trace_dF','Q_trace_dF','W_trace_dF'};

taste         = {'S', 'M', 'CA','Q','W'};
thr = 3;
for j = 1:length(taste)
    for i = 1:length(data)
        ind = find(data(i).framT>-2 & data(i).framT<-1);  % choose [-2,-1]  as baseline; [-1,0] was used to calculate dF
        baseline.(taste{j}) = data(i).(trace{j})(:,ind);
        m = mean(baseline.(taste{j}),2);
        indT = find(data(i).T>0 & data(i).T<3.5); % [0 3.5] after taste delivery was used to for stats
        tasteResAvg = data(i).(tasteResponse{j})(:,indT);
        [M, I] = max(tasteResAvg,[],2); % get the peak for each trial
        for z = 1:length(I) % loop through each trial
            if I(z)<=2; % the peak is the first two frames after the taste
                indTest(z,:) =indT(1:5);
            else
                indTest(z,:) = (I(z)-2+indT(1)-1):1:(I(z)+2 +indT(1)-1); % 5 frames around the peak
            end
            tasteTest(z,:) = data(i).(tasteResponse{j})(z,indTest(z,:)); % get 5 frames around the peak.
        end
        tasteTestavg = mean(mean(tasteTest,2));
        
        
        if isempty(find(tasteTestavg>mean(m)+thr*std(m))) % only test the excitatory response
            resp(i).(taste{j}) = 0;
        else
            resp(i).(taste{j}) = 1;
        end
        clear indTest
    end
end
%%
stats.taste = cell2mat(squeeze(struct2cell(resp))');
stats.cue   = [all.CueRes];
stats.lick  = [all.LickRes];
stats.tasteSum = sum(stats.taste,2);
for i = 1:length(stats.cue)
    if stats.cue(i) == 1
        stats2.cue(i) = 1
        stats2.lick(i) =0;
    elseif stats.cue(i) == 0 &stats.lick(i) == 1
        stats2.cue(i) = 0;
        stats2.lick(i) = 1;
    else
        stats2.cue(i)  =0;
        stats2.lick(i) =0;
    end
end
for i = 1:length(stats.lick)
    if stats.lick(i) == 1
        stats2.taste(i) =0;
    elseif stats.cue(i) == 0 & stats.tasteSum(i) >= 1
        stats2.taste(i) = 1;
        
    else
        stats2.taste(i) =0;
    end
end
% a=[];
% for i = 1:length(stats.lick)
%     if stats.lick(i) ==1 & stats.tasteSum(i) >=1
%         a = [a,i];
%     end
% end


%% tuning curve
% ind = find([stats2.taste]==1);
% taste = stats.taste(ind,:);
% 
% tuning = sum(taste,2);
% 
% figure;
% bar(sum(taste)./length(taste))
% legend('S','M','CA','Q','W')
% ylim([0,1])
% for i = 1:5
%     res(i) = length(find(tuning ==i))
% end
% figure;
% plot(res./length(taste),'-o')
% legend('3 SD')
% xlim([0.5,5.5])
%%
% data = all;
% tasteResponse = {'S_Taste_dF','M_Taste_dF','CA_Taste_dF','Q_Taste_dF','W_Taste_dF'};
% taste         = {'S', 'M', 'CA','Q','W'};
% thr = 3;
% for j = 1:length(taste)
%     for i = 1:length(data)
%         ind = find(data(i).T>-1 & data(i).T<0);
%         baseline.(taste{j}) = data(i).(tasteResponse{j})(:,ind);
%         m = mean(baseline.(taste{j}),2);
%         indT = find(data(i).T>0 & data(i).T<3.5);
%         tasteResAvg = mean(data(i).(tasteResponse{j})(:,indT),1);
%         if isempty(find(tasteResAvg>mean(m)+thr*std(m)))
%             resp(i).(taste{j}) = 0;
%         else
%             resp(i).(taste{j}) = 1;
%         end
%     end
% end