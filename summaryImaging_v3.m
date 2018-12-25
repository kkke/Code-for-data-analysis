%%
summaryImaging_v2
[stats,stats2] = TasteTestDev(all,3.5)
%%
ind      = find([stats2.taste]==1);
ind_cue  = find([stats2.cue]==1);
ind_lick = find([stats2.lick]==1);
c = intersect(ind,ind_cue);
taste = stats.taste(ind,:);

tuning = sum(taste,2);

% figure;
% bar(sum(taste)./length(taste))
% legend('S','M','CA','Q','W')
% ylim([0,1])
for i = 1:5
    res(i) = length(find(tuning ==i))
end
figure;
plot(res./length(taste),'-^')
legend('2 SD')
xlim([0.5,5.5])

%%