function barplot_equal(control)
m1 = mean(control(:,1));
m2 = mean(control(:,2));
m3 = mean(control(:,3));

name = {[],'Central lick','Left Lick','Right Lick',[]}
sem1 = std(control(:,1))./sqrt(length(control(:,3)));
sem2 = std(control(:,2))./sqrt(length(control(:,3)));
sem3 = std(control(:,3))./sqrt(length(control(:,3)));

figure;
bar(1,m1,'EdgeColor',[0 0 0],'FaceColor',[0.5,0.5,1],'LineWidth',1);hold on
bar(2,m2,'EdgeColor',[0 0 0],'FaceColor',[1,0.5,0.25],'LineWidth',1);hold on
bar(3,m3,'EdgeColor',[0 0 0],'FaceColor',[1,0.5,0.25],'LineWidth',1);hold on
set(gca,'XTick',0:4)
set(gca,'xticklabel',name)
ylabel('Licking Rate (Hz)')
hold on;
e=errorbar([m1, m2, m3], [sem1, sem2, sem3],'LineWidth',1.5);
% plot([1 1],[mean(a_m) mean(a_m)+a_sem],'r');hold on;
e.Color ='k';
e.LineStyle = 'none';