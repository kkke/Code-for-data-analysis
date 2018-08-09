function plot_dF_lick(i,neuron)
% i = 25;
figure
for j = 1:size(neuron(i).Licktrace_dF,1)
    plot(neuron(i).TLick,neuron(i).Licktrace_dF(j,:),'Color',[0.5,0.5,0.5,0.5])
    hold on
end
plot(neuron(i).TLick,mean(neuron(i).Licktrace_dF,1),'r')
xlim([-2,4])
xlabel('Time (s)')
ylabel('dF/F')
title(['Lick response Neuron #',num2str(i)])

% figure
% for j = 1:size(neuron(i).M_Taste_dF,1)
%     plot(neuron(i).T,neuron(i).M_Taste_dF(j,:),'Color',[0.5,0.5,0.5,0.5])
%     hold on
% end
% plot(neuron(i).T,mean(neuron(i).M_Taste_dF,1))
% xlim([-2,4])
% xlabel('Time (s)')
% ylabel('dF/F')
% title('Maltose')

% figure
% for j = 1:size(neuron(i).CA_Taste_dF,1)
%     plot(neuron(i).T,neuron(i).CA_Taste_dF(j,:),'Color',[0.5,0.5,0.5,0.5])
%     hold on
% end
% plot(neuron(i).T,mean(neuron(i).CA_Taste_dF,1))
% xlim([-2,4])
% xlabel('Time (s)')
% ylabel('dF/F')
% title('Quinine')

% figure
% for j = 1:size(neuron(i).Q_Taste_dF,1)
%     plot(neuron(i).T,neuron(i).Q_Taste_dF(j,:),'Color',[0.5,0.5,0.5,0.5])
%     hold on
% end
% plot(neuron(i).T,mean(neuron(i).Q_Taste_dF,1))
% xlim([-2,4])
% xlabel('Time (s)')
% ylabel('dF/F')
% title('Cycloheximide')

% figure
% for j = 1:size(neuron(i).W_Taste_dF,1)
%     plot(neuron(i).T,neuron(i).W_Taste_dF(j,:),'Color',[0.5,0.5,0.5,0.5])
%     hold on
% end
% plot(neuron(i).T,mean(neuron(i).W_Taste_dF,1))
% xlim([-2,4])
% xlabel('Time (s)')
% ylabel('dF/F')
% title('Water')
% 
% figure
% plot(neuron(i).T,mean(neuron(i).S_Taste_dF,1))
% hold on
% plot(neuron(i).T,mean(neuron(i).M_Taste_dF,1))
% % plot(neuron(i).T,mean(neuron(i).CA_Taste_dF,1))
% % plot(neuron(i).T,mean(neuron(i).Q_Taste_dF,1))
% plot(neuron(i).T,mean(neuron(i).W_Taste_dF,1))
% % legend('S','M','Q','Cyx','W')
% legend('S','M','W')
% xlim([-2,4])
% xlabel('Time (s)')
% ylabel('dF/F')
% title('Taste response')


