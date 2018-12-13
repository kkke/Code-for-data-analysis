function p = tasteResponse5(neuron,trial)
t = 0.05;
rw = 3;
for j = 1:length(neuron)
    % j = 15;
    idx = find(neuron(j).T>0 & neuron(j).T <rw);
    S_Taste    = mean(neuron(j).S_Taste_dF(:,idx),2);
    M_Taste    = mean(neuron(j).M_Taste_dF(:,idx),2);
    CA_Taste   = mean(neuron(j).CA_Taste_dF(:,idx),2);
    Q_Taste    = mean(neuron(j).Q_Taste_dF(:,idx),2);
    W_Taste    = mean(neuron(j).W_Taste_dF(:,idx),2);
    Taste      = [S_Taste; M_Taste; CA_Taste; Q_Taste; W_Taste];
    group      = [zeros(size(S_Taste)); ones(size(M_Taste)); 2*ones(size(CA_Taste))...
                  ;3*ones(size(Q_Taste)); 4*ones(size(W_Taste))];
    [p(j),tbl,stats] = kruskalwallis(Taste,group);

end
%% cue response with signed rank
