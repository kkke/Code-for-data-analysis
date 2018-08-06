function resp = tasteResponse2(neuron)
t = 0.05;
rw = 3;
for j = 1:length(neuron)
    % j = 15;
    idx = find(neuron(j).T>-1 & neuron(j).T <0);
    T_idx1 = find(neuron(j).T>0 & neuron(j).T <1);
    T_idx2 = find(neuron(j).T>1 & neuron(j).T <2);
    T_idx3 = find(neuron(j).T>2 & neuron(j).T <3);
    S_baseline    = mean(neuron(j).S_Taste_dF(:,idx),2);
    S_Taste_1     = mean(neuron(j).S_Taste_dF(:,T_idx1),2);
    S_Taste_2     = mean(neuron(j).S_Taste_dF(:,T_idx2),2);
    S_Taste_3     = mean(neuron(j).S_Taste_dF(:,T_idx3),2);
    S_Taste       = [S_Taste_1, S_Taste_2, S_Taste_3];
    for i = 1:size(S_Taste,2)
        [p(i),h(i)] = ranksum(S_baseline,S_Taste(:,i));
        if p(i)>0.05/3 || mean(S_Taste(:,i))<0 || mean(S_Taste(:,i))< mean(S_baseline) 
            h(i) = 0;
        end
    end
    if sum(h)>0
        s_resp = 1;
    else
        s_resp = 0;
    end
   
    M_baseline    = mean(neuron(j).M_Taste_dF(:,idx),2); % 2nd taste
    M_Taste_1     = mean(neuron(j).M_Taste_dF(:,T_idx1),2);
    M_Taste_2     = mean(neuron(j).M_Taste_dF(:,T_idx2),2);
    M_Taste_3     = mean(neuron(j).M_Taste_dF(:,T_idx3),2);
    M_Taste       = [M_Taste_1, M_Taste_2, M_Taste_3];
    for i = 1:size(M_Taste,2)
        [p(i),h(i)] = ranksum(M_baseline,M_Taste(:,i));
        if p(i)>0.05/3 || mean(M_Taste(:,i))<0 || mean(M_Taste(:,i))< mean(M_baseline) 
            h(i) = 0;
        end
    end
    if sum(h)>0
        m_resp = 1;
    else
        m_resp = 0;
    end
    
    CA_baseline    = mean(neuron(j).CA_Taste_dF(:,idx),2); % 3rd taste
    CA_Taste_1     = mean(neuron(j).CA_Taste_dF(:,T_idx1),2);
    CA_Taste_2     = mean(neuron(j).CA_Taste_dF(:,T_idx2),2);
    CA_Taste_3     = mean(neuron(j).CA_Taste_dF(:,T_idx3),2);
    CA_Taste       = [CA_Taste_1, CA_Taste_2, CA_Taste_3];
    for i = 1:size(CA_Taste,2)
        [p(i),h(i)] = ranksum(CA_baseline,CA_Taste(:,i));
        if p(i)>0.05/3 || mean(CA_Taste(:,i))<0 || mean(CA_Taste(:,i))< mean(CA_baseline) 
            h(i) = 0;
        end
    end
    if sum(h)>0
        ca_resp = 1;
    else
        ca_resp = 0;
    end
    
    Q_baseline    = mean(neuron(j).Q_Taste_dF(:,idx),2); % 4th taste
    Q_Taste_1     = mean(neuron(j).Q_Taste_dF(:,T_idx1),2);
    Q_Taste_2     = mean(neuron(j).Q_Taste_dF(:,T_idx2),2);
    Q_Taste_3     = mean(neuron(j).Q_Taste_dF(:,T_idx3),2);
    Q_Taste       = [Q_Taste_1, Q_Taste_2, Q_Taste_3];
    for i = 1:size(Q_Taste,2)
        [p(i),h(i)] = ranksum(Q_baseline,Q_Taste(:,i));
        if p(i)>0.05/3 || mean(Q_Taste(:,i))<0 || mean(Q_Taste(:,i))< mean(Q_baseline) 
            h(i) = 0;
        end
    end
    if sum(h)>0
        q_resp = 1;
    else
        q_resp = 0;
    end   
    
    W_baseline    = mean(neuron(j).W_Taste_dF(:,idx),2); % 4th taste
    W_Taste_1     = mean(neuron(j).W_Taste_dF(:,T_idx1),2);
    W_Taste_2     = mean(neuron(j).W_Taste_dF(:,T_idx2),2);
    W_Taste_3     = mean(neuron(j).W_Taste_dF(:,T_idx3),2);
    W_Taste       = [W_Taste_1, W_Taste_2, W_Taste_3];
    for i = 1:size(W_Taste,2)
        [p(i),h(i)] = ranksum(W_baseline,W_Taste(:,i));
        if p(i)>0.05/3 || mean(W_Taste(:,i))<0 || mean(W_Taste(:,i))< mean(W_baseline) 
            h(i) = 0;
        end
    end
    if sum(h)>0
        w_resp = 1;
    else
        w_resp = 0;
    end  
    resp(j).Sres = s_resp;
    resp(j).Mres = m_resp;
    resp(j).CAres = ca_resp;
    resp(j).Qres = q_resp;
    resp(j).Wres = w_resp;
end