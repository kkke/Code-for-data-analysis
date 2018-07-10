%% this code is to reorganize the trial structure into neuron structure

function neuron = trial2neuron5tastant(trial);
for n = 1:size(trial(1).trace,1)
    
 it = 1;
for i = 1: length(trial)
    if ~isnan(trial(i).S)
        S_trace_dF(it,:) = trial(i).traceSmooth_dF(n,:);
        S_licks{it}    = trial(i).licks - trial(i).S(1);
        it =1+it;
    end
end

% n =12;
it = 1;
for i = 1: length(trial)
    if ~isnan(trial(i).N)
        M_trace_dF(it,:) = trial(i).traceSmooth_dF(n,:);
        M_licks{it}    = trial(i).licks - trial(i).N(1);
        it =1+it;
    end
end

% n =12;
it = 1;
for i = 1: length(trial)
    if ~isnan(trial(i).CA)
        CA_trace_dF(it,:) = trial(i).traceSmooth_dF(n,:);
        CA_licks{it}    = trial(i).licks - trial(i).CA(1);
        it =1+it;
    end
end

% n =12;
it = 1;
for i = 1: length(trial)
    if ~isnan(trial(i).Q)
        Q_trace_dF(it,:) = trial(i).traceSmooth_dF(n,:);
        Q_licks{it}    = trial(i).licks - trial(i).Q(1);
        it =1+it;
    end
end

it = 1;
for i = 1: length(trial)
    if ~isnan(trial(i).W)
        W_trace_dF(it,:) = trial(i).traceSmooth_dF(n,:);
        W_licks{it}    = trial(i).licks - trial(i).W(1);
        it =1+it;
    end
end

it = 1;
for i = 1: length(trial)
    if ~isnan(trial(i).S)
        S_Taste_dF(it,:) = trial(i).Taste(n,:);
        it =1+it;
    end
end

% n =12;
it = 1;
for i = 1: length(trial)
    if ~isnan(trial(i).N)
        M_Taste_dF(it,:) = trial(i).Taste(n,:);
        it =1+it;
    end
end

% n =12;
it = 1;
for i = 1: length(trial)
    if ~isnan(trial(i).CA)
        CA_Taste_dF(it,:) = trial(i).Taste(n,:);
        it =1+it;
    end
end

% n =12;
it = 1;
for i = 1: length(trial)
    if ~isnan(trial(i).Q)
        Q_Taste_dF(it,:) = trial(i).Taste(n,:);
        it =1+it;
    end
end

it = 1;
for i = 1: length(trial)
    if ~isnan(trial(i).W)
        W_Taste_dF(it,:) = trial(i).Taste(n,:);
        it =1+it;
    end
end

for i = 1: length(trial)
    trace_dF(i,:) = trial(i).traceSmooth_dF(n,:);
end

neuron(n).trial = trial;
neuron(n).framT = trial(1).framT;
neuron(n).T     = trial(1).Tpro;
neuron(n).trace_dF   = trace_dF;
neuron(n).S_trace_dF = S_trace_dF;
neuron(n).M_trace_dF = M_trace_dF;
neuron(n).CA_trace_dF = CA_trace_dF;
neuron(n).Q_trace_dF = Q_trace_dF;
neuron(n).W_trace_dF = W_trace_dF;


neuron(n).S_Taste_dF = S_Taste_dF;
neuron(n).M_Taste_dF = M_Taste_dF;
neuron(n).CA_Taste_dF = CA_Taste_dF;
neuron(n).Q_Taste_dF = Q_Taste_dF;
neuron(n).W_Taste_dF = W_Taste_dF;

neuron(n).S_licks    = S_licks;
neuron(n).M_licks    = M_licks;
neuron(n).CA_licks    = CA_licks;
neuron(n).Q_licks    = Q_licks;
neuron(n).W_licks    = W_licks;






end
%%

%
%% statistical test

