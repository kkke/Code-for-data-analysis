%% this code is to reorganize the trial structure into neuron structure

function neuron = trial2neuron(trial);
for n = 1:size(trial(1).trace,1)
    for i = 1: length(trial)
        Licktrace_dF(i,:) = trial(i).LickTrace(n,:);
    end
    
    for i = 1: length(trial)
        trace_dF(i,:) = trial(i).traceSmooth_dF(n,:);
    end

neuron(n).trial = trial;
neuron(n).framT = trial(1).framT;
neuron(n).TLick = trial(1).TLick;
neuron(n).trace_dF     = trace_dF;
neuron(n).Licktrace_dF = Licktrace_dF;
end
%%

%
%% statistical test

