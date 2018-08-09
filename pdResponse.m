function resp = pdResponse(neuron,trial,test)
%% cue response with signed rank
switch test
    case 1  % signed rank
        resp = pdResponse4(neuron,trial);
    case 2
        %% cue response with rank sum
         resp = pdResponse2(neuron,trial);
        %%
    case 3
        %% stats with bootstrap
        resp = pdResponse3(neuron,trial);
end
