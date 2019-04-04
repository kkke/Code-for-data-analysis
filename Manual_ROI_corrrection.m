clear;close all
load('data_CNMF')
[CC,jsf,im] = plot_contours(A_keep,Cn,options,1);
gui2p_Calman3(CC,Cn,jsf,F_dff);
%%
% figure;
% plot(F_dff(104,:))
% plot_components_GUI(data,A_keep,C_keep,b,f,Cn,options);