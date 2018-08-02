% summaryData is the variable loaded from load('Summary.mat') which is the
% summary file for each mice trained in the decision rig.

for j = 1:length(summaryData)
    for i = 1:size(summaryData(j).raw_decision,1)
        if isempty(summaryData(j).raw_decision{i,8})
            itiL(i) = 0;
        else
            iti = diff(summaryData(j).raw_decision{i,8});
            iti = iti(find(iti>0.1));
            iti = iti(find(iti<0.5));
            itiL(i) = mean(iti);
            if isnan(itiL(i))
                itiL(i) = 0;
            end
        end
        
        if isempty(summaryData(j).raw_decision{i,9})
            itiR(i) = 0;
        else
            iti = diff(summaryData(j).raw_decision{i,9});
            iti = iti(find(iti>0.1));
            iti = iti(find(iti<0.5));
            itiR(i) = mean(iti);
            if isnan(itiR(i))
                itiR(i) = 0;
            end
        end

    end
        itiLavg(j) = mean(itiL(find(itiL>0)));
        itiRavg(j) = mean(itiR(find(itiR>0)));
        clear itiL itiR iti
        
end
barplot_test(1./itiLavg, 1./itiRavg)
%%
for j = 1:length(summaryData)
    for i = 1:size(summaryData(j).raw_decision,1)
        if isempty(summaryData(j).raw_decision{i,6})
            itiL(i) = 0;
        else
            iti = diff(summaryData(j).raw_decision{i,6});
            iti = iti(find(iti>0.1));
            iti = iti(find(iti<0.5));
            itiC(i) = mean(iti);
            if isnan(itiC(i))
                itiC(i) = 0;
            end
        end
        
      end
        itiCavg(j) = mean(itiC(find(itiC>0)));
        clear itiC iti
        
end
CF = 1./itiCavg;
LF = 1./itiLavg;
RF = 1./itiRavg;
CF = mean(CF');
LF = mean(LF');
RF = mean(RF');
% barplot_test(1./itiLavg, 1./itiRavg)