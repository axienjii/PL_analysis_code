function [r p]=bj_plot_coef_expo_perf(coefperf,sampleContrast,testContrast,ind)
%Written on 08/01/13 by Xing.
%Plots coefficient related to slope of exponential fitted curve, against
%absolute value of difference between test and sample contrasts.

xlowertestdiff=fliplr(sampleContrast-testContrast(1:ind));
ylowertestdiff=flipud(coefperf(1:ind,1));
plot(xlowertestdiff,ylowertestdiff,'ok','LineStyle',':');hold on%test of lower contrast than sample
plot(testContrast(ind+1:end)-sampleContrast,coefperf(ind+1:end,1)*-1,'ok','MarkerFaceColor','k','LineStyle','-');hold on
[r(1),p(1)] = corr(xlowertestdiff',ylowertestdiff,'type','Spearman');
[r(2),p(2)] = corr(testContrast(ind+1:end)'-sampleContrast,coefperf(ind+1:end,1)*-1,'type','Spearman');
