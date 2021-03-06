function [chiselinear chise coefperflinear coefperf adjustedRSlinear adjustedRS]=bj_linearexpo_fitting(testContrast,allMeanPerf,i,startpoint5,analysisType,plotLinear,startPoint,taskColsLight)
%Written by Xing 20/12/12 to calculate chi-squared error associated with
%linear polynomial and exponential curve fits. Returns normalised value of
%chi squared error (summed across sessions).
%Fourth input arg indicates whether: 
%1:input data are for performance in terms of proportion correct/proportion
%of reports that test contrast is higher, or
%0: a different measure of performance, such as PSE or slope.

if nargin<=5||isempty(plotLinear)
    plotLinear=0;
end
if nargin<=6||isempty(startPoint)
    startPoint=0;
end
% colmapText=[colormap(winter(size(testContrast,2)/2));colormap(cool(size(testContrast,2)/2))];
if size(testContrast,2)>7
    colmapText=colormap(jet(size(testContrast,2)));
else
    colmapText=colormap(jet(14));
end
colmapText=[colmapText(1,:);132/255 22/255 216/255;202/255 65/255 223/255;colmapText(3:7,:);157/255 212/255 61/255;colmapText(10:12,:);178/255 111/255 12/255;colmapText(end,:)];
if i==0
    j=1;
    colmapText=[0 0 0];
else
    j=i;
end
lineText='-';
xTemp=1:size(allMeanPerf,1);xTemp=xTemp';
if strcmp(analysisType,'ROC')
    xplot=(1:size(allMeanPerf,1))';
elseif strcmp(analysisType,'NVP')
    xplot=testContrast;
end
yplot=allMeanPerf;
%linear fit:
seedvalslin=[1 1];
if startpoint5==0
    f = fittype('a+b*x');
elseif startpoint5==1
    f = fittype('50+b*x');
    seedvalslin=1;
end
fit3 = fit(xTemp,allMeanPerf,f,'StartPoint',seedvalslin,'Robust','on');
% plot(fit3,'b:');
[clinear,goflinear] = fit(xplot,yplot,fit3)    
coefperflinear=coeffvalues(clinear);
if startpoint5==1
    coefperflinear=[50 coefperflinear];
end
if plotLinear==1
    xcurvevals=1:0.01:size(allMeanPerf,1);
    ycurvevals=coefperflinear(1)+coefperflinear(2)*xcurvevals;
    if ~isempty(taskColsLight)
        plot(xcurvevals+startPoint,ycurvevals,'Color',taskColsLight);
    else
        plot(xcurvevals+startPoint,ycurvevals,'Color',colmapText(j,:));
    end
end
expected=coefperflinear(1)+allMeanPerf*coefperflinear(2);
observed=allMeanPerf;
chiselinear=sum((expected-observed).^2);
chiselinear=chiselinear./(size(allMeanPerf,1)-size(coefperflinear,2)-1);%normalised value of chi squared error
sst=sum((observed-mean(allMeanPerf)).^2);
adjustedRSlinear=1-(chiselinear*(size(allMeanPerf,1)-1)/(sst*(size(allMeanPerf,1)-size(coefperflinear,2))));%adjusted R-square value
sselinear=getfield(goflinear,'sse');
rmselinear=getfield(goflinear,'rmse');
%exponential fit:
if strcmp(analysisType,'ROC')
    if i>7||i==0%if test contrast is higher than sample, or if values to be plotted are averaged across conditions
        seedvals=[-0.1 0.2 1];
    elseif i<=7
        seedvals=[0.1 0 0.01];
    end
elseif strcmp(analysisType,'NVP')
    seedvals=[-0.1 0 0.01];
end
if startpoint5==0
    s= fitoptions('Method','NonlinearLeastSquares','Lower',[-Inf,0,-Inf],...
        'Upper',[Inf,1,Inf],...
        'Startpoint',seedvals);
    f2 = fittype('a*exp(-x)^b+c','options',s);
elseif startpoint5==1
    seedvals=seedvals(1:2);
    s= fitoptions('Method','NonlinearLeastSquares','Lower',[-Inf,0],...
        'Upper',[Inf,1],...
        'Startpoint',seedvals);
    f2 = fittype('50-a+a*exp(-x)^b','options',s);
end
[c2,gof2] = fit(xplot,yplot,f2)
coefperf=coeffvalues(c2);
if startpoint5==1
    coefperf(3)=0;
end
expected=coefperf(1).*exp(-1.*allMeanPerf).^coefperf(2)+coefperf(3);
observed=allMeanPerf;
chise=sum((expected-observed).^2);
chise=chise./(size(allMeanPerf,1)-size(coefperf,2)-1);%normalised value of chi squared error
sst=sum((observed-mean(allMeanPerf)).^2);
adjustedRS=1-(chise*(size(allMeanPerf,1)-1)/(sst*(size(allMeanPerf,1)-size(coefperf,2))));%adjusted R-square value
sse=getfield(gof2,'sse');
rmse=getfield(gof2,'rmse');
if strcmp(analysisType,'ROC')
    xcurvevals=1:0.01:size(allMeanPerf,1);
elseif strcmp(analysisType,'NVP')
    xcurvevals=1:0.01:testContrast(end);
end
if startpoint5==0
    ycurvevals=coefperf(1).*exp(-1.*xcurvevals).^coefperf(2)+coefperf(3);
elseif startpoint5==1
    ycurvevals=50-coefperf(1)+coefperf(1).*exp(-1.*xcurvevals).^coefperf(2)+coefperf(3);
end
if plotLinear==0
    if nargin<=7||~isempty(taskColsLight)
        plot(xcurvevals+startPoint,ycurvevals,'Color',taskColsLight,'LineStyle',lineText);
    else
        plot(xcurvevals,ycurvevals,'Color',colmapText(j,:),'LineStyle',lineText);
    end
end
hold on%,'LineWidth',1%exponential fit