function plot_example_DICAF_channel
%load example ch 6 from Jack V4, show changes in raw DICAF values with
%training (colour coded by condition).
%created for SF March 2014 ppt presentation.
load('F:\PL\ROC\jack\v4_1\ROC_Ch6_30_1024_to_1536.mat')
testContrast=[10 15 20 25 27 28 29 31 32 33 35 40 50 60];
sessions=[];
colmapText=colormap(jet(size(testContrast,2)));
colmapText=[colmapText(1,:);132/255 22/255 216/255;202/255 65/255 223/255;colmapText(3:7,:);157/255 212/255 61/255;colmapText(10:12,:);178/255 111/255 12/255;colmapText(end,:)];
for sessNum=1:size(ROCmat,1)
    sessions=[sessions ROCmat{sessNum,1}];
end
[sessionSorted sortInd]=sort(sessions);
ROCmatSorted=[];
ROCvals=[];
for sessNum=1:size(ROCmat,1)
    ROCmatSorted=[ROCmatSorted;ROCmat(sortInd(sessNum),:)];
    ROCvals=[ROCvals;ROCmat{sessNum,3}];
end
markerText='x';
markerS=8;
figure
subplot(2,2,1)
xlim([0 size(ROCvals,1)+1]);
for i=1:size(ROCvals,2)
    plot(1:size(ROCvals,1),ROCvals(:,i),'Color',colmapText(i,:),'LineStyle','-','Marker',markerText,'MarkerFaceColor',colmapText(i,:),'MarkerEdgeColor',colmapText(i,:),'MarkerSize',markerS);hold on
end
xlim([0 size(ROCvals,1)+1]);
for i=1:size(ROCvals,2)
    yLimVals=get(gca,'ylim');
    xLimVals=get(gca,'xlim');
    unitSpace=(yLimVals(2)-yLimVals(1))/30;
    text('Position',[xLimVals(2)+(xLimVals(2)-xLimVals(1))/25 yLimVals(1)+unitSpace*i*2],'FontSize',9,'String',[markerText,'  ',num2str(testContrast(i)),'%'],'Color',colmapText(i,:));
end
xlabel('session number');
ylabel('DICAF');

figure
copperCols1=[];
copperCols2=[];
for colMapInd=1:ceil(size(ROCmat,1)/2)
    copperCols1(colMapInd,:)=[1 0 (colMapInd-1)/ceil((size(ROCmat,1))/2)];
end
for colMapInd=1:size(ROCmat,1)-size(ROCmat,1)/2
    copperCols2(colMapInd,:)=[1-colMapInd/floor((size(ROCmat,1))/2) 0 1];
end
copperCols=[copperCols1;copperCols2];
subplotInd=0;
for i=[1 9 16 25]
    subplotInd=subplotInd+1;
    subplot(1,4,subplotInd);
    for j=1:length(testContrast)
        plot(testContrast(j),ROCvals(i,j),'Color',colmapText(j,:),'LineStyle','-','Marker',markerText,'MarkerFaceColor',colmapText(j,:),'MarkerEdgeColor',colmapText(j,:),'MarkerSize',markerS);hold on
    end
    datavals=ROCvals(i,:);
    if sum(datavals(1:3))<=sum(datavals(end-2:end))
        X0=[2 30 0.2 0.1];
    elseif sum(datavals(1:3))>sum(datavals(end-2:end))
        X0=[-2 30 0.2 0.1];
    end
    options = optimset('Display','off','MaxFunEvals',10^4,'MaxIter',10^4,'TolFun',1.0E-6,'TolX',1.0E-6);
    X1=fminsearch(@fit_weibull,X0,options,testContrast,datavals,[],'least_square',[1 1 1 0],[10 100 1 0],[1 1 0 0],[-20 0 0 0]);
    %         X=fminsearch(@fit_weibull,X1,options,testContrast,datavals,[],'mle',[1 1 1 1],[2 100 1 0.2],[1 1 0 0],[-20 0 0 0]);
    X=fminsearch(@fit_weibull,X1,options,testContrast,datavals,[],'mle',[1 1 1 0],[10 100 1 1],[1 1 0 0],[-10 0 0 0]);
    fitted_yvals=1-X(4)-X(3).*exp(-(testContrast./X(2)).^X(1));
    residuals=datavals-fitted_yvals;
    sseCRF=sum(residuals.^2);
%     chSSE(chInd,1:2)=[chInd sseCRF];
%     if sseCRF>0.1%if the fit seems poor, try a variety of values for the upper and lower limits
%         coefEsts2 = zeros(6,4);
%         InitVar=1;
%         for upperMax=[X1(3) X1(3)+0.1]
%             for lowerMin=[X1(4) X1(4)+0.1]
%                 [coefEsts2(InitVar,:)]=fminsearch(@fit_weibull,X1,options,testContrast,datavals,[],'mle',[1 1 1 0],[20 100 upperMax 0],[1 1 0 1],[-20 0 0 lowerMin]);
%                 fitted_yvals=1-coefEsts2(InitVar,4)-coefEsts2(InitVar,3).*exp(-(testContrast./coefEsts2(InitVar,2)).^coefEsts2(InitVar,1));
%                 residuals=datavals-fitted_yvals;
%                 sseCRFtemp(InitVar)=sum(residuals.^2);
%                 InitVar=InitVar+1;
%             end
%         end
%         [minSSE I]=min(sseCRFtemp);
%         if minSSE<sseCRF
%             X=coefEsts2(I,:);
%             chSSE(chInd,1:2)=[chInd minSSE];
%         end
%     end
    hold on
    xvals=testContrast(1):1:testContrast(end);
    yvals=1-X(4)-X(3).*exp(-(xvals./X(2)).^X(1));
%     plot(xvals,yvals,'Color',copperCols(i,:));
    plot(xvals,yvals,'Color','k');
    title(num2str(i));
    ylim([0 1]);
    xlim([5 70]);
end
