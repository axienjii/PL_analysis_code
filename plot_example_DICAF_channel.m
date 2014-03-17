function plot_example_DICAF_channel
%load example ch 6 from Jack V4, show changes in raw DICAF values with
%training (colour coded by condition).
load('F:\PL\ROC\jack\v4_1\ROC_Ch6_30_1024_to_1536.mat')
testContrast=[10 15 20 25 27 28 29 31 32 33 35 40 50 60];
sessions=[];
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