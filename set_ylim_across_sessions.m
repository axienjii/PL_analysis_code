function set_ylim_across_sessions(chNum,numsessions,sampleContrast,testContrast,c50,chSSE,yLimCRF,fileText)
if ischar(chNum)
    titleText=chNum;
else
    titleText=num2str(chNum);
end
%standardises y-limits across all subplots
for i=1:numsessions
    subplot(ceil(numsessions/5),5,i);
    ylim([min(yLimCRF(:,1)) max(yLimCRF(:,2))]);
    line(c50(1,i),0:max(yLimCRF(:,2))/100:max(yLimCRF(:,2)),'Color','r');
    line(sampleContrast,0:max(yLimCRF(:,2))/100:max(yLimCRF(:,2)));
    if i==1
        crfname=[titleText,' ',num2str(sampleContrast),' ',fileText];
        ptext=sprintf('%s',crfname);
        orient landscape
        yLimVals=get(gca,'YLim');
        text('Position',[-10 yLimVals(2)+0.2*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
    end
    ptext=sprintf('SSE: %.3f',chSSE(i,2));
    text('Position',[testContrast(end)/2 min(yLimCRF(:,1))+(max(yLimCRF(:,2))-min(yLimCRF(:,1)))/4],'FontSize',9,'String',ptext);
end