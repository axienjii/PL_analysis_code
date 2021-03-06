function set_ylim_across_sessions(chNum,numsessions,sampleContrast,testContrast,c50,chSSE,yLimData,fileText)
if ischar(chNum)
    titleText=chNum;
else
    titleText=num2str(chNum);
end
%standardises y-limits across all subplots
for i=1:numsessions
    subplot(ceil(numsessions/5),5,i);
    ylim([min(yLimData(:,1)) max(yLimData(:,2))]);
    if strcmp(fileText,'ROC_diff2')
        line(c50(1,i),0:max(yLimData(:,2))/100:max(yLimData(:,2)),'Color','b');
    else
        line(c50(1,i),0:max(yLimData(:,2))/100:max(yLimData(:,2)),'Color','r');
    end
    line(sampleContrast,0:max(yLimData(:,2))/100:max(yLimData(:,2)));
    if i==1
        if strcmp(fileText,'ROC_diff2')
        dataName=[titleText,' ',num2str(sampleContrast),' ROCdiff'];
        else
        dataName=[titleText,' ',num2str(sampleContrast),' ',fileText];
        end
        ptext=sprintf('%s',dataName);
        orient landscape
        yLimVals=get(gca,'YLim');
        text('Position',[-10 yLimVals(2)+0.2*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
    end
    ptext=sprintf('SSE: %.3f',chSSE(i,2));
    if ~strcmp(fileText,'ROC_diff2')&&~strcmp(fileText,'ROC_diff')
        text('Position',[testContrast(end)/2 min(yLimData(:,1))+(max(yLimData(:,2))-min(yLimData(:,1)))/4],'FontSize',9,'String',ptext);
    end
end