function plot_neurometric_coefs(animal,area,chNum,appendText,startEndTime,slopeNeuro,c50,plotDiffC50_30,diffc50,minRate,maxRate,sessionSorted1,analysisTypeText,example_ch_54,excludeSessHighSSE,excludeOutliers,writeCoefs,slSigmaMultiple,c50SigmaMultiple)
% Handles plotting of correlation coefficients
if plotDiffC50_30==1
    fig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
    set(fig,'PaperUnits','centimeters','PaperType','A4','PaperOrientation', 'portrait', 'PaperPosition', [0.63452 0.63452 21 28.41]);
else
    fighandle2=figure('Color',[1,1,1],'Units','Normalized','Position',[0.14, 0.46, 0.8, 0.4]);
    set(fighandle2, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
    orient landscape
end
excludeSessions=[26 50 306 312 316 322:328 342];
loadText=['load F:\PL\psycho_data\',animal,'\',area,'_psycho_all_sessions.mat psychoAll'];
eval(loadText)
for excludeCount=1:length(excludeSessions)
    ind=(psychoAll(:,1)==excludeSessions(excludeCount));
    slopePsycho=psychoAll(~ind,17);
    sessionSorted2=psychoAll(~ind,1);
end
if plotDiffC50_30==1
    subplot(2,3,1);
else
    subplot(1,3,1);
end
if example_ch_54==1
    xvals=0:20;
    %plot slope vals against ordinal, not actual, session numbers:
    plot(xvals,slopeNeuro,'sb','MarkerFaceColor','b');hold on
    brob = robustfit(xvals,slopeNeuro);
    plot(xvals,brob(1)+brob(2)*xvals,'k','LineWidth',4);hold on
    title('neurometric slope vs time')
    if plotDiffC50_30==1
        subplot(2,3,4);
        plot(xvals,c50,'sb');
        title('c50 vs time')
        subplot(2,2,5);
        plot(xvals,diffc50,'sb');
        title('abs(c50-30) vs time')
    else
        subplot(1,3,2);
        plot(xvals,c50,'sb');
        title('c50 vs time')
    end
    if plotDiffC50_30==1
        subplot(2,2,2);
    else
        subplot(1,3,3);
    end
    matchPsycho=zeros(1,length(slopeNeuro));%find slopePsycho for sessions where data from slopeNeuro is available
    for i=1:length(slopeNeuro)
        ind=find(sessionSorted1(i)==sessionSorted2);
        matchPsycho(i)=slopePsycho(ind(1));%if 2 cells present for 1 channel, simply duplicates value of slopePsycho for that session
    end
    plot(matchPsycho,slopeNeuro,'sb');
    title('neurometric slope vs psychometric slope')
else
    %     plot(sessionSorted1,slopeNeuro,'ok');
    plot(1:length(sessionSorted1),slopeNeuro,'ok');
    title('neurometric slope vs time')
    if plotDiffC50_30==1
        subplot(2,3,4);
        plot(1:length(c50),c50,'ok');
        title('c50 vs time')
        subplot(2,3,5);
        plot(1:length(sessionSorted1),diffc50,'ok');
        title('abs(c50-30) vs time')
    else
        subplot(1,3,2);
        plot(1:length(sessionSorted1),c50,'ok');
        title('c50 vs time')
    end
    %     plot(sessionSorted1,c50,'ok');
    if plotDiffC50_30==1
        subplot(2,3,2);
    else
        subplot(1,3,3);
    end
    matchPsycho=zeros(1,length(slopeNeuro));%find slopePsycho for sessions where data from slopeNeuro is available
    for i=1:length(slopeNeuro)
        ind=find(sessionSorted1(i)==sessionSorted2);
        if ~isempty(ind)
            matchPsycho(i)=slopePsycho(ind(1));%if 2 cells present for 1 channel, simply duplicates value of slopePsycho for that session
        end
    end
    plot(matchPsycho,slopeNeuro,'ok');
    title('neurometric slope vs psychometric slope')
    if plotDiffC50_30==1    
        subplot(2,3,6);
        plot(1:length(sessionSorted1),minRate,'ok');
        if strcmp(analysisTypeText,'ROC')
            title('minimum ROC value vs time')
        else
            title('minimum firing rate vs time')
        end
        subplot(2,3,3);
        plot(1:length(sessionSorted1),maxRate,'ok');
        if strcmp(analysisTypeText,'ROC')
            title('maximum ROC value vs time')
        else
            title('maximum firing rate vs time')
        end
    end
end
% slopeNeuro
% slopePsycho
% matchPsycho
% sessionSorted1
% sessionSorted2
% session1
% session2
a=[sessionSorted1' slopeNeuro'];
c50=real(c50);
b=[sessionSorted1' c50'];
c=[slopeNeuro' matchPsycho'];
[coefficients1 p1]=corrcoef(a);
[coefficients2 p2]=corrcoef(b);
[coefficients3 p3]=corrcoef(c);
if plotDiffC50_30==1
    diffc50=real(diffc50);
    d=[sessionSorted1' diffc50'];
    [coefficients4 p4]=corrcoef(d);
end
e=[sessionSorted1' minRate'];
f=[sessionSorted1' maxRate'];
[coefficients5 p5]=corrcoef(e);
[coefficients6 p6]=corrcoef(f);
if length(coefficients1)>1
    if plotDiffC50_30==1
        coefficients(1,1:6)=[coefficients1(2) coefficients2(2) coefficients3(2) coefficients4(2) coefficients5(2) coefficients6(2)];
        coefficients(2,1:6)=[p1(2) p2(2) p3(2) p4(2) p5(2) p6(2)]
        subplot(2,3,1);
        ptext=sprintf('r= %f  p= %f',coefficients(1,1),coefficients(2,1));
        yLimVals=get(gca,'YLim');
        xLimVals=get(gca,'XLim');
        text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
        subplot(2,3,4);
        ptext=sprintf('r= %f  p= %f',coefficients(1,2),coefficients(2,2));
        yLimVals=get(gca,'YLim');
        xLimVals=get(gca,'XLim');
        text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
        subplot(2,3,2);
        ptext=sprintf('r= %f  p= %f',coefficients(1,3),coefficients(2,3));
        yLimVals=get(gca,'YLim');
        xLimVals=get(gca,'XLim');
        text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
        subplot(2,3,5);
        ptext=sprintf('r= %f  p= %f',coefficients(1,4),coefficients(2,4));
        yLimVals=get(gca,'YLim');
        xLimVals=get(gca,'XLim');
        text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
        subplot(2,3,3);
        ptext=sprintf('r= %f  p= %f',coefficients(1,5),coefficients(2,5));
        yLimVals=get(gca,'YLim');
        xLimVals=get(gca,'XLim');
        text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
        subplot(2,3,6);
        ptext=sprintf('r= %f  p= %f',coefficients(1,6),coefficients(2,6));
        yLimVals=get(gca,'YLim');
        xLimVals=get(gca,'XLim');
        text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
    else
        coefficients(1,1:5)=[coefficients1(2) coefficients2(2) coefficients3(2) coefficients5(2) coefficients6(2)];
        coefficients(2,1:5)=[p1(2) p2(2) p3(2) p5(2) p6(2)]
        subplot(1,3,1);
        ptext=sprintf('r= %f  p= %f',coefficients(1,1),coefficients(2,1));
        yLimVals=get(gca,'YLim');
        xLimVals=get(gca,'XLim');
        text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
        subplot(1,3,2);
        ptext=sprintf('r= %f  p= %f',coefficients(1,2),coefficients(2,2));
        yLimVals=get(gca,'YLim');
        xLimVals=get(gca,'XLim');
        text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
        subplot(1,3,3);
        ptext=sprintf('r= %f  p= %f',coefficients(1,3),coefficients(2,3));
        yLimVals=get(gca,'YLim');
        xLimVals=get(gca,'XLim');
        text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
        subplot(1,3,4);
        ptext=sprintf('r= %f  p= %f',coefficients(1,4),coefficients(2,4));
        yLimVals=get(gca,'YLim');
        xLimVals=get(gca,'XLim');
        text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
    end
    set(gcf,'PaperPositionMode','auto')
    if isempty(chNum)
        chText='mean_across_channels';
    else
        chText=num2str(chNum);
    end
    subFolderName=[analysisTypeText,'_coef_images'];
    if excludeSessHighSSE==0
            coefImagename=[chText,appendText,startEndTime,'_',analysisTypeText,'_coefs_',area];
    elseif excludeSessHighSSE==1
        if excludeOutliers==0
            coefImagename=[chText,appendText,startEndTime,'_',analysisTypeText,'_coefs_',area,'goodSSE'];
        elseif excludeOutliers==1
            coefImagename=[chText,appendText,startEndTime,'_',analysisTypeText,'_coefs_',area,'goodSSE_no_outliers_sl',num2str(slSigmaMultiple),'_C50',num2str(c50SigmaMultiple)];
        end
    end
    coefFolder=fullfile('F:','PL',analysisTypeText,animal,subFolderName);
    if ~exist(coefFolder,'dir')
        mkdir(coefFolder)
    end
    coefPathname=fullfile(coefFolder,coefImagename);
    printtext=sprintf('print -dpng %s.png',coefPathname);
    eval(printtext);
    if writeCoefs==1
        %write correlation coefficients, session numbers, and slope values
        %to file:
        if excludeSessHighSSE==0
            coefMatname=[chText,appendText,startEndTime,'_',analysisTypeText,'_coefs_',area];
        elseif excludeSessHighSSE==1
            if excludeOutliers==0
                coefMatname=[chText,appendText,startEndTime,'_',analysisTypeText,'_coefs_',area,'_goodSSE'];
            elseif excludeOutliers==1
                coefMatname=[chText,appendText,startEndTime,'_',analysisTypeText,'_coefs_',area,'_goodSSE_no_outliers_sl',num2str(slSigmaMultiple),'_C50',num2str(c50SigmaMultiple)];
            end
        end
        subFolderName=[analysisTypeText,'_coef_mat'];
        coefMatFolder=fullfile('F:','PL',analysisTypeText,animal,subFolderName);
        coefMatPathname=fullfile(coefMatFolder,coefMatname);
        if ~exist(coefMatFolder,'dir')
            mkdir(coefMatFolder)
        end
        if plotDiffC50_30==1
            saveText=['save ',coefMatPathname,'.mat coefficients sessionSorted1 slopeNeuro sessionSorted2 slopePsycho matchPsycho c50 diffc50 minRate maxRate'];
        else
            saveText=['save ',coefMatPathname,'.mat coefficients sessionSorted1 slopeNeuro sessionSorted2 slopePsycho matchPsycho c50 minRate maxRate'];
        end
        eval(saveText)
    end
end