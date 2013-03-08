function plot_neurometric_coefs(animal,area,chNum,appendText,startEndTime,slopeNeuro,c50,plotDiffC50_30,diffc50,minRate,maxRate,sessionSorted1,analysisTypeText,example_ch_54,excludeSessHighSSE,excludeOutliers,writeCoefs)

excludeSessions=[26 50 306 312 316 322:328 342];
loadText=['load F:\PL\psycho_data\',animal,'\',area,'_psycho_all_sessions.mat psychoAll'];
eval(loadText)
for excludeCount=1:length(excludeSessions)
    ind=(psychoAll(:,1)==excludeSessions(excludeCount));
    slopePsycho=psychoAll(~ind,17);
    sessionSorted2=psychoAll(~ind,1);
end
if plotDiffC50_30==1
    subplot(2,2,1);
else
    subplot(1,3,1);
end
if example_ch_54==1
    xvals=0:20;
    plot(xvals,slopeNeuro,'sb','MarkerFaceColor','b');hold on
    brob = robustfit(xvals,slopeNeuro);
    plot(xvals,brob(1)+brob(2)*xvals,'k','LineWidth',4);hold on
    title('neurometric slope vs time')
    if plotDiffC50_30==1
        subplot(1,3,3);
        plot(xvals,c50,'sb');
        title('c50 vs time')
        subplot(2,2,4);
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
        subplot(2,2,3);
        plot(1:length(c50),c50,'ok');
        title('c50 vs time')
        subplot(2,2,4);
        plot(1:length(sessionSorted1),diffc50,'ok');
        title('abs(c50-30) vs time')
    else
        subplot(1,3,2);
        plot(1:length(sessionSorted1),c50,'ok');
        title('c50 vs time')
    end
    %     plot(sessionSorted1,c50,'ok');
    if plotDiffC50_30==1
        subplot(2,2,2);
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
        subplot(2,2,1);
        ptext=sprintf('r= %f  p= %f',coefficients(1,1),coefficients(2,1));
        yLimVals=get(gca,'YLim');
        xLimVals=get(gca,'XLim');
        text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
        subplot(2,2,3);
        ptext=sprintf('r= %f  p= %f',coefficients(1,2),coefficients(2,2));
        yLimVals=get(gca,'YLim');
        xLimVals=get(gca,'XLim');
        text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
        subplot(2,2,2);
        ptext=sprintf('r= %f  p= %f',coefficients(1,3),coefficients(2,3));
        yLimVals=get(gca,'YLim');
        xLimVals=get(gca,'XLim');
        text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
        subplot(2,2,4);
        ptext=sprintf('r= %f  p= %f',coefficients(1,4),coefficients(2,4));
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
    end
    set(gcf,'PaperPositionMode','auto')
    if isempty(chNum)
        chText='mean_across_channels';
    else
        chText=num2str(chNum);
    end
    subFolderName=[analysisTypeText,'_coef_images'];
    if excludeSessHighSSE==0
            crfCoefImagename=[chText,appendText,startEndTime,'_',analysisTypeText,'_coefs'];
    elseif excludeSessHighSSE==1
        if excludeOutliers==0
            crfCoefImagename=[chText,appendText,startEndTime,'_,analysisTypeText,_coefs_goodSSE'];
        elseif excludeOutliers==1
            crfCoefImagename=[chText,appendText,startEndTime,'_,analysisTypeText,_coefs_goodSSE_no_outliers_sl',num2str(slSigmaMultiple),'_C50',num2str(c50SigmaMultiple)];
        end
    end
    crfCoefFolder=fullfile('F:','PL',analysisTypeText,animal,subFolderName);
    if ~exist(crfCoefFolder,'dir')
        mkdir(crfCoefFolder)
    end
    crfCoefPathname=fullfile(crfCoefFolder,crfCoefImagename);
    printtext=sprintf('print -dpng %s.png',crfCoefPathname);
    eval(printtext);
    if writeCoefs==1
        %write correlation coefficients, session numbers, and slope values
        %to file:
        if excludeSessHighSSE==0
            crfCoefMatname=[chText,appendText,startEndTime,'_',analysisTypeText,'_coefs_',area];
        elseif excludeSessHighSSE==1
            if excludeOutliers==0
                crfCoefMatname=[chText,appendText,startEndTime,'_',analysisTypeText,'_coefs_',area,'_goodSSE'];
            elseif excludeOutliers==1
                crfCoefMatname=[chText,appendText,startEndTime,'_',analysisTypeText,'_coefs_',area,'_goodSSE_no_outliers_sl',num2str(slSigmaMultiple),'_C50',num2str(c50SigmaMultiple)];
            end
        end
        subFolderName=[analysisTypeText,'_coef_mat'];
        crfCoefMatFolder=fullfile('F:','PL',analysisTypeText,animal,subFolderName);
        crfCoefMatPathname=fullfile(crfCoefMatFolder,crfCoefMatname);
        if ~exist(crfCoefMatFolder,'dir')
            mkdir(crfCoefMatFolder)
        end
    if plotDiffC50_30==1
        saveText=['save ',crfCoefMatPathname,'.mat coefficients sessionSorted1 slopeNeuro sessionSorted2 slopePsycho matchPsycho c50 diffc50 minRate maxRate'];
    else
        saveText=['save ',crfCoefMatPathname,'.mat coefficients sessionSorted1 slopeNeuro sessionSorted2 slopePsycho matchPsycho c50 minRate maxRate'];
    end
    eval(saveText)
    end
end