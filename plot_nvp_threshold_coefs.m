function plot_nvp_threshold_coefs(animal,area,chNum,appendText,startEndTime,threshold82lower,threshold82higher,sessionSorted1,analysisTypeText,excludeSessHighSSE,excludeOutliers,writeCoefs,threshSigmaMultiple)
% Handles plotting of correlation coefficients
tLNeuro=threshold82lower;%lower neurometric threshold
tHNeuro=threshold82higher;
fighandle2=figure('Color',[1,1,1],'Units','Normalized','Position',[0.14, 0.46, 0.8, 0.4]);
set(fighandle2, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
orient landscape
psychoThresholdMatName=['F:\PL\psycho_data\',animal,'\',area,'_psyThreshold.mat'];
loadText=['load ',psychoThresholdMatName,' sessionSorted2 threshold82lower threshold82higher'];
eval(loadText)
tLPsycho=threshold82lower;%lower psychometric threshold
tHPsycho=threshold82higher;
% for excludeCount=1:length(excludeSessions)
%     ind=(psychoAll(:,1)==excludeSessions(excludeCount));
%     slopePsycho=psychoAll(~ind,17);
%     sessionSorted2=psychoAll(~ind,1);
% end
subplot(2,2,1);
plot(1:length(sessionSorted1),tLNeuro,'ok');hold on
title('lower contrast threshold vs time')%neurometric: black; psychometric: red
subplot(2,2,2);
plot(1:length(sessionSorted1),tHNeuro,'ok');hold on
title('higher contrast threshold vs time')
subplot(2,2,3);
matchPsychoLower=zeros(1,length(tLNeuro));%find slopePsycho for sessions where data from slopeNeuro is available
for i=1:length(tLNeuro)
    ind=find(sessionSorted1(i)==sessionSorted2);
    if ~isempty(ind)
        matchPsychoLower(i)=tLPsycho(ind(1));%if 2 cells present for 1 channel, simply duplicates value of slopePsycho for that session
    end
end
plot(matchPsychoLower,tLNeuro,'ok');
title('lower contrast neurometric vs psychometric threshold')
subplot(2,2,4);
matchPsychoHigher=zeros(1,length(tHNeuro));%find slopePsycho for sessions where data from slopeNeuro is available
for i=1:length(tHNeuro)
    ind=find(sessionSorted1(i)==sessionSorted2);
    if ~isempty(ind)
        matchPsychoHigher(i)=tHPsycho(ind(1));%if 2 cells present for 1 channel, simply duplicates value of slopePsycho for that session
    end
end
plot(matchPsychoHigher,tHNeuro,'ok');
title('lower contrast neurometric vs psychometric threshold')
subplot(2,2,1);
plot(1:length(sessionSorted1),matchPsychoLower,'or');
subplot(2,2,2);
plot(1:length(sessionSorted1),matchPsychoHigher,'or');

a=[sessionSorted1' tLNeuro'];
b=[sessionSorted1' tHNeuro'];
c=[sessionSorted1' matchPsychoLower'];
d=[sessionSorted1' matchPsychoHigher'];
e=[threshold82lower' matchPsychoLower'];
f=[threshold82lower' matchPsychoHigher'];
[coefficients1 p1]=corrcoef(a);
[coefficients2 p2]=corrcoef(b);
[coefficients3 p3]=corrcoef(c);
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
        saveText=['save ',coefMatPathname,'.mat coefficients sessionSorted1 tLNeuro tLPsycho sessionSorted2 tHNeuro tHPsycho matchPsychoLower matchPsychoHigher'];
        eval(saveText)
    end
end