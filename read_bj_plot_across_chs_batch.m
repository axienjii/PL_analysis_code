function read_bj_plot_across_chs_batch
%Written on 05/03/13.

allChCoefs=cell(1,4);
areas=[{'v4_1'} {'v1_1'}];
animals=[{'blanco'} {'jack'}];
plotFig=0;
totalSigChs=0;
for animalInd=1:2
    for areaInd=1:2
        animal=animals{animalInd};
        area=areas{areaInd};
        if strcmp(animal,'jack')
            if strcmp(area,'v4_1')
                testContrasts=[10 15 20 25 27 28 29 31 32 33 35 40 50 60];
                sampleContrasts=30;
            elseif strcmp(area,'v1_1')
                testContrasts=[5 10 15 20 22 25 28 32 35 40 45 50 60 90];
                sampleContrasts=30;
            elseif strcmp(area,'v1_2')
                testContrasts=[5 10 12 15 18 22 25 28 35 45 63 90;5 10 15 22 25 28 32 35 38 45 60 90;5 10 15 25 32 35 38 42 45 50 60 90];
                sampleContrasts=[20 30 40];
            end
        elseif strcmp(animal,'blanco')
            if strcmp(area,'v4_1')
                testContrasts=[10 15 20 25 27 28 29 31 32 33 35 40 50 60];
                sampleContrasts=30;
            elseif strcmp(area,'v1_1')
                testContrasts=[5 10 15 20 22 25 28 32 35 40 45 50 60 90];
                sampleContrasts=30;
            elseif strcmp(area,'v1_2')
                testContrasts=[5 10 12 15 18 22 25 28 35 45 63 90;5 10 15 22 25 28 32 35 38 45 60 90;5 10 15 25 32 35 38 42 45 50 60 90];
                sampleContrasts=[20 30 40];
            end
        end        
        if plotFig==1
            fighandle1=  figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
            set(fighandle1, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
            fighandle2=  figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
            set(fighandle2, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
            fighandle3=  figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
            set(fighandle3, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
            fighandle4=  figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
            set(fighandle4, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
        end
        channels = main_channels(animal,area);
        for i=1:length(channels)
            chNum=channels(i);
            for j=1:length(sampleContrasts)
                sampleContrast=sampleContrasts(j);
                testContrast=testContrasts(j,:);
                coefsFileName=[num2str(chNum),'_',num2str(sampleContrast),'_roc_coefs_',area,'_goodSSE_no_outliers_sl3_C503','.mat'];
                coefsPathName=fullfile('F:','PL','ROC_coef_mat',animal,coefsFileName);
                loadText=['load ',coefsPathName];
                eval(loadText);
                if find(coefficients(2,:)<0.05)
                    totalSigChs=totalSigChs+1
                end
                figureHandles=[{'fighandle1'} {'fighandle2'} {'fighandle3'} {'fighandle4'}];
                xVals=[{1:length(slopeNeuro)} {1:length(c50)} {matchPsycho} {1:length(diffc50)}];
                yVals=[{slopeNeuro} {c50} {slopeNeuro} {diffc50}];
                for k=1:size(coefficients,2)
                    if plotFig==1
                        figText=sprintf('figure(%s)',figureHandles{k});
                        eval(figText)
                        subplot(ceil(length(channels)/5),5,i);
                        plot(xVals{k},yVals{k},'ok');
                        yLimVals=get(gca,'YLim');
                        xLimVals=get(gca,'XLim');
                        if k==1||k==4
                            if xLimVals(1)>=0
                                xlim([0 xLimVals(2)]);
                            end
                            if yLimVals(1)>=0
                                ylim([0 yLimVals(2)]);
                            end
                        end
                        ptext=sprintf('R: %.3f p: %.4f',coefficients(1,k),coefficients(2,k));
                        if coefficients(2,k)<0.05
                            colText='g';
                        else
                            colText='k';
                        end
                        yLimVals=get(gca,'YLim');
                        xLimVals=get(gca,'XLim');
                        title(num2str(chNum));
                        text('Position',[xLimVals(1)+(xLimVals(2)-xLimVals(1))/5 yLimVals(1)+(yLimVals(2)-yLimVals(1))/6],'FontSize',9,'String',ptext,'Color',colText);
                    end
                    allChCoefs{k}=[allChCoefs{k};{animal},{area},chNum,coefficients(1,k),coefficients(2,k)];
                end
            end
        end
    end
end
comparison=[{'slope_vs_session'} {'C50_vs_session'} {'neuro_vs_psycho'} {'diffC50_vs_session'}];
for k=1:size(coefficients,2)
    figText=sprintf('figure(%s)',figureHandles{k});
    eval(figText)
    orient landscape
    rocImagename=[area,'_allChROC_',comparison{k}];
    allChROCPathname=fullfile('F:','PL','neurometric_ROC',animal,rocImagename);
    printtext=sprintf('print -dpng %s.png',allChROCPathname);
    eval(printtext);
    printtext=sprintf('print -deps %s.eps',allChROCPathname);
    eval(printtext);
end

%combine ROC values across channels:
loadText=['load ',folder,'\allChROC.mat allChROC'];
eval(loadText);
sessions=unique(allChROC(:,2));
fighandle2=figure('Color',[1,1,1],'Units','Normalized','Position',[0.14, 0.46, 0.8, 0.4]);
set(fighandle2, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
orient landscape
allChSlopes=[];%neurometric and session number
allChC50s=[];
allChSlopesNP=[];%neurometric and psychometric slopes
for i=1:length(sessions)
    session=sessions(i);
    rowInd=find(allChROC(:,2)==session);
    chs=allChROC(rowInd,1);
    if chs~=unique(chs)
        sprintf('duplicate channel data for session %s!',num2str(session))
        pause
    end
    allChvals=allChROC(rowInd,3:end);%compile ROC values, neurometric slopes, and C50s, across channels for each session
%     if length(chs)==length(channels)
%         allChSlopes=[allChSlopes allChvals(:,end-1)];
%         allChC50s=[allChC50s allChvals(:,end)];
%     else
%         padding=NaN(length(channels)-length(chs),1);
%         slopePadded=[allChvals(:,end-1);padding];
%         allChSlopes=[allChSlopes slopePadded];
%         c50sPadded=[allChvals(:,end);padding];
%         allChC50s=[allChC50s c50sPadded];
%     end
    for j=1:length(chs)
        allChSlopes=[allChSlopes;session allChvals(j,end-1)];
        allChC50s=[real(allChC50s);session allChvals(j,end)];
        allChSlopesNP=[allChSlopesNP;allChvals(j,end-1) matchPsycho(i)];
    end
    avSlopes(i)=mean(allChvals(:,end-1));
    avC50s(i)=real(mean(allChvals(:,end)));
    stdSlopes(i)=std(allChvals(:,end-1));
    stdC50s(i)=std(allChvals(:,end));
    test=subplot(1,3,1);
    plot(session,allChvals(:,end-1),'ok');hold on
    title('neurometric slope vs time')
    test=subplot(1,3,2);
    plot(session,allChvals(:,end),'ok');hold on
    title('c50 vs time')
    test=subplot(1,3,3);
    plot(matchPsycho(i),allChvals(:,end-1),'ok');hold on
    title('neurometric slope vs psychometric slope')
end
[RHO,PVAL]=corr(allChSlopes)
test=subplot(1,3,1);
ptext=sprintf('r= %f  p= %f',RHO(1,2),PVAL(1,2));
yLimVals=get(gca,'YLim');
xLimVals=get(gca,'XLim');
text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
[RHO,PVAL]=corr(allChC50s)
test=subplot(1,3,2);
ptext=sprintf('r= %f  p= %f',RHO(1,2),PVAL(1,2));
yLimVals=get(gca,'YLim');
xLimVals=get(gca,'XLim');
text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
[RHO,PVAL]=corr(allChSlopesNP)
test=subplot(1,3,3);
ptext=sprintf('r= %f  p= %f',RHO(1,2),PVAL(1,2));
yLimVals=get(gca,'YLim');
xLimVals=get(gca,'XLim');
text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
%V1 non-roving: slope correlation R=0.4670, p=0; C50 correlation R=0.1020,
%p=0.0177 neuro versus psycho slope correlation R=0.4045, p=0.
coefFilename=['all_roc_neuro_psycho_',area];
coefPathname=fullfile('F:','PL','ROC_coef_images',animal,coefFilename);
set(gcf,'PaperPositionMode','auto')
printtext=sprintf('print -dpng %s',coefPathname);
eval(printtext);
% fighandle2=figure('Color',[1,1,1],'Units','Normalized','Position',[0.14, 0.46, 0.8, 0.4]);
% set(fighandle2, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
% orient landscape
% test=subplot(1,3,1);
% plot(sessions,avSlopes,'ko');hold on
% errorbar(sessions,avSlopes,stdSlopes);hold on
% title('mean neurometric slope across channels vs time')
% a=[avSlopes' sessions];
% b=[avC50s' sessions];
% [coefficients1 p1]=corrcoef(a); 
% [coefficients2 p2]=corrcoef(b); 
% coefficients(1,1:2)=[coefficients1(2) coefficients2(2)];
% coefficients(2,1:2)=[p1(2) p2(2)]
% ptext=sprintf('r= %f  p= %f',coefficients(1,1),coefficients(2,1));
% yLimVals=get(gca,'YLim');
% xLimVals=get(gca,'XLim');
% text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
% test=subplot(1,3,2);
% plot(sessions,avC50s,'ko');hold on
% errorbar(sessions,avC50s,stdSlopes);hold on
% title('mean C50 across channels vs time')
% ptext=sprintf('r= %f  p= %f',coefficients(1,2),coefficients(2,2));
% yLimVals=get(gca,'YLim');
% xLimVals=get(gca,'XLim');
% text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
% coeffname=[folder,'/mean_roc_neuro_psycho'];
% set(gcf,'PaperPositionMode','auto')
% printtext=sprintf('print -dpng %s',coeffname);
% eval(printtext);
    

