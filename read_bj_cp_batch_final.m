function read_bj_cp_batch_final(animal,area,roving)
%Written by Xing 08/06/14
%Modified from read_bj_cp_batch, plots mean CPs with multi-day running average (specified by numSmoothSessions),
%calculates statistics for early vs late sessions.
%to calculate and write choice probability correlation coefficients
useRobustfit=1;
onExternalHD=0;
numSmoothSessions=5;%number of sessions over which to calculate moving average
if onExternalHD==1
    rootFolder='G:\PL_backup_060413';
else
    rootFolder='F:';
end
if nargin<7||isempty(plotLeastSquares)
    plotLeastSquares=[];
end
if nargin<1||isempty(animal)
    animals=[{'blanco'} {'jack'}];
end
if nargin<2||isempty(area)
    if roving==0
        areas=[{'v4_1'} {'v1_1'}];
    elseif roving==1
        areas=[{'v1_2_1'} {'v1_2_2'}];
        areas=[{'v1_2_1'}];
    end
end
sigChsRs=[];
statsObserved=cell(length(animals),length(areas));
analysisType='CP';
MarkerType='.';
markerTexts='+x';
markerText=MarkerType;
excludeSessions=[26 50 306 312 316 322:328 342 398 451];
excludeSessions=[26 50 342 398 451];
test_epochs={0 512 512*2 512*3};durSpon=150;
if roving==0
    figCPconds=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
    set(figCPconds, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
    figAve=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.5, 0.8]); %
    set(figAve, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
    figMovingAve=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.25, 0.4]); %
    set(figMovingAve, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
    ttestTable=cell(length(animals),length(areas));
    ttestTableChs=cell(length(animals),length(areas));
    ttestTableChRs=cell(length(animals),length(areas));
elseif roving==1
    for animalInd=1:2
        for areaInd=1:2
            figCPconds(animalInd+(areaInd-1)*2)=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
            set(figCPconds(animalInd+(areaInd-1)*2), 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
        end
    end
    ttestTable=cell(length(animals),3);
    ttestTableChs=cell(length(animals),3);
    ttestTableChRs=cell(length(animals),3);
end
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        sessionNums = main_raw_sessions_final(animal,area,[],0);
        channels=main_channels(animal,area);
        if strcmp(animal,'blanco')&&strcmp(area,'v4_1')
            channels=channels(channels~=4);
        end
        [sampleContrasts testContrasts]=area_metadata(area);
        colmapText=colormap(jet(size(testContrasts,2)));
        colmapText=[colmapText(1,:);132/255 22/255 216/255;202/255 65/255 223/255;colmapText(3:7,:);157/255 212/255 61/255;colmapText(10:12,:);178/255 111/255 12/255;colmapText(end,:)];
        for sampleContrastsInd=1:length(sampleContrasts)
            spearmanTableCh=[];
            spearmanTableChR=[];
            spearmanTableCh=[];
            spearmanTableChR=[];
            sampleContrast=sampleContrasts(sampleContrastsInd);
            testContrast=testContrasts(sampleContrastsInd,:);
            if strncmp(area,'v1',2)
                closeCond=find(abs(testContrast-sampleContrast)<=5);
            elseif strncmp(area,'v4',2)
                closeCond=find(abs(testContrast-sampleContrast)<5);
            end
            for epoch=4:4
                periods=[test_epochs{epoch-1} test_epochs{epoch}(1)];
                for subPeriod=1:length(periods)-1
                    startEndTime=['_',num2str(periods(subPeriod)),'_to_',num2str(periods(subPeriod+1))];
                    figChCPconds=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
                    set(figChCPconds, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                    allChCP=cell(length(sessionNums),length(closeCond));
                    for i=1:length(channels)
                        matName=[analysisType,'_Ch',num2str(channels(i)),'_',num2str(sampleContrast),startEndTime,'.mat'];
                        matPath=fullfile('F:','PL',analysisType,animal,area,matName);
                        loadText=['load ',matPath,' ',analysisType,'mat'];
                        eval(loadText);
                        sessionList=cell2mat(CPmat(:,1));
                        [sessionList ind]=sort(sessionList);
                        CPmat=CPmat(ind,:);
                        saveText=['save ',matPath,' ',analysisType,'mat'];
                        eval(saveText);
                        dataArray=CPmat;
                        includeMatch=[];
                        for includeInd=1:size(dataArray,1)
                            if sum(excludeSessions==cell2mat(dataArray(includeInd,1)))==0
                                includeMatch=[includeMatch includeInd];
                            end
                        end
                        dataArray=dataArray(includeMatch,:);
                        %figure combining data across channels:
                        if roving==0
                            figure(figCPconds);
                        elseif roving==1
                            figure(figCPconds(animalInd+(areaInd-1)*2));
                        end
                        for cond=1:length(closeCond)
                            chCPs=[];
                            if roving==0
                                subplot(length(closeCond),2*length(areas),1+((cond-1)*2*length(areas))+animalInd-1+2*(areaInd-1));
                                %subplot(2*length(areas),length(closeCond),cond+length(closeCond)*(animalInd+2*(areaInd-1)-1));
                            elseif roving==1
                                subplot(length(sampleContrasts),length(closeCond),cond+length(closeCond)*(sampleContrastsInd-1));
                            end
                            for sessionInd=1:size(dataArray,1)
                                if length(dataArray{sessionInd,3})==length(closeCond)
                                    plot(sessionInd,dataArray{sessionInd,3}(cond),'LineStyle','none','Marker',MarkerType,'Color',colmapText(closeCond(cond),:));hold on
                                    allChCP{sessionInd,cond}=[allChCP{sessionInd,cond} dataArray{sessionInd,3}(cond)];
                                    chCPs=[chCPs dataArray{sessionInd,3}(cond)];
                                end
                            end
                            xlim([0 length(sessionNums)+1]);
                            ylim([0 1]);
                            title([num2str(testContrast(closeCond(cond))),'%']);
                            if roving==0
                                if cond==1
%                                     title(['Monkey ',num2str(animalInd)]);
                                    if animalInd==1&&areaInd==1&&cond==1
                                        xlabel('Session number');
                                        ylabel('Choice probability');
                                    end
                                end
                            elseif roving==1
                                if cond==1&&sampleContrastsInd==1
%                                     title(['Monkey ',num2str(animalInd)]);
                                    xlabel('Session number');
                                    ylabel('Choice probability');
                                end
                            end
                            [rho p]=corr(chCPs',[1:length(sessionNums)]','type','Spearman');
                            spearmanTableCh(i,cond)=p;
                            spearmanTableChR(i,cond)=rho;
                        end
                        if animalInd==2&&cond==1&&i==1
                            for cond=1:length(closeCond)
                                yLimVals=get(gca,'ylim');
                                xLimVals=get(gca,'xlim');
                                unitSpace=(yLimVals(2)-yLimVals(1))/20;
                                text('Position',[xLimVals(2)+(xLimVals(2)-xLimVals(1))/25 yLimVals(2)-unitSpace*cond*2],'FontSize',9,'String',[markerText,'  ',num2str(testContrast(closeCond(cond))),'%'],'Color',colmapText(closeCond(cond),:));
                            end
                        end
                    end
                    meanCPsessions=[];%mean CP across channels for each session and condition
                    sdCPsessions=[];%SD of CP across channels for each session and condition
                    meanCPmovingAve=[];%mean of mean CP across three sessions
                    for cond=1:length(closeCond)
                        for sessionInd=1:size(allChCP,1)
                            meanCPsessions(sessionInd,cond)=mean(allChCP{sessionInd,cond});
                            sdCPsessions(sessionInd,cond)=std(allChCP{sessionInd,cond});
                        end
                        figure(figAve)
                        subplot(length(animals),length(areas),animalInd+(areaInd-1)*2);
                        plot(1:size(meanCPsessions,1),meanCPsessions(:,cond),'Color',colmapText(closeCond(cond),:),'Marker','o','MarkerFaceColor',colmapText(closeCond(cond),:),'MarkerSize',3);hold on
%                         errorbar(1:size(meanCPsessions,1),meanCPsessions(:,cond),sdCPsessions(:,cond),'Color',colmapText(cond,:));hold on
                        for sessionInd=1:size(allChCP,1)-(numSmoothSessions-1)
                            meanCPmovingAve(sessionInd,cond)=mean(meanCPsessions(sessionInd:sessionInd+2,cond));
                        end                        
                        figure(figMovingAve)
                        subplot(length(animals),length(areas),animalInd+(areaInd-1)*2);
                        plot(1:size(meanCPmovingAve,1),meanCPmovingAve(:,cond),'Color',colmapText(closeCond(cond),:),'Marker','o','MarkerFaceColor',colmapText(closeCond(cond),:),'MarkerSize',3);hold on
                        set(gca,'box','off');
                        axis square
                    end
                    xlim([0 size(allChCP,1)-1]);
                    numSess=size(dataArray,1);
                    figHist=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.07, 0.4]); %
                    set(figHist, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                    for cond=1:length(closeCond)
                        earlyCPs=[];
                        lateCPs=[];
                        for numEarlyLate=1:5
                            earlyCPs=[earlyCPs allChCP{numEarlyLate,cond}];
                            lateCPs=[lateCPs allChCP{end-(numEarlyLate-1),cond}];
                        end
                        if cond<=length(closeCond)/2
                            [h p ci stats]=ttest2(earlyCPs,lateCPs,[],'right');%unpaired one-sided t-test
                        elseif cond>length(closeCond)/2
                            [h p ci stats]=ttest2(earlyCPs,lateCPs,[],'left');
                        end
                        if roving==1
                            ttestTable{animalInd,sampleContrastsInd}=[ttestTable{animalInd,areaInd};stats.df stats.tstat p];%t-test df, t value, and p-val
                        elseif roving==0
                            ttestTable{animalInd,areaInd}=[ttestTable{animalInd,areaInd};stats.df stats.tstat p];
                        end
                        figure(figHist)
                        subplot(length(closeCond),1,length(closeCond)-(cond-1));
                        [n,x]=hist(earlyCPs,0.2:0.05:0.8);
%                         stairs([x(1)-(x(2)-x(1))/2 x-(x(2)-x(1))/2 x(length(x))+(x(2)-x(1))/2],[0 n 0],'Color',colmapText(closeCond(cond),:))
                        stairs([x(1)-(x(2)-x(1))/2 x-(x(2)-x(1))/2 x(length(x))+(x(2)-x(1))/2],[0 n 0],'Color','k')
                        hold on
                        hist(lateCPs(:),0.2:0.05:0.8,'EdgeColor','none');
                        h = findobj(gca,'Type','patch');
%                         set(h,'FaceColor',colmapText(closeCond(cond),:),'facealpha',0.4,'EdgeColor','none');
                        set(h,'FaceColor',colmapText(closeCond(cond),:),'facealpha',0.8,'EdgeColor','none');
                        title(sprintf('%.5f',p));
%                         set(gca, 'box', 'on');
%                         if roving==0
%                             figure(figCPconds);
%                             subplot(length(closeCond),2*length(areas),1+((cond-1)*2*length(areas))+animalInd-1+2*(areaInd-1));
%                             %subplot(2*length(areas),length(closeCond),cond+length(closeCond)*(animalInd+2*(areaInd-1)-1));
%                         elseif roving==1
%                             figure(figCPconds(animalInd+(areaInd-1)*2));
%                             subplot(length(sampleContrasts),length(closeCond),cond+length(closeCond)*(sampleContrastsInd-1));
%                         end
%                         yFitted=0:numSess+1;
%                         if useRobustfit==1
%                             yFitted=b(1)+b(2)*yFitted;
%                         elseif useRobustfit==0
%                         end
%                         plot([0 numSess+1],[0.5 0.5],'k');
% %                         if stats.p(2)<.05%/length(closeCond)%Bonferroni correction for number of contrast conditions
%                             plot(0:numSess+1,yFitted,'r');
% %                         end
%                         goodChCParray=chCParray(~isnan(chCParray));
%                         if roving==0
%                             [pPointfive{animalInd,areaInd}(cond),h,stats]=signrank(goodChCParray,0.5);%test whether median differs significantly from 0.5
%                             statsPointfive{animalInd,areaInd}(cond)=stats;
%                             meanCP{animalInd,areaInd}(cond)=mean(goodChCParray);
%                         elseif roving==1
%                             [pPointfive{animalInd,sampleContrastsInd}(cond),h,stats]=signrank(goodChCParray,0.5);%test whether median differs significantly from 0.5
%                             statsPointfive{animalInd,sampleContrastsInd}(cond)=stats;
%                             meanCP{animalInd,sampleContrastsInd}(cond)=mean(goodChCParray);
%                         end
%                         for bootstrapInd=1:1000
%                             randSess=randperm(length(sessionNums));
%                             randSessArr=[];
%                             for sessionInd=1:size(chCParray,2)
%                                 randSessArr=[randSessArr randSess'];  
%                             end
%                             [bBoot(bootstrapInd,:),statsBoot{bootstrapInd}]=robustfit(randSessArr,chCParray);                            
%                         end
%                         meanBoot=mean(bBoot(:,2));
%                         varBoot=std(bBoot(:,2));
%                         zObserved=(bObserved(2)-meanBoot)/varBoot;%converted observed slope into z-score
                    end
                end
            end
        end
%         if roving==0
%             spearmanTableChs{areaInd,animalInd}=spearmanTableCh;
%             spearmanTableChRs{areaInd,animalInd}=spearmanTableChR;
%         elseif roving==1
%             spearmanTableChs{animalInd,sampleContrastsInd}=spearmanTableCh;
%             spearmanTableChRs{animalInd,sampleContrastsInd}=spearmanTableChR;            
%         end
%         for i=1:size(spearmanTableCh,1)
%             if sum(spearmanTableCh(i,:)<.05)==6||sum(spearmanTableCh(i,:)<.05/6*5)==5||sum(spearmanTableCh(i,:)<.05/6*4)==4||sum(spearmanTableCh(i,:)<.05/6*3)==3||sum(spearmanTableCh(i,:)<.05/6*2)==2||sum(spearmanTableCh(i,:)<.05/6)==1
%                 sigChsRs=[sigChsRs;{animal} {area} channels(i) spearmanTableCh(i,:) spearmanTableChR(i,:)];
%             end
%         end
    end
end
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        subplot(length(animals),length(areas),animalInd+(areaInd-1)*2);
        ylabel('choice probability');
        xlabel('sessions');
        title([animal,' ',area]);
    end
end
figure(figHist)
plot2svg
% if roving==0
%     bTableStatsV4=[];
%     bTableStatsV1=[];
%     for animalInd=1:2
%         bTableRow=[];
%         for cond=1:6
%             bTableRow=[bTableRow bObserved{animalInd,1}(cond) pSlope{animalInd,1}(cond)];
%         end
%         bTableStatsV4=[bTableStatsV4;bTableRow];
%     end
%     for animalInd=1:2
%         bTableRow=[];
%         for cond=1:4
%             bTableRow=[bTableRow bObserved{animalInd,2}(cond) pSlope{animalInd,2}(cond)];
%         end
%         bTableStatsV1=[bTableStatsV1;bTableRow];
%     end
%     sTableStatsV4=[];%signed rank test
%     sTableStatsV1=[];
%     for animalInd=1:2
%         sTableRow=[];
%         for cond=1:6
%             sTableRow=[sTableRow;meanCP{animalInd,1}(cond) statsPointfive{animalInd,1}(cond).zval pPointfive{animalInd,1}(cond)];
%         end
%         sTableStatsV4=[sTableStatsV4;sTableRow];
%     end
%     for animalInd=1:2
%         sTableRow=[];
%         for cond=1:4
%             sTableRow=[sTableRow;meanCP{animalInd,2}(cond) statsPointfive{animalInd,2}(cond).zval pPointfive{animalInd,2}(cond)];
%         end
%         sTableStatsV1=[sTableStatsV1;sTableRow];
%     end
% elseif roving==1    
%     bTableStatsV1=[];
%     for animalInd=1:2
%         bTableRow=[];
%         for cond=1:4
%             bTableRow=[bTableRow bObserved{animalInd,2}(cond) pSlope{animalInd,2}(cond)];
%         end
%         bTableStatsV1=[bTableStatsV1;bTableRow];
%     end
%     sTableStatsV1=[];%signed rank test
%     for animalInd=1:2
%         sTableRow=[];
%         for cond=1:4
%             sTableRow=[sTableRow;meanCP{animalInd,1}(cond) statsPointfive{animalInd,1}(cond).zval pPointfive{animalInd,1}(cond) meanCP{animalInd,2}(cond) statsPointfive{animalInd,2}(cond).zval pPointfive{animalInd,2}(cond) meanCP{animalInd,3}(cond) statsPointfive{animalInd,3}(cond).zval pPointfive{animalInd,3}(cond)];
%         end
%         sTableStatsV1=[sTableStatsV1;sTableRow];
%     end
% end
folderpathname=fullfile(rootFolder,'PL',analysisType);
if ~exist(folderpathname,'dir')
    mkdir(folderpathname);
end
if roving==0
    figure(figMovingAve);
    imagename='CPs_across_chs_nonroving_runningmean';
    pathname=fullfile(rootFolder,'PL',analysisType,imagename);
    printtext=sprintf('print -dpng %s.png',pathname);
    set(gcf,'PaperPositionMode','auto')
    eval(printtext);
    matname=['CPs_across_chs_nonroving','_',num2str(sampleContrast),'_runningmean'];
elseif roving==1
    for animalInd=1:length(animals)
        animal=animals{animalInd};
        for areaInd=1:length(areas)
            area=areas{areaInd};
            figure(animalInd+(areaInd-1)*2);
            imagename=['CPs_across_chs_',animal,'_',area,'_roving_runningmean'];
            pathname=fullfile(rootFolder,'PL',analysisType,imagename);
            printtext=sprintf('print -dpng %s.png',pathname);
            set(gcf,'PaperPositionMode','auto')
            eval(printtext);
        end
    end
    matname=['CPs_across_chs_roving_runningmean'];
end
pathname=fullfile(rootFolder,'PL',analysisType,matname);
saveText=['save ',pathname,'.mat statsObserved'];
eval(saveText);