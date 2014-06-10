function read_bj_cp_batch(animal,area,roving)
%Written by Xing 27/06/13
%to calculate and write choice probability correlation coefficients
useRobustfit=1;
onExternalHD=0;
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
    spearmanTable=cell(length(animals),length(areas));
    spearmanTableChs=cell(length(animals),length(areas));
    spearmanTableChRs=cell(length(animals),length(areas));
elseif roving==1
    for animalInd=1:2
        for areaInd=1:2
            figCPconds(animalInd+(areaInd-1)*2)=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
            set(figCPconds(animalInd+(areaInd-1)*2), 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
        end
    end
    spearmanTable=cell(length(animals),3);
    spearmanTableChs=cell(length(animals),3);
    spearmanTableChRs=cell(length(animals),3);
end
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        sessionNums = main_raw_sessions_final(animal,area,[],0);
        channels=main_channels(animal,area);
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
                        %figure with data for each channel in separate plots:
                        figure(figChCPconds);
                        subplot(ceil(length(channels)/5),5,i);
                        for cond=1:length(closeCond)
                            for sessionInd=1:size(dataArray,1)
                                if length(dataArray{sessionInd,3})==length(closeCond)
                                    plot(sessionInd,dataArray{sessionInd,3}(cond),'LineStyle','none','Marker',MarkerType,'Color',colmapText(closeCond(cond),:));hold on
                                end
                            end
                        end
                        xlim([0 length(sessionNums)+1]);
                        title(num2str(channels(i)));
                    end
                    if roving==0
                        imagename=['CPs_chs_',animal,'_',area,'_',num2str(sampleContrast)];
                    elseif roving==1
                        imagename=['CPs_chs_',animal,'_',area,'_',num2str(sampleContrast)];
                    end
                    folderpathname=fullfile(rootFolder,'PL',analysisType);
                    if ~exist(folderpathname,'dir')
                        mkdir(folderpathname);
                    end
                    pathname=fullfile(rootFolder,'PL',analysisType,imagename);
                    printtext=sprintf('print -dpng %s.png',pathname);
                    set(gcf,'PaperPositionMode','auto')
                    eval(printtext);
                    numSess=size(dataArray,1);
                    for cond=1:length(closeCond)
                        chCParray=[];
                        for sessionInd=1:size(dataArray,1)
                            if length(allChCP{sessionInd,cond})==length(channels)
                                chCParray=[chCParray;allChCP{sessionInd,cond}];
                            end
                        end
                        allChCParray{cond}=chCParray;
                        sessionArr=[];
                        for sessionInd=1:size(chCParray,2)
                            sessionArr=[sessionArr [1:size(chCParray,1)]'];
                        end
                        chCParray=reshape(chCParray,1,size(chCParray,1)*size(chCParray,2));
                        sessionArr=reshape(sessionArr,1,size(sessionArr,1)*size(sessionArr,2));
                        [b,stats]=robustfit(sessionArr,chCParray);
                        bObserved{animalInd,areaInd}(cond)=b(2);
                        statsObserved{animalInd,areaInd}=[statsObserved{animalInd,areaInd};b(2) stats.dfe stats.t(2) stats.p(2)];
                        pSlope{animalInd,areaInd}(cond)=stats.p(2);
                        [rho p]=corr(chCParray',sessionArr','type','Spearman');
                        if roving==1
                            spearmanTable{animalInd,sampleContrastsInd}=[spearmanTable{animalInd,areaInd};length(chCParray)-2 rho p];
                        elseif roving==0
                            spearmanTable{animalInd,areaInd}=[spearmanTable{animalInd,areaInd};length(chCParray)-2 rho p];
                        end
                        if roving==0
                            figure(figCPconds);
                            subplot(length(closeCond),2*length(areas),1+((cond-1)*2*length(areas))+animalInd-1+2*(areaInd-1));
                            %subplot(2*length(areas),length(closeCond),cond+length(closeCond)*(animalInd+2*(areaInd-1)-1));
                        elseif roving==1
                            figure(figCPconds(animalInd+(areaInd-1)*2));
                            subplot(length(sampleContrasts),length(closeCond),cond+length(closeCond)*(sampleContrastsInd-1));
                        end
                        yFitted=0:numSess+1;
                        if useRobustfit==1
                            yFitted=b(1)+b(2)*yFitted;
                        elseif useRobustfit==0
                        end
                        plot([0 numSess+1],[0.5 0.5],'k');
%                         if stats.p(2)<.05%/length(closeCond)%Bonferroni correction for number of contrast conditions
                            plot(0:numSess+1,yFitted,'r');
%                         end
                        goodChCParray=chCParray(~isnan(chCParray));
                        if roving==0
                            [pPointfive{animalInd,areaInd}(cond),h,stats]=signrank(goodChCParray,0.5);%test whether median differs significantly from 0.5
                            statsPointfive{animalInd,areaInd}(cond)=stats;
                            meanCP{animalInd,areaInd}(cond)=mean(goodChCParray);
                        elseif roving==1
                            [pPointfive{animalInd,sampleContrastsInd}(cond),h,stats]=signrank(goodChCParray,0.5);%test whether median differs significantly from 0.5
                            statsPointfive{animalInd,sampleContrastsInd}(cond)=stats;
                            meanCP{animalInd,sampleContrastsInd}(cond)=mean(goodChCParray);
                        end
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
        if roving==0
            spearmanTableChs{areaInd,animalInd}=spearmanTableCh;
            spearmanTableChRs{areaInd,animalInd}=spearmanTableChR;
        elseif roving==1
            spearmanTableChs{animalInd,sampleContrastsInd}=spearmanTableCh;
            spearmanTableChRs{animalInd,sampleContrastsInd}=spearmanTableChR;            
        end
        for i=1:size(spearmanTableCh,1)
            if sum(spearmanTableCh(i,:)<.05)==6||sum(spearmanTableCh(i,:)<.05/6*5)==5||sum(spearmanTableCh(i,:)<.05/6*4)==4||sum(spearmanTableCh(i,:)<.05/6*3)==3||sum(spearmanTableCh(i,:)<.05/6*2)==2||sum(spearmanTableCh(i,:)<.05/6)==1
                sigChsRs=[sigChsRs;{animal} {area} channels(i) spearmanTableCh(i,:) spearmanTableChR(i,:)];
            end
        end
    end
end
if roving==0
    bTableStatsV4=[];
    bTableStatsV1=[];
    for animalInd=1:2
        bTableRow=[];
        for cond=1:6
            bTableRow=[bTableRow bObserved{animalInd,1}(cond) pSlope{animalInd,1}(cond)];
        end
        bTableStatsV4=[bTableStatsV4;bTableRow];
    end
    for animalInd=1:2
        bTableRow=[];
        for cond=1:4
            bTableRow=[bTableRow bObserved{animalInd,2}(cond) pSlope{animalInd,2}(cond)];
        end
        bTableStatsV1=[bTableStatsV1;bTableRow];
    end
    sTableStatsV4=[];%signed rank test
    sTableStatsV1=[];
    for animalInd=1:2
        sTableRow=[];
        for cond=1:6
            sTableRow=[sTableRow;meanCP{animalInd,1}(cond) statsPointfive{animalInd,1}(cond).zval pPointfive{animalInd,1}(cond)];
        end
        sTableStatsV4=[sTableStatsV4;sTableRow];
    end
    for animalInd=1:2
        sTableRow=[];
        for cond=1:4
            sTableRow=[sTableRow;meanCP{animalInd,2}(cond) statsPointfive{animalInd,2}(cond).zval pPointfive{animalInd,2}(cond)];
        end
        sTableStatsV1=[sTableStatsV1;sTableRow];
    end
elseif roving==1    
    bTableStatsV1=[];
    for animalInd=1:2
        bTableRow=[];
        for cond=1:4
            bTableRow=[bTableRow bObserved{animalInd,2}(cond) pSlope{animalInd,2}(cond)];
        end
        bTableStatsV1=[bTableStatsV1;bTableRow];
    end
    sTableStatsV1=[];%signed rank test
    for animalInd=1:2
        sTableRow=[];
        for cond=1:4
            sTableRow=[sTableRow;meanCP{animalInd,1}(cond) statsPointfive{animalInd,1}(cond).zval pPointfive{animalInd,1}(cond) meanCP{animalInd,2}(cond) statsPointfive{animalInd,2}(cond).zval pPointfive{animalInd,2}(cond) meanCP{animalInd,3}(cond) statsPointfive{animalInd,3}(cond).zval pPointfive{animalInd,3}(cond)];
        end
        sTableStatsV1=[sTableStatsV1;sTableRow];
    end
end
folderpathname=fullfile(rootFolder,'PL',analysisType);
if ~exist(folderpathname,'dir')
    mkdir(folderpathname);
end
if roving==0
    figure(figCPconds);
    imagename='CPs_across_chs_nonroving';
    pathname=fullfile(rootFolder,'PL',analysisType,imagename);
    printtext=sprintf('print -dpng %s.png',pathname);
    set(gcf,'PaperPositionMode','auto')
    eval(printtext);
    matname=['CPs_across_chs_nonroving','_',num2str(sampleContrast)];
elseif roving==1
    for animalInd=1:length(animals)
        animal=animals{animalInd};
        for areaInd=1:length(areas)
            area=areas{areaInd};
            figure(animalInd+(areaInd-1)*2);
            imagename=['CPs_across_chs_',animal,'_',area,'_roving'];
            pathname=fullfile(rootFolder,'PL',analysisType,imagename);
            printtext=sprintf('print -dpng %s.png',pathname);
            set(gcf,'PaperPositionMode','auto')
            eval(printtext);
        end
    end
    matname=['CPs_across_chs_roving'];
end
pathname=fullfile(rootFolder,'PL',analysisType,matname);
saveText=['save ',pathname,'.mat statsObserved'];
eval(saveText);