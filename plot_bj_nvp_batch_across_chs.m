function plot_bj_nvp_batch_across_chs(roving,areas,psychoType)

%written 29/03/13, modified from plot_bj_crf_or_roc_batch_across_chs.
%Reads neuro- and psychometric threshold values for all conditions, to calculate averages and
%examine correlation between threshold and session.
%input arg for roving data: areas={'v1_2_1'} or {'v1_2_2'} or {'v1_2_3'}   
%Input arg 'psychoType' is either 'psycho' (Weibull function curve fitting
%performed without 'lambda' associational learning/ attention modulation
%parameter, or 'psycho_param' (Weibull curve fitting with lambda
%parameter).
onExternalHD=0;
if onExternalHD==1
    rootFolder='G:\PL_backup_060413';
else
    rootFolder='F:';
end
plotPsychoOnly=1;
analysisType='NVP';
animalTexts=[{'subject B'} {'subject J'}];
animals=[{'blanco'} {'jack'}];
if roving==0
    areaTexts=[{'V4'} {'V1'}];
    areas=[{'v4_1'} {'v1_1'}];
elseif roving==1
    areaTexts={'20' '30' '40'};
end
calculateStats=1;
fig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.55, 0.8]);
set(fig,'PaperUnits','centimeters','PaperType','A4','PaperOrientation', 'portrait', 'PaperPosition', [0.63452 0.63452 21/8*5.5 28.41]);
allcReshape=[];
allcondReshape=[];
allsessHalfReshape=[];
allsubjectReshape=[];
allareaReshape=[];
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        channels = main_channels(animal,area);
        sessions = main_raw_sessions_final(animal,area,[],0);
        if strcmp(psychoType,'psycho_param')
            sessions = main_raw_sessions_final_psycho(animal,area,[],0);
        end
        if ~strcmp(psychoType,'psycho_param')
            includeMatch=[];
            excludeSessions=[26 50 306 312 316 322:328 342];
            for includeInd=1:length(sessions)
                if sum(excludeSessions==sessions(includeInd))==0
                    includeMatch=[includeMatch includeInd];
                end
            end
            sessions=sessions(includeMatch);
        end
        [sampleContrasts testContrasts]=area_metadata(area);
        if strcmp(area,'v4_1')||strcmp(area,'v1_1')
            test_epochs={0 529 529*2 529*3};durSpon=150;
        elseif strncmp(area,'v1_2',4)
            test_epochs={0 512 512*2 512*3};durSpon=150;
        end
        
        psychoname=['psycho_constants_',area];
        psychoPathname=fullfile(rootFolder,'PL','psycho_data',animal,psychoname);
        subplotTitleText={'neurometric threshold lower contrast' 'neurometric threshold higher contrast' 'psychometric threshold lower contrast' 'psychometric threshold higher contrast'};
        
        %read in psycho data
        %
        % if strcmp(folder,'F:\blanco\v4_1_roc_analysis')
        %     markerCol='r';
        %     divisor1=17;
        % elseif strcmp(folder,'F:\blanco\v1_roc_analysis\multiple_time_periods')
        %     markerCol='b';
        %     divisor1=17;
        % else
        %     markerCol='b';
        %     divisor1=15;
        % end
        setFill='none';
        markerSh='s';
        currcol = colormap(winter);
        margin=7;
        markerCols='rbrb';
        figure(fig);
        cellCategory=zeros(length(channels),2)-2;
        cellCategory(:,1)=channels;
        for sampCond=1:length(sampleContrasts)
            sampleContrast=sampleContrasts(sampCond);
            testContrast=testContrasts(sampCond,:);
            overThreshold=[sampleContrast 100-sampleContrast sampleContrast 100-sampleContrast];
            appendText=['_',num2str(sampleContrast)];
            for epoch=1:size(test_epochs,2)
                if strcmp(analysisType,'NVP')&&epoch==4
                    if epoch==1
                        periods=[-durSpon 0];
                    else
                        periods=[test_epochs{epoch-1} test_epochs{epoch}(1)];
                    end
                    for subPeriod=1:length(periods)-1
                        startEndTime=['_',num2str(periods(subPeriod)),'_to_',num2str(periods(subPeriod+1))];
                        psychoThresholdMatName=[area,appendText,'wholetrial_psyThreshold'];
                        psychoThresholdMatFolder=fullfile(rootFolder,'PL',psychoType,animal,'psyThreshold_mat');
                        psychoThresholdMatPathname=fullfile(psychoThresholdMatFolder,psychoThresholdMatName);
                        loadText=['load ',psychoThresholdMatPathname,' sessionSorted2 threshold82lower threshold82higher'];
                        eval(loadText)
                        thresholdsNL=zeros(1,length(sessionSorted2));
                        thresholdsNH=zeros(1,length(sessionSorted2));
                        thresholdsPL=zeros(1,length(sessionSorted2));
                        thresholdsPH=zeros(1,length(sessionSorted2));
                        sessionTally=zeros(1,length(sessionSorted2));
                        allThresholdsNL=cell(1,length(sessionSorted2));
                        allThresholdsNH=cell(1,length(sessionSorted2));
                        allThresholdsPL=cell(1,length(sessionSorted2));
                        allThresholdsPH=cell(1,length(sessionSorted2));
                        Ncells=0;
                        thresholdsPL=threshold82lower;
                        thresholdsPH=threshold82higher;
                        for i=1:length(channels)
                            chNum=channels(i);
                            thMatname=[num2str(chNum),appendText,startEndTime,'_nvpThreshold_',area];
                            thMatFolder=fullfile(rootFolder,'PL',analysisType,animal,'nvpThreshold_mat');
                            thMatPathname=fullfile(thMatFolder,thMatname);
                            loadText=['load ',thMatPathname,' sessionSorted1 threshold82lower threshold82higher'];
                            eval(loadText)
                            if ~strcmp(psychoType,'psycho_param')
                                Ncells=Ncells+1;
                                includeMatch=[];
                                for includeInd=1:length(sessionSorted1)
                                    if sum(excludeSessions==sessionSorted1(includeInd))==0
                                        includeMatch=[includeMatch includeInd];
                                    end
                                end
                                sessionSorted1=sessionSorted1(:,includeMatch);
                            end
                            for sessionInd=1:length(sessionSorted2)%for each good psycho session, see if corresponding neuro session exists
                                colInd=find(sessionSorted2(sessionInd)==sessionSorted1);
                                if ~isempty(colInd)
                                    sessionTally(1,sessionInd)=sessionTally(1,sessionInd)+1;
                                    thresholdsNL(1,sessionInd)=thresholdsNL(1,sessionInd)+threshold82lower(colInd);%neurometric
                                    thresholdsNH(1,sessionInd)=thresholdsNH(1,sessionInd)+threshold82higher(colInd);
                                    
                                    allThresholdsNL{1,sessionInd}=[allThresholdsNL{1,sessionInd} threshold82lower(colInd)];
                                    allThresholdsNH{1,sessionInd}=[allThresholdsNH{1,sessionInd} threshold82higher(colInd)];
                                    allThresholdsPL{1,sessionInd}=thresholdsPL(sessionInd);
                                    allThresholdsPH{1,sessionInd}=thresholdsPH(sessionInd);
                                end
                            end
                        end
                        thresholdsNL=thresholdsNL./sessionTally;
                        thresholdsNH=thresholdsNH./sessionTally;
                        allData=[];
                        allSessionData=[];
                        allXData=[];
                        allData=[allData;{thresholdsNL} {thresholdsNH} {thresholdsPL} {thresholdsPH}];
                        allSessionData=[allSessionData;{sessions} {sessions} {sessions} {sessions}];
                        allXData=[allXData;{1:length(sessions)} {1:length(sessions)} {1:length(sessions)} {1:length(sessions)}];
                        startComparison=1;
                        if strcmp(psychoType,'psycho_param')
                            startComparison=3;
                        end
                        for comparison=startComparison:4
                            temp=isnan(allData{1,comparison});
                            allData{1,comparison}=allData{1,comparison}(~temp);
                            allSessionData{1,comparison}=allSessionData{1,comparison}(~temp);
                            allXData{1,comparison}=allXData{1,comparison}(~temp);
                            c=[(allSessionData{1,comparison})' (allData{1,comparison})'];%session in column 1, data in column 2
                            [a b CI,stats]=corrcoef(c);
                            coefficientsCat(animalInd+2*(areaInd-1),comparison)=a(2);%1st row: B V4, 2nd: J V4, 3rd: B V1, 4th: J V1; 1st column: NL, 2nd: NH, 3rd: PL, 4th: PH
                            pCat(animalInd+2*(areaInd-1),comparison)=b(2);
                            statsCat{animalInd+2*(areaInd-1),comparison}=stats;
                            [a b]=corr(c,'type','Spearman');
                            rho(animalInd+2*(areaInd-1),comparison)=a(2);
                            pSpearman(animalInd+2*(areaInd-1),comparison)=b(2);
                        end
                        c=[(allData{1,3})' (allData{1,4})'];%session in column 1, data in column 2
                        condList=[zeros(size((allData{1,comparison})',1),size((allData{1,comparison})',2))+1 zeros(size((allData{1,comparison})',1),size((allData{1,comparison})',2))+2];
                        sessHalfList=[zeros(ceil(size(c,1)/2),size(c,2))+1;zeros(floor(size(c,1)/2),size(c,2))+2];%divide sessions into first half and second half
                        cReshape=reshape(c,1,size(c,1)*size(c,2));
                        condReshape=reshape(condList,1,size(condList,1)*size(condList,2));%1:lower; 2: higher
                        sessHalfReshape=reshape(sessHalfList,1,size(sessHalfList,1)*size(sessHalfList,2));
                        [p,table,stats]=anovan(cReshape,{condReshape,sessHalfReshape},'model','interaction')%check for differences in threshold based on condition type (higher or lower test than sample contrast) and on whether session belongs to first or second half of training
                        pCondHalfANOVA{areaInd,animalInd}=p;
                        statsCondHalf{areaInd,animalInd}=stats;
                        if roving==1
                            for tableInd=1:4
                                grandCorrTable((sampCond-1)*2+animalInd,3+(tableInd-1)*4)=coefficientsCat(animalInd,tableInd);
                                grandCorrTable((sampCond-1)*2+animalInd,1+(tableInd-1)*4)=pCat(animalInd,tableInd);
                            end
                        end
                        subjectList=zeros(1,length(cReshape))+animalInd;
                        areaList=zeros(1,length(cReshape))+areaInd;
                        allcReshape=[allcReshape cReshape];
                        allcondReshape=[allcondReshape condReshape];
                        allsessHalfReshape=[allsessHalfReshape sessHalfReshape];
                        allsubjectReshape=[allsubjectReshape subjectList];
                        allareaReshape=[allareaReshape areaList];
                        if animalInd==length(animals)&&areaInd==length(areas)
                            [allp,alltable,allstats]=anovan(allcReshape,{allcondReshape,allsessHalfReshape,allsubjectReshape,allareaReshape},'model','full')%check for differences in threshold based on condition type (higher or lower test than sample contrast) and on whether session belongs to first or second half of training
                        end
                        subplotRemap=[1 2 3 4];
                        figure(fig);
                        if roving==0
                            subplot(length(areas),2,animalInd+2*(areaInd-1));
                        elseif roving==1
                            subplot(length(sampleContrasts),2,animalInd+2*(sampCond-1));
                        end
                        for subplotInd=1:length(subplotRemap)
                            %                     subplot(2,2,subplotRemap(subplotInd));
                            if ~isempty(allData{1,subplotInd})
                                if subplotInd<=2
                                    if plotPsychoOnly==0
                                        plot(allXData{1,subplotInd},allData{1,subplotInd},'Marker','o','Color',markerCols(subplotInd),'LineStyle','none','MarkerFaceColor','none');hold on
                                    end
                                else
                                    for sessionNum=1:length(allData{1,subplotInd})
                                        if allData{1,subplotInd}(sessionNum)==overThreshold(subplotInd)                                            
                                            plot(allXData{1,subplotInd}(sessionNum),allData{1,subplotInd}(sessionNum),'Marker','o','Color',markerCols(subplotInd),'LineStyle','none','MarkerFaceColor','none');hold on
                                        else
                                            plot(allXData{1,subplotInd}(sessionNum),allData{1,subplotInd}(sessionNum),'Marker','o','Color',markerCols(subplotInd),'LineStyle','none','MarkerFaceColor',markerCols(subplotInd));hold on
                                        end
                                    end
                                end
                            end
                            xlim([0 length(allData{1,subplotInd})+1]);
                        end
                        if strcmp(animal,'blanco')&&strcmp(area,'v4_1')
                            subplot(2,2,1);
                            plot(3,9.96,'Marker','o','Color','b','LineStyle','none','MarkerFaceColor','b')
                            ylim([0 10]);
                            set(gca,'YTick',[0:10],'YTickLabel',[{'0'} {'1'} {'2'} {'3'} {'4'} {'5'} {'6'} {'7'} {'8'} {''} {'23'}]);
                            text(0,9,'//','FontSize',15)
                        end
                        if strcmp(animal,'jack')&&strcmp(area,'v1_1')
                            subplot(2,2,4);
                            plot(1,30.0,'Marker','o','Color','b','LineStyle','none','MarkerFaceColor','none')
                            ylim([0 30]);
                            set(gca,'YTick',[0 5 10 15 20 25 30],'YTickLabel',[{'0'} {'5'} {'10'} {'15'} {'20'} {''} {'70'}]);
                            text(0,25,'//','FontSize',15)
                        end
                        for subplotInd=1:length(subplotRemap)
                            %                     subplot(2,2,subplotRemap(subplotInd));
                            %                     if subplotRemap(subplotInd)==1
                            %                         ptext=sprintf('N = %s channels',num2str(Ncells));
                            %                         yLimVals=get(gca,'YLim');
                            %                         text('Position',[-2 yLimVals(2)+0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
                            %                     end
                            %                     if pCat(1,subplotInd)<0.05
                            %                         colText='g';
                            %                     else
                            %                         colText=markerCols(1);
                            %                     end
                            %                     ptext=sprintf('r= %f  p= %f',coefficientsCat(1,subplotInd),pCat(1,subplotInd));
                            %                     yLimVals=get(gca,'YLim');
                            %                     xLimVals=get(gca,'XLim');
                            %                     remapText=[3 2 1];
                            %                     text('Position',[xLimVals(1)+(xLimVals(2)-xLimVals(1))/2 yLimVals(1)+0.05*remapText(1)*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext,'Color',colText);
                            %                     title(subplotTitleText{subplotInd});
                            axis square
                            ylimits=get(gca,'YLim');
                            ylim([0 ylimits(2)]);
                            if animalInd==1
                                if roving==0
                                    ptext=areaTexts{areaInd};
                                elseif roving==1
                                    ptext=areaTexts{sampCond};
                                end
                                xlimits=get(gca,'XLim');
                                text('Position',[xlimits(1)-(xlimits(2)-xlimits(1))/4 (ylimits(2)-ylimits(1))/2],'FontSize',9,'String',ptext);
                                ylabel('threshold');
                            end
                            if areaInd==1&&sampCond==1
                                title(animalTexts{animalInd});
                                if length(areas)==1
                                    xlabel('session');
                                end
                            else
                                xlabel('session');
                            end
                        end
                    end
                end
            end
            if calculateStats==1
                if plotPsychoOnly==0
                    anovaVals=[thresholdsNL thresholdsNH thresholdsPL thresholdsPH];
                    nvpFactor=[zeros(1,length(thresholdsNL))+1 zeros(1,length(thresholdsNH))+1 zeros(1,length(thresholdsPL))+2 zeros(1,length(thresholdsPH))+2];
                    highlowFactor=[zeros(1,length(thresholdsNL))+1 zeros(1,length(thresholdsNH))+2 zeros(1,length(thresholdsPL))+1 zeros(1,length(thresholdsPH))+2];
                    %         sessionFactor=[zeros(1,(floor(length(thresholdsNL)/2)))+1 zeros(1,(ceil(length(thresholdsNL)/2)))+2 zeros(1,(floor(length(thresholdsNH)/2)))+1 zeros(1,(ceil(length(thresholdsNH)/2)))+2 zeros(1,(floor(length(thresholdsPL)/2)))+1 zeros(1,(ceil(length(thresholdsPL)/2)))+2 zeros(1,(floor(length(thresholdsPH)/2)))+1 zeros(1,(ceil(length(thresholdsPH)/2)))+2];
                    %         sessionFactor=[zeros(1,(floor(length(thresholdsNL)/3)))+1 zeros(1,(floor(length(thresholdsNL)/3)))+2 zeros(1,length(thresholdsNL)-(2*floor(length(thresholdsNL)/3)))+3 zeros(1,(floor(length(thresholdsNH)/3)))+1 zeros(1,(floor(length(thresholdsNH)/3)))+2 zeros(1,length(thresholdsNH)-(2*floor(length(thresholdsNH)/3)))+3 zeros(1,(floor(length(thresholdsPL)/3)))+1 zeros(1,(floor(length(thresholdsPL)/3)))+2 zeros(1,length(thresholdsPL)-(2*floor(length(thresholdsPL)/3)))+3 zeros(1,(floor(length(thresholdsPH)/3)))+1 zeros(1,(floor(length(thresholdsPH)/3)))+2 zeros(1,length(thresholdsPH)-(2*floor(length(thresholdsPH)/3)))+3];
                    %         [p,table,stats]=anovan(anovaVals,{nvpFactor highlowFactor sessionFactor},'model','full');anovap{animalInd+2*(areaInd-1)}=p;anovatable{animalInd+2*(areaInd-1)}=table;anovastats{animalInd+2*(areaInd-1)}=stats;
                    [p,table,stats]=anovan(anovaVals,{nvpFactor highlowFactor},'model','full');anovap{sampCond,animalInd+2*(areaInd-1)}=p;anovatable{sampCond,animalInd+2*(areaInd-1)}=table;anovastats{sampCond,animalInd+2*(areaInd-1)}=stats;
                    figure;
                    [c,m,h,nms]=multcompare(stats)
                    %%1st row: B V4, 2nd: B V1, 3rd: J V4, 4th: J V1; 1st column: NL vs NH, 2nd: PL vs PH, 3rd: NL vs PL, 4th: NH vs PH
                    [h,p,ci,stats]=ttest(thresholdsNL,thresholdsNH);allh(1,animalInd+2*(areaInd-1))=h;allp(1,animalInd+2*(areaInd-1))=p;allci{1,animalInd+2*(areaInd-1)}=ci;allstats{1,animalInd+2*(areaInd-1)}=stats;
                    [h,p,ci,stats]=ttest(thresholdsPL,thresholdsPH);allh(2,animalInd+2*(areaInd-1))=h;allp(2,animalInd+2*(areaInd-1))=p;allci{2,animalInd+2*(areaInd-1)}=ci;allstats{2,animalInd+2*(areaInd-1)}=stats;
                    [h,p,ci,stats]=ttest(thresholdsNL,thresholdsPL);allh(3,animalInd+2*(areaInd-1))=h;allp(3,animalInd+2*(areaInd-1))=p;allci{3,animalInd+2*(areaInd-1)}=ci;allstats{3,animalInd+2*(areaInd-1)}=stats;
                    [h,p,ci,stats]=ttest(thresholdsNH,thresholdsPH);allh(4,animalInd+2*(areaInd-1))=h;allp(4,animalInd+2*(areaInd-1))=p;allci{4,animalInd+2*(areaInd-1)}=ci;allstats{4,animalInd+2*(areaInd-1)}=stats;
                else
                    [h,p,ci,stats]=ttest(thresholdsPL,thresholdsPH);allhttest(2,animalInd+2*(areaInd-1))=h;allpttest(2,animalInd+2*(areaInd-1))=p;allcittest{2,animalInd+2*(areaInd-1)}=ci;allstatsttest{2,animalInd+2*(areaInd-1)}=stats;
                    thresholdsPH1=thresholdsPH(1:ceil(length(thresholdsPH)/2));
                    thresholdsPH2=thresholdsPH(end-floor(length(thresholdsPH)/2)+1:end);
                    thresholdsPL1=thresholdsPL(1:ceil(length(thresholdsPL)/2));
                    thresholdsPL2=thresholdsPL(end-floor(length(thresholdsPL)/2)+1:end);
                    diffThresholds=thresholdsPH1-thresholdsPL1;
                    orderedDiff=sort(diffThresholds);
                    medianDiff1(animalInd+2*(areaInd-1))=orderedDiff(ceil(length(diffThresholds)/2));
                    meanDiff1(animalInd+2*(areaInd-1))=mean(diffThresholds);
                    diffThresholds=thresholdsPH2-thresholdsPL2;
                    orderedDiff=sort(diffThresholds);
                    medianDiff2(animalInd+2*(areaInd-1))=orderedDiff(ceil(length(diffThresholds)/2));
                    meanDiff2(animalInd+2*(areaInd-1))=mean(diffThresholds);
                    [p,h,stats]=signrank(thresholdsPL1,thresholdsPH1);allh(1,animalInd+2*(areaInd-1))=h;allp(1,animalInd+2*(areaInd-1))=p;
                    if isfield(stats,'zval')
                        allstats(1,animalInd+2*(areaInd-1))=stats.zval;
                    end
                    [p,h,stats]=signrank(thresholdsPL2,thresholdsPH2);allh(2,animalInd+2*(areaInd-1))=h;allp(2,animalInd+2*(areaInd-1))=p;
                    if isfield(stats,'zval')
                        allstats(2,animalInd+2*(areaInd-1))=stats.zval;
                    end
                    [p,h,stats]=ranksum(thresholdsPL1,thresholdsPL2);allhranksum(1,animalInd+2*(areaInd-1))=h;allpranksum(1,animalInd+2*(areaInd-1))=p;
                    if isfield(stats,'zval')
                        allstatsranksum(1,animalInd+2*(areaInd-1))=stats.zval;
                    end
                    [p,h,stats]=ranksum(thresholdsPH1,thresholdsPH2);allhranksum(2,animalInd+2*(areaInd-1))=h;allpranksum(2,animalInd+2*(areaInd-1))=p;
                    if isfield(stats,'zval')
                        allstatsranksum(2,animalInd+2*(areaInd-1))=stats.zval;
                    end
                    [p,h,stats]=ranksum(thresholdsPL1,thresholdsPH1);allhranksum(3,animalInd+2*(areaInd-1))=h;allpranksum(3,animalInd+2*(areaInd-1))=p;
                    if isfield(stats,'zval')
                        allstatsranksum(3,animalInd+2*(areaInd-1))=stats.zval;
                    end
                    [p,h,stats]=ranksum(thresholdsPL2,thresholdsPH2);allhranksum(4,animalInd+2*(areaInd-1))=h;allpranksum(4,animalInd+2*(areaInd-1))=p;
                    if isfield(stats,'zval')
                        allstatsranksum(4,animalInd+2*(areaInd-1))=stats.zval;
                    end
                    [p,h,stats]=ranksum(thresholdsPL1,thresholdsPH2);allhranksum(5,animalInd+2*(areaInd-1))=h;allpranksum(5,animalInd+2*(areaInd-1))=p;
                    if isfield(stats,'zval')
                        allstatsranksum(5,animalInd+2*(areaInd-1))=stats.zval;
                    end
                    [p,h,stats]=ranksum(thresholdsPL2,thresholdsPH1);allhranksum(6,animalInd+2*(areaInd-1))=h;allpranksum(6,animalInd+2*(areaInd-1))=p;
                    if isfield(stats,'zval')
                        allstatsranksum(6,animalInd+2*(areaInd-1))=stats.zval;
                    end
                end
            end
            if roving==1
                grandANOVATable(1:3,(sampCond*2)-1+6*(animalInd-1))=anovap{sampCond,animalInd};
            end
        end
    end
end
if plotPsychoOnly==1
    tableRanksum=[allstatsranksum(:,1)';allpranksum(:,1)';allstatsranksum(:,3)';allpranksum(:,3)';allstatsranksum(:,2)';allpranksum(:,2)';allstatsranksum(:,4)';allpranksum(:,4)'];
end

if roving==0
    saveCoefImageName=[areaTexts{1},'_',areaTexts{2},'_',analysisType,'_coefs',startEndTime];
elseif roving==1
    saveCoefImageName=[areaTexts{1},'_',analysisType,'_coefs',startEndTime];
end
saveCoefImageFolder=fullfile(rootFolder,'PL',analysisType);
saveCoefImagePath=fullfile(saveCoefImageFolder,saveCoefImageName);
printtext=sprintf('print -dpng %s.png',saveCoefImagePath);
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperPositionMode', 'auto');
eval(printtext);
