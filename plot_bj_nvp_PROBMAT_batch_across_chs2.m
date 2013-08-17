function plot_bj_nvp_PROBMAT_batch_across_chs2(roving,areas,psychoType)

%Written 15/08/13, modified from plot_bj_nvp_PROBMAT_batch_across_chs (function that plotted both neuro- and psychometric thresholds on same figure).
%Reads and plots either just the neurometric or just the psychometric threshold values for all conditions, based
%on PROBMAT values for population-wide data. Converts threshold values according to sample contrast.
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
plotSameGraph=1;%plot flankerless and flanker data on same subplot
flankerlessOnly=0;
sglroc3IndividualChs=1;
excludeSuppressed=0;
normalize=0;
useISI=0;
cutoff=1;
plotPsychoOnly=0;%set to 0 to plot neuro, set to 1 to plot psycho
analysisType='NVP';
animalTexts=[{'Monkey 1'} {'Monkey 2'}];
animals=[{'blanco'} {'jack'}];
if roving==0
    areaTexts=[{'V4'} {'V1'}];
    areas=[{'v4_1'} {'v1_1'}];
elseif roving==1
    areaTexts={'20' '30' '40'};
    areas=[{'v1_2_1'} {'v1_2_2'}];
    areas=[{'v1_2_2'}];
end
calculateStats=0;
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
        if roving==0
            if animalInd==1&&areaInd==1
                fig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.55, 0.8]);
                set(fig,'PaperUnits','centimeters','PaperType','A4','PaperOrientation', 'portrait', 'PaperPosition', [0.63452 0.63452 21/8*5.5 28.41]);
            end
        elseif roving==1
            if animalInd==1
                if plotSameGraph==0
                    fig(areaInd)=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.4, 0.8]);
                    set(fig(areaInd),'PaperUnits','centimeters','PaperType','A4','PaperOrientation', 'portrait', 'PaperPosition', [0.63452 0.63452 21/8*5.5 28.41]);
                elseif plotSameGraph==1%plot different stages of roving sessions adjacent to each other
                    if areaInd==1
                        fig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.55, 0.8]);
                        set(fig,'PaperUnits','centimeters','PaperType','A4','PaperOrientation', 'portrait', 'PaperPosition', [0.63452 0.63452 21/8*5.5 28.41]);
                    end
                end
            else
                if plotSameGraph==0
                    figure(fig(areaInd));
                else
                    figure(fig);
                end
            end
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
        if roving==0
            figure(fig);
        elseif roving==1
            if plotSameGraph==0
                figure(fig(areaInd));
            else
                figure(fig);
            end
        end
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
                        allData=[];
                        allSessionData=[];
                        allXData=[];
                        startEndTime=['_',num2str(periods(subPeriod)),'_to_',num2str(periods(subPeriod+1))];
                        psychoThresholdMatName=[area,appendText,'wholetrial_psyThreshold_dualThresholds'];
                        psychoThresholdMatFolder=fullfile(rootFolder,'PL',psychoType,animal,'psyThreshold_mat');
                        psychoThresholdMatPathname=fullfile(psychoThresholdMatFolder,psychoThresholdMatName);
                        loadText=['load ',psychoThresholdMatPathname,' sessionSorted2 threshold82lower threshold82higher'];
                        eval(loadText)
                        sessPsyInd=[];
                        for sessInd=1:length(sessions)
                            sessPsyInd=[sessPsyInd find(sessions(sessInd)==sessionSorted2)];
                        end
                        thresholdsNL=zeros(1,length(sessPsyInd));
                        thresholdsNH=zeros(1,length(sessPsyInd));
                        thresholdsPL=zeros(1,length(sessPsyInd));
                        thresholdsPH=zeros(1,length(sessPsyInd));
                        sessionTally=zeros(1,length(sessPsyInd));
                        allThresholdsNL=cell(1,length(sessPsyInd));
                        allThresholdsNH=cell(1,length(sessPsyInd));
                        allThresholdsPL=cell(1,length(sessPsyInd));
                        allThresholdsPH=cell(1,length(sessPsyInd));
                        Ncells=0;
                        for sessionInd=1:length(threshold82lower)
                            thresholdsPL(sessionInd)=sampleContrast-threshold82lower(sessionInd);%find difference between threshold and sample contrast
                            thresholdsPH(sessionInd)=threshold82higher(sessionInd)-sampleContrast;%find difference between threshold and sample contrast
                        end
                        psychoThresholds=[thresholdsPL;thresholdsPH];
                        isnanListPL=zeros(1,length(threshold82higher));
                        isnanListPH=zeros(1,length(threshold82lower));
                        for sessPsyInd=1:length(threshold82lower)
                            if isnan(threshold82lower(sessPsyInd))
                                thresholdsPL(sessPsyInd)=sampleContrast;%set to max if threshold cannot be obtained
                                isnanListPL(sessPsyInd)=1;
                            else
                                thresholdsPL(sessPsyInd)=sampleContrast-threshold82lower(sessPsyInd);%find difference between threshold and sample contrast
                            end
                            if isnan(threshold82higher(sessPsyInd))
                                thresholdsPH(sessPsyInd)=100-sampleContrast;%set to max if threshold cannot be obtained
                                isnanListPH(sessPsyInd)=1;
                            else
                                thresholdsPH(sessPsyInd)=threshold82higher(sessPsyInd)-sampleContrast;%find difference between threshold and sample contrast
                            end
                        end
                        allSessionData=[allSessionData;{sessionSorted2} {sessionSorted2}];
                        allXData=[allXData;{1:length(sessionSorted2)} {1:length(sessionSorted2)}];
                        if sglroc3IndividualChs==1
                            subFolder='new_vs_old_sglroc3acrosschannels';
                        elseif sglroc3IndividualChs==0
                            if useISI==0
                                subFolder='new_vs_old_sglrocmeanchannels';
                            elseif useISI==1
                                subFolder='new_ROC_useISI_meanchannels';
                            end
                        end
                        if excludeSuppressed
                            subFolder=[subFolder,'_excludeSuppressed'];
                        end
                        if normalize
                            subFolder=[subFolder,'_normalised'];
                        end
                        matname=['cumulative_ROCs_old_new_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10)];
                        pathname=fullfile(rootFolder,'PL','ROC_zero_one',animal,subFolder,matname);
                        loadText=['load ',pathname,' threshold82lower threshold82higher'];
                        eval(loadText)
                        maxThreshold=[sampleContrast 100-sampleContrast];%to draw filled markers representing uncalculable thresholds
                        for sessionInd=1:length(threshold82lower)
                            thresholdsNL(sessionInd)=sampleContrast-threshold82lower(sessionInd);%find difference between threshold and sample contrast
                            thresholdsNH(sessionInd)=threshold82higher(sessionInd)-sampleContrast;%find difference between threshold and sample contrast
                        end
                        neuroThresholds=[thresholdsNL;thresholdsNH];
                        isnanListNL=zeros(1,length(threshold82higher));
                        isnanListNH=zeros(1,length(threshold82lower));
                        for sessionInd=1:length(threshold82lower)
                            if isnan(threshold82lower(sessionInd))
                                thresholdsNL(sessionInd)=sampleContrast;%set to max if threshold cannot be obtained
                                isnanListNL(sessionInd)=1;
                            else
                                thresholdsNL(sessionInd)=sampleContrast-threshold82lower(sessionInd);%find difference between threshold and sample contrast
                            end
                            if isnan(threshold82higher(sessionInd))
                                thresholdsNH(sessionInd)=100-sampleContrast;%set to max if threshold cannot be obtained
                                isnanListNH(sessionInd)=1;
                            else
                                thresholdsNH(sessionInd)=threshold82higher(sessionInd)-sampleContrast;%find difference between threshold and sample contrast
                            end
                        end
                        allData=[allData;{thresholdsNL} {thresholdsNH} {thresholdsPL} {thresholdsPH}];
                        allisNanList=[{isnanListNL} {isnanListNH} {isnanListPL} {isnanListPH}];
                        allSessionData=[{sessions} {sessions} allSessionData];
                        allXData=[{1:length(sessions)} {1:length(sessions)} allXData];
                        startComparison=1;
                        if strcmp(psychoType,'psycho_param')
                            startComparison=3;
                        end
                        for comparison=startComparison:4
%                             temp=isnan(allData{1,comparison});
%                             allData{1,comparison}=allData{1,comparison}(~temp);
%                             allSessionData{1,comparison}=allSessionData{1,comparison}(~temp);
%                             allXData{1,comparison}=allXData{1,comparison}(~temp);
                            c=[(allSessionData{1,comparison})' (allData{1,comparison})'];%session in column 1, data in column 2
                            if size(c,1)>floor(length(sessions)/2)%only analyse if at least half of session have measurable threshold
                                [a b CI,stats]=corrcoef(c);
                                coefficientsCat(animalInd+2*(areaInd-1),comparison)=a(2);%1st row: B V4, 2nd: J V4, 3rd: B V1, 4th: J V1; 1st column: NL, 2nd: NH, 3rd: PL, 4th: PH
                                pCat(animalInd+2*(areaInd-1),comparison)=b(2);
                                statsCat{animalInd+2*(areaInd-1),comparison}=stats;
                                [a b]=corr(c,'type','Spearman');
                                if roving==0
                                    rho(animalInd+2*(areaInd-1),comparison)=a(2);
                                    pSpearman(animalInd+2*(areaInd-1),comparison)=b(2);
                                    dfs(animalInd+2*(areaInd-1),comparison)=size(c,1)-2;
                                elseif roving==1
                                    rho{sampCond}(animalInd+2*(areaInd-1),comparison)=a(2);
                                    pSpearman{sampCond}(animalInd+2*(areaInd-1),comparison)=b(2);
                                    dfs{sampCond}(animalInd+2*(areaInd-1),comparison)=size(c,1)-2;
                                end
                            end
                        end
                        if length((allData{1,3}))==length((allData{1,4}))
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
                        end
                        if plotPsychoOnly==0%neurometric data
                            subplotRemap=[1 2];
                            thresholdInterest=neuroThresholds;
                        elseif plotPsychoOnly==1%psychometric data
                            subplotRemap=[3 4];
                            thresholdInterest=psychoThresholds;
                        end
                        if roving==0
                            figure(fig);
                            subplot(length(areas),2,animalInd+2*(areaInd-1));
                        elseif roving==1
                            if plotSameGraph==0
                                figure(fig(areaInd));
                                subplot(length(sampleContrasts),2,animalInd+2*(sampCond-1));
                            else
                                figure(fig);
                                subplot(length(sampleContrasts),2,animalInd+2*(sampCond-1));
                            end
                        end
                        for subplotInd=subplotRemap(1):subplotRemap(2)
                            %                     subplot(2,2,subplotRemap(subplotInd));
                            if ~isempty(allData{1,subplotInd})
                                if plotSameGraph==0
                                    startPoint=0;
                                elseif plotSameGraph==1&&areaInd==2&&subplotInd==1%plot different stages of roving sessions adjacent to each other
                                    [startPoint]=get(gca,'xlim');
                                    startPoint=startPoint(2);
                                    plot([startPoint-0.5 startPoint-0.5],[0 100],'k--');
                                elseif plotSameGraph==1&&areaInd==1
                                    startPoint=0;
                                end
                                if subplotInd<=2
                                    if plotPsychoOnly==0
                                        for sessionNum=1:size(thresholdInterest,2)
                                            if allisNanList{subplotInd}(sessionNum)
                                            %if isnan(thresholdInterest(subplotInd,sessionNum))
                                                plot(sessionNum+startPoint,maxThreshold(subplotInd),'Marker','o','Color',markerCols(subplotInd),'LineStyle','none','MarkerFaceColor','none');hold on
                                            else
                                                plot(sessionNum+startPoint,thresholdInterest(subplotInd,sessionNum),'Marker','o','Color',markerCols(subplotInd),'LineStyle','none','MarkerFaceColor',markerCols(subplotInd));hold on
                                            end
                                        end
                                    end
                                else
                                end
                            end
                            if plotSameGraph==0
                                xlim([0 length(allData{1,subplotInd})+1]);
                            elseif roving==1
                                if strcmp(area,'v1_2_1')
                                    animalsXLim=[34 15];
                                elseif strcmp(area,'v1_2_2')
                                    animalsXLim=[50 40];
                                end
                                xlim([0 animalsXLim(animalInd)+1]);
                            end
                        end
                        if strcmp(animal,'blanco')&&strcmp(area,'v4_1')
                            subplot(2,2,1);
%                             plot(3,9.96,'Marker','o','Color','b','LineStyle','none','MarkerFaceColor','b')
%                             ylim([0 10]);
%                             set(gca,'YTick',[0:10],'YTickLabel',[{'0'} {'1'} {'2'} {'3'} {'4'} {'5'} {'6'} {'7'} {'8'} {''} {'23'}]);
%                             text(0,9,'//','FontSize',15)
                        end
                        if strcmp(animal,'jack')&&strcmp(area,'v1_1')
                            subplot(2,2,4);
%                             plot(1,30.0,'Marker','o','Color','b','LineStyle','none','MarkerFaceColor','none')
%                             ylim([0 30]);
%                             set(gca,'YTick',[0 5 10 15 20 25 30],'YTickLabel',[{'0'} {'5'} {'10'} {'15'} {'20'} {''} {'70'}]);
%                             text(0,25,'//','FontSize',15)
                        end
                        for subplotInd=1:length(subplotRemap)
                            if plotSameGraph==0
                                axis square
                            end
                            ylimits=get(gca,'YLim');
                            ylim([0 ylimits(2)]);
                            if animalInd==1
                                if roving==0
                                    ptext=areaTexts{areaInd};
                               elseif roving==1
                                    ptext=areaTexts{sampCond};
                                end
                                xlimits=get(gca,'XLim');
%                                 text('Position',[xlimits(1)-(xlimits(2)-xlimits(1))/4 (ylimits(2)-ylimits(1))/2],'FontSize',9,'String',ptext);
                            end
                            if roving==0
                                if areaInd==1&&sampCond==1
                                    title(animalTexts{animalInd});
                                    if length(areas)==1
                                        xlabel('session');
                                        ylabel('threshold');
                                    end
                                else
                                    xlabel('session');
                                    ylabel('threshold');
                                end
                            elseif roving==1
                                if animalInd==1&&sampCond==1
                                    xlabel('session');
                                    ylabel('threshold')
                                end
                                if sampCond==1
                                    title(animalTexts{animalInd});
                                end
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
                if roving==1
                    grandANOVATable(1:3,(sampCond*2)-1+6*(animalInd-1))=anovap{sampCond,animalInd};
                end
            end
        end
    end
end
if flankerlessOnly==1%plot just flankerless data
    figure(fig);
    subplot(3,2,1);
    xlim([0 34]);
    subplot(3,2,3);
    xlim([0 34]);
    subplot(3,2,5);
    xlim([0 34]);
    subplot(3,2,2);
    xlim([0 15]);
    subplot(3,2,4);
    xlim([0 15]);
    subplot(3,2,6);
    xlim([0 15]);
    subplot(3,2,6);
    text(1,3,'no flankers');
    text(12,3,'flankers');
    subplot(3,2,4);
    plot(3,8,'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r');
    plot(1,8,'Marker','o','MarkerEdgeColor','b','MarkerFaceColor','b');
    plot(1,4,'Marker','o','MarkerEdgeColor','b','MarkerFaceColor','none');
    plot(3,4,'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','none');
    text(3.5,7,'N_H')
    text(1.5,7,'P_H')
    text(1.5,3,'P_L')
    text(3.5,3,'N_L')
end
if roving==0
    rho2=rho';
    pSpearman2=pSpearman';
    dfs2=dfs';
    allCorrStats=[dfs2(:,1) rho2(:,1) pSpearman2(:,1) dfs2(:,3) rho2(:,3) pSpearman2(:,3);dfs2(:,2) rho2(:,2) pSpearman2(:,2) dfs2(:,4) rho2(:,4) pSpearman2(:,4)];
    allCorrStats=[allCorrStats(1,:);allCorrStats(3,:);allCorrStats(2,:);allCorrStats(4,:);allCorrStats(5,:);allCorrStats(7,:);allCorrStats(6,:);allCorrStats(8,:)];
elseif roving==1
    for sampCond=1:3
        rho2=rho{sampCond}';
        pSpearman2=pSpearman{sampCond}';
        dfs2=dfs{sampCond}';
        if length(areas)==2
            allCorrStats=[dfs2(:,1) rho2(:,1) pSpearman2(:,1) dfs2(:,3) rho2(:,3) pSpearman2(:,3);dfs2(:,2) rho2(:,2) pSpearman2(:,2) dfs2(:,4) rho2(:,4) pSpearman2(:,4)];
            allCorrStats=[allCorrStats(1,:);allCorrStats(3,:);allCorrStats(2,:);allCorrStats(4,:);allCorrStats(5,:);allCorrStats(7,:);allCorrStats(6,:);allCorrStats(8,:)];
        elseif length(areas)==1
            allCorrStats=[dfs2(:,1) rho2(:,1) pSpearman2(:,1) dfs2(:,2) rho2(:,2) pSpearman2(:,2)];
            allCorrStats=[allCorrStats(1,:);allCorrStats(3,:);allCorrStats(2,:);allCorrStats(4,:)];
        end
        allCorrSamplesStats{sampCond}=allCorrStats;
    end
    subplot(3,2,6);
    text(1,3,'no flankers');
    text(34,3,'flankers');
    subplot(3,2,4);
    plot(10,8,'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r');
    plot(2,8,'Marker','o','MarkerEdgeColor','b','MarkerFaceColor','b');
    plot(2,4,'Marker','o','MarkerEdgeColor','b','MarkerFaceColor','none');
    plot(10,4,'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','none');
    text(11,7,'N_H')
    text(3,7,'P_H')
    text(3,3,'P_L')
    text(11,3,'N_L')
    subplot(3,2,1);xlim([0 34.5]);
    subplot(3,2,3);xlim([0 34.5]);
    subplot(3,2,5);xlim([0 34.5]);
    subplot(3,2,2);xlim([0 15.5]);
    subplot(3,2,4);xlim([0 15.5]);
    subplot(3,2,6);xlim([0 15.5]);
end
if plotPsychoOnly==1
    tableRanksum=[allstatsranksum(:,1)';allpranksum(:,1)';allstatsranksum(:,3)';allpranksum(:,3)';allstatsranksum(:,2)';allpranksum(:,2)';allstatsranksum(:,4)';allpranksum(:,4)'];
end

if roving==0
    saveCoefImageName=[areaTexts{1},'_',areaTexts{2},'_',analysisType,'_coefs',startEndTime,'_populationPROBMAT'];
elseif roving==1
    saveCoefImageName=[areaTexts{1},'_',analysisType,'_coefs',startEndTime,'_populationPROBMAT'];
end
saveCoefImageFolder=fullfile(rootFolder,'PL',analysisType);
saveCoefImagePath=fullfile(saveCoefImageFolder,saveCoefImageName);
printtext=sprintf('print -dpng %s.png',saveCoefImagePath);
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperPositionMode', 'auto');
eval(printtext);
