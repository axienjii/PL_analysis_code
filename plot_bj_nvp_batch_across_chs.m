function plot_bj_nvp_batch_across_chs(animal,area,analysisType)

%written 29/03/13, modified from plot_bj_crf_or_roc_batch_across_chs.
%Reads neuro- and psychometric threshold values for all conditions, to calculate averages and
%examine correlation between threshold and session.
analysisType='NVP';
animalTexts=[{'subject B'} {'subject J'}];
animals=[{'blanco'} {'jack'}];
roving=0;
if roving==0
    areaTexts=[{'V4'} {'V1'}];
    areas=[{'v4_1'} {'v1_1'}];
elseif roving==1
    areaTexts={'V1 roving data'};
    areas={'v1_2'};
end
fig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.55, 0.8]);
set(fig,'PaperUnits','centimeters','PaperType','A4','PaperOrientation', 'portrait', 'PaperPosition', [0.63452 0.63452 21/8*5.5 28.41]);
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        channels = main_channels(animal,area);
        sessions = main_raw_sessions_final(animal,area,[],0);
        includeMatch=[];
        excludeSessions=[26 50 306 312 316 322:328 342];
        for includeInd=1:length(sessions)
            if sum(excludeSessions==sessions(includeInd))==0
                includeMatch=[includeMatch includeInd];
            end
        end
        sessions=sessions(includeMatch);
        [sampleContrasts testContrasts]=area_metadata(area);
        test_epochs={0 529 529*2 529*3};durSpon=150;
        
        psychoname=['psycho_constants_',area];
        psychoPathname=fullfile('F:','PL','psycho_data',animal,psychoname);
        subplotTitleText={'neurometric threshold lower contrast' 'neurometric threshold higher contrast' 'psychometric threshold lower contrast' 'psychometric threshold higher contrast'};
        
        includeOnlySigC50Chs=1;
        if includeOnlySigC50Chs==1
            includeSigOnlyText='_onlySigC50ChsIncluded';
        else
            includeSigOnlyText=[];
        end
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
        subplotInd=subplot(length(areas),2,animalInd+2*(areaInd-1));
        cellCategory=zeros(length(channels),2)-2;
        cellCategory(:,1)=channels;
        for sampCond=1:length(sampleContrasts)
            sampleContrast=sampleContrasts(sampCond);
            testContrast=testContrasts(sampCond,:);
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
                        psychoThresholdMatFolder=fullfile('F:','PL','psycho',animal,'psyThreshold_mat');
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
                            thMatname=[num2str(chNum),appendText,startEndTime,'_nvpThreshold'];
                            thMatFolder=fullfile('F:','PL',analysisType,animal,'nvpThreshold_mat');
                            thMatPathname=fullfile(thMatFolder,thMatname);
                            loadText=['load ',thMatPathname,' sessionSorted1 threshold82lower threshold82higher'];
                            eval(loadText)
                            Ncells=Ncells+1;
                            includeMatch=[];
                            for includeInd=1:length(sessionSorted1)
                                if sum(excludeSessions==sessionSorted1(includeInd))==0
                                    includeMatch=[includeMatch includeInd];
                                end
                            end
                            sessionSorted1=sessionSorted1(:,includeMatch);
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
                        for comparison=1:4
                            temp=isnan(allData{1,comparison});
                            allData{1,comparison}=allData{1,comparison}(~temp);
                            allSessionData{1,comparison}=allSessionData{1,comparison}(~temp);
                            allXData{1,comparison}=allXData{1,comparison}(~temp);
                            c=[(allSessionData{1,comparison})' (allData{1,comparison})'];
                            [a b]=corrcoef(c);%slope
                            coefficientsCat(animalInd+2*(areaInd-1),comparison)=a(2);%1st row: B V4, 2nd: B V1, 3rd: J V4, 4th: J V1; 1st column: NL, 2nd: NH, 3rd: PL, 4th: PH
                            pCat(animalInd+2*(areaInd-1),comparison)=b(2);
                        end
                        subplotRemap=[1 2 3 4];
                        for subplotInd=1:length(subplotRemap)
                            %                     subplot(2,2,subplotRemap(subplotInd));
                            if ~isempty(allData{1,subplotInd})
                                if subplotInd<=2
                                    plot(allXData{1,subplotInd},allData{1,subplotInd},'Marker','o','Color',markerCols(subplotInd),'LineStyle','none','MarkerFaceColor','none');hold on
                                else
                                    plot(allXData{1,subplotInd},allData{1,subplotInd},'Marker','o','Color',markerCols(subplotInd),'LineStyle','none','MarkerFaceColor',markerCols(subplotInd));hold on
                                end
                            end
                            xlim([0 length(allData{1,subplotInd})+1]);
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
                            if animalInd==1
                                ptext=areaTexts{areaInd};
                                xlimits=get(gca,'XLim');
                                ylimits=get(gca,'YLim');
                                text('Position',[xlimits(1)-(xlimits(2)-xlimits(1))/5 (ylimits(2)-ylimits(1))/2],'FontSize',9,'String',ptext);
                                ylabel('threshold');
                            end
                            if areaInd==1
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
        end
%         saveCoefImageName=['all_',animal,'_',area,'_',analysisType,'_coefs',startEndTime,includeSigOnlyText];
%         subFolder=[analysisType,'_coef_images'];
%         saveCoefImageFolder=fullfile('F:','PL',analysisType,animal,subFolder);
%         saveCoefImagePath=fullfile(saveCoefImageFolder,saveCoefImageName);
%         printtext=sprintf('print -dpng %s.png',saveCoefImagePath);
%         set(gcf, 'PaperOrientation', 'landscape');
%         set(gcf, 'PaperPositionMode', 'auto');
%         eval(printtext);
        anovaVals=[thresholdsNL thresholdsNH thresholdsPL thresholdsPH];
        nvpFactor=[zeros(1,length(thresholdsNL))+1 zeros(1,length(thresholdsNH))+1 zeros(1,length(thresholdsPL))+2 zeros(1,length(thresholdsPH))+2];
        highlowFactor=[zeros(1,length(thresholdsNL))+1 zeros(1,length(thresholdsNH))+2 zeros(1,length(thresholdsPL))+1 zeros(1,length(thresholdsPH))+2];
%         sessionFactor=[zeros(1,(floor(length(thresholdsNL)/2)))+1 zeros(1,(ceil(length(thresholdsNL)/2)))+2 zeros(1,(floor(length(thresholdsNH)/2)))+1 zeros(1,(ceil(length(thresholdsNH)/2)))+2 zeros(1,(floor(length(thresholdsPL)/2)))+1 zeros(1,(ceil(length(thresholdsPL)/2)))+2 zeros(1,(floor(length(thresholdsPH)/2)))+1 zeros(1,(ceil(length(thresholdsPH)/2)))+2];
%         sessionFactor=[zeros(1,(floor(length(thresholdsNL)/3)))+1 zeros(1,(floor(length(thresholdsNL)/3)))+2 zeros(1,length(thresholdsNL)-(2*floor(length(thresholdsNL)/3)))+3 zeros(1,(floor(length(thresholdsNH)/3)))+1 zeros(1,(floor(length(thresholdsNH)/3)))+2 zeros(1,length(thresholdsNH)-(2*floor(length(thresholdsNH)/3)))+3 zeros(1,(floor(length(thresholdsPL)/3)))+1 zeros(1,(floor(length(thresholdsPL)/3)))+2 zeros(1,length(thresholdsPL)-(2*floor(length(thresholdsPL)/3)))+3 zeros(1,(floor(length(thresholdsPH)/3)))+1 zeros(1,(floor(length(thresholdsPH)/3)))+2 zeros(1,length(thresholdsPH)-(2*floor(length(thresholdsPH)/3)))+3];
%         [p,table,stats]=anovan(anovaVals,{nvpFactor highlowFactor sessionFactor},'model','full');anovap{animalInd+2*(areaInd-1)}=p;anovatable{animalInd+2*(areaInd-1)}=table;anovastats{animalInd+2*(areaInd-1)}=stats;
        [p,table,stats]=anovan(anovaVals,{nvpFactor highlowFactor},'model','full');anovap{animalInd+2*(areaInd-1)}=p;anovatable{animalInd+2*(areaInd-1)}=table;anovastats{animalInd+2*(areaInd-1)}=stats;
        [c,m,h,nms]=multcompare(stats)
         %%1st row: B V4, 2nd: B V1, 3rd: J V4, 4th: J V1; 1st column: NL vs NH, 2nd: PL vs PH, 3rd: NL vs PL, 4th: NH vs PH
        [h,p,ci,stats]=ttest(thresholdsNL,thresholdsNH);allh(1,animalInd+2*(areaInd-1))=h;allp(1,animalInd+2*(areaInd-1))=p;allci{1,animalInd+2*(areaInd-1)}=ci;allstats{1,animalInd+2*(areaInd-1)}=stats;
        [h,p,ci,stats]=ttest(thresholdsPL,thresholdsPH);allh(2,animalInd+2*(areaInd-1))=h;allp(2,animalInd+2*(areaInd-1))=p;allci{2,animalInd+2*(areaInd-1)}=ci;allstats{2,animalInd+2*(areaInd-1)}=stats;
        [h,p,ci,stats]=ttest(thresholdsNL,thresholdsPL);allh(3,animalInd+2*(areaInd-1))=h;allp(3,animalInd+2*(areaInd-1))=p;allci{3,animalInd+2*(areaInd-1)}=ci;allstats{3,animalInd+2*(areaInd-1)}=stats;
        [h,p,ci,stats]=ttest(thresholdsNH,thresholdsPH);allh(4,animalInd+2*(areaInd-1))=h;allp(4,animalInd+2*(areaInd-1))=p;allci{4,animalInd+2*(areaInd-1)}=ci;allstats{4,animalInd+2*(areaInd-1)}=stats;
    end
end
if roving==0
    saveCoefImageName=[areaTexts{1},'_',areaTexts{2},'_',analysisType,'_coefs',startEndTime];
elseif roving==1
    saveCoefImageName=[areaTexts{1},'_',analysisType,'_coefs',startEndTime];
end
saveCoefImageFolder=fullfile('F:','PL',analysisType);
saveCoefImagePath=fullfile(saveCoefImageFolder,saveCoefImageName);
printtext=sprintf('print -dpng %s.png',saveCoefImagePath);
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperPositionMode', 'auto');
eval(printtext);
