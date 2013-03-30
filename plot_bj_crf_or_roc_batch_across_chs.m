function plot_bj_crf_or_roc_batch_across_chs(animal,area,analysisType)

%written 03/11/11, modified from plot_blanco_V1_roc_batch2_across_chs. 
%Batch file to read the ROC values for the 3 contrasts immediately 
%flanking the sample contrast, plots ROC value across sessions, for each of 
%these 6 contrast conditions. 
%Also reads all 14 ROC values for 14 conditions, to calculate averages and
%examine correlation between ROC and session, c50 and session, and ROC
%slope and psychometric slope.

channels = main_channels(animal,area);
sessions = main_raw_sessions_final(animal,area);
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
subplotTitleText={'slope vs session' 'C50 vs session' 'neuro vs psycho' 'diff C50 vs session' 'minRate vs session' 'maxRate vs session'};

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
markerCols='mkc';

fig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
set(fig,'PaperUnits','centimeters','PaperType','A4','PaperOrientation', 'portrait', 'PaperPosition', [0.63452 0.63452 21 28.41]);
cellCategory=zeros(length(channels),2)-2;
cellCategory(:,1)=channels;
for sampCond=1:length(sampleContrasts)   
    sampleContrast=sampleContrasts(sampCond);
    testContrast=testContrasts(sampCond,:);
    for epoch=1:size(test_epochs,2)
        if strcmp(analysisType,'CRF')||strcmp(analysisType,'ROC')&&epoch==4
            if epoch==1
                periods=[-durSpon 0];
            else
                periods=[test_epochs{epoch-1} test_epochs{epoch}(1)];
            end
            for subPeriod=1:length(periods)-1
                startEndTime=['_',num2str(periods(subPeriod)),'_to_',num2str(periods(subPeriod+1))];                
                slopes=zeros(3,length(sessions));
                c50s=zeros(3,length(sessions));
                diffc50s=zeros(3,length(sessions));
                minRates=zeros(3,length(sessions));
                maxRates=zeros(3,length(sessions));
                sessionTally=zeros(3,length(sessions));
                allSlopes=cell(3,length(sessions));
                allC50s=cell(3,length(sessions));
                allDiffc50s=cell(3,length(sessions));
                allMinRates=cell(3,length(sessions));
                allMaxRates=cell(3,length(sessions));
                Ncells=0;
                for i=1:length(channels)
                    chNum=channels(i);
                    if strcmp(analysisType,'ROC')
                        coefsFileName=[num2str(chNum),'_',num2str(sampleContrast),startEndTime,'_',analysisType,'_coefs_',area,'_goodSSE_no_outliers_sl3_C503'];%,'_goodSSE_no_outliers_sl3_C503','.mat'];
                    else
                        coefsFileName=[num2str(chNum),'_',num2str(sampleContrast),startEndTime,'_',analysisType,'_coefs_',area];%,'_goodSSE_no_outliers_sl3_C503','.mat'];
                    end
                    subFolder=[analysisType,'_coef_mat'];
                    coefsPathName=fullfile('F:','PL',analysisType,animal,subFolder,coefsFileName);
                    loadText=['load ',coefsPathName];
                    eval(loadText);
                    if includeOnlySigC50Chs==1&&coefficients(2,2)<0.05||includeOnlySigC50Chs==0
                        Ncells=Ncells+1;
                        includeMatch=[];
                        for includeInd=1:length(sessionSorted2)
                            if sum(excludeSessions==sessionSorted2(includeInd))==0
                                includeMatch=[includeMatch includeInd];
                            end
                        end
                        sessionSorted2=sessionSorted2(includeMatch,:);
                        slopePsycho=slopePsycho(includeMatch,:);
                        if mean(c50(1:5))>sampleContrast+margin%cell belongs to group where average c50 across first 5 sessions is more than sample contrast plus a certain value
                            cellCategory(i,2)=1;
                        elseif mean(c50(1:5))<sampleContrast-margin%average c50 across first 5 sessions is less than sample contrast minus a certain value
                            cellCategory(i,2)=2;
                        elseif mean(c50(1:5))<sampleContrast+margin&&mean(c50(1:5))>30-margin%average c50 across first 5 sessions is less than sample contrast minus a certain value
                            cellCategory(i,2)=3;
                        end
                        for sessionInd=1:length(sessionSorted1)
                            colInd=find(sessionSorted1(sessionInd)==sessionSorted2);
                            sessionTally(cellCategory(i,2),colInd)=sessionTally(cellCategory(i,2),colInd)+1;
                            slopes(cellCategory(i,2),colInd)=slopes(cellCategory(i,2),colInd)+slopeNeuro(sessionInd);
                            c50s(cellCategory(i,2),colInd)=c50s(cellCategory(i,2),colInd)+c50(sessionInd);
                            diffc50s(cellCategory(i,2),colInd)=diffc50s(cellCategory(i,2),colInd)+diffc50(sessionInd);
                            minRates(cellCategory(i,2),colInd)=minRates(cellCategory(i,2),colInd)+minRate(sessionInd);
                            maxRates(cellCategory(i,2),colInd)=maxRates(cellCategory(i,2),colInd)+maxRate(sessionInd);
                            
                            allSlopes{cellCategory(i,2),colInd}=[allSlopes{cellCategory(i,2),colInd} slopeNeuro(sessionInd)];
                            allC50s{cellCategory(i,2),colInd}=[allC50s{cellCategory(i,2),colInd} c50(sessionInd)];
                            allDiffc50s{cellCategory(i,2),colInd}=[allDiffc50s{cellCategory(i,2),colInd} diffc50(sessionInd)];
                            allMinRates{cellCategory(i,2),colInd}=[allMinRates{cellCategory(i,2),colInd} minRate(sessionInd)];
                            allMaxRates{cellCategory(i,2),colInd}=[allMaxRates{cellCategory(i,2),colInd} maxRate(sessionInd)];
                        end
                    end
                end
                sum(cellCategory(:,2)==1)
                sum(cellCategory(:,2)==2)
                sum(cellCategory(:,2)==3)
                slopes=slopes./sessionTally;
                c50s=c50s./sessionTally;
                diffc50s=diffc50s./sessionTally;
                minRates=minRates./sessionTally;
                maxRates=maxRates./sessionTally;
                allData=[];
                allSessionData=[];
                allXData=[];
                for cellCat=1:3
                    allData=[allData;{slopes(cellCat,:)} {c50s(cellCat,:)} {slopes(cellCat,:)} {diffc50s(cellCat,:)} {minRates(cellCat,:)} {maxRates(cellCat,:)}];
                    allSessionData=[allSessionData;{sessions} {sessions} {slopePsycho'} {sessions} {sessions} {sessions}];
                    allXData=[allXData;{1:length(sessions)} {1:length(sessions)} {slopePsycho'} {1:length(sessions)} {1:length(sessions)} {1:length(sessions)}];
                end
                for cellCat=1:3                    
                    for comparison=1:6
                        temp=isnan(allData{cellCat,comparison});
                        allData{cellCat,comparison}=allData{cellCat,comparison}(~temp);
                        allSessionData{cellCat,comparison}=allSessionData{cellCat,comparison}(~temp);
                        allXData{cellCat,comparison}=allXData{cellCat,comparison}(~temp);
                        c=[(allSessionData{cellCat,comparison})' (allData{cellCat,comparison})'];
                        [a b]=corrcoef(c);%slope
                        coefficientsCat(cellCat,comparison)=a(2);
                        pCat(cellCat,comparison)=b(2);
                    end
                    subplotRemap=[1 4 2 5 3 6];
                    for subplotInd=1:6
                        subplot(2,3,subplotRemap(subplotInd));
                        if ~isempty(allData{cellCat,subplotInd})
                            plot(allXData{cellCat,subplotInd},allData{cellCat,subplotInd},'Marker','o','Color',markerCols(cellCat),'LineStyle','none');hold on
                        end
                    end
                end           
                for cellCat=1:3  
                    for subplotInd=1:6
                        subplot(2,3,subplotRemap(subplotInd));
                        if subplotRemap(subplotInd)==1
                            ptext=sprintf('N = %s channels',num2str(Ncells));
                            yLimVals=get(gca,'YLim');
                            text('Position',[-2 yLimVals(2)+0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
                        end
                        if pCat(cellCat,subplotInd)<0.05
                            colText='g';
                        else
                            colText=markerCols(cellCat);
                        end
                        ptext=sprintf('r= %f  p= %f',coefficientsCat(cellCat,subplotInd),pCat(cellCat,subplotInd));
                        yLimVals=get(gca,'YLim');
                        xLimVals=get(gca,'XLim');
                        remapText=[3 2 1];
                        text('Position',[xLimVals(1)+(xLimVals(2)-xLimVals(1))/2 yLimVals(1)+0.05*remapText(cellCat)*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext,'Color',colText);
                        title(subplotTitleText{subplotInd});
                    end
                end
            end
        end
    end
end
saveCoefImageName=['all_',animal,'_',area,'_',analysisType,'_coefs_C50_categories',startEndTime,includeSigOnlyText];
subFolder=[analysisType,'_coef_images'];
saveCoefImageFolder=fullfile('F:','PL',analysisType,animal,subFolder);
saveCoefImagePath=fullfile(saveCoefImageFolder,saveCoefImageName);
printtext=sprintf('print -dpng %s.png',saveCoefImagePath);
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperPositionMode', 'auto');
eval(printtext);
% if x==1
%     all_c50=[];%c50 list
%     c50=[];
%     sessions=[];
%     all_sessions=[];
%     cellCategory=zeros(1,length(channels))-2;
%     X0=[30 2];
%     options = optimset('Display','off','MaxFunEvals',10^6,'MaxIter',10^6,'TolFun',1.0E-6,'TolX',1.0E-6);
%     for i=1:length(channels)
%         ind=find(all_roc14(:,1)==channels(i));
%         if length(ind)>4
%             for c50count=1:length(ind)
%                 perfvals=all_roc14(ind(c50count),3:length(testContrast)+2);%neurometric performance for each session
%                 [X,fval]=fminsearch('weib_sim_min_max',X0,options,testContrast,perfvals);
%                 yvals=max(perfvals)-(max(perfvals)-min(perfvals)).*exp(-((testContrast./X(1)).^X(2)));
%                 [coeffit pfit]=corrcoef([yvals' perfvals']);
%                 if coeffit>0.5
%                     if X(1)>10&&X(1)<80
%                         all_roc14(ind(c50count),17)=X(2);%slope for each channel
%                         c50=[c50 X(1)];%store the c50s for that channel
%                         all_roc14(ind(c50count),18)=X(1);
%                         sessions=[sessions all_roc14(ind(c50count),2)];
%                     end
%                 end
%             end
%             all_c50=[all_c50 c50];
%             all_sessions=[all_sessions sessions];
%             if isempty(find(c50<10))&&isempty(find(c50>80))&&length(c50)>5
%                 if mean(c50(1:5))>sampleContrast+margin%cell belongs to group where average c50 across first 5 sessions is more than sample contrast plus a certain value
%                     cellCategory(i)=1;
%                 elseif mean(c50(1:5))<sampleContrast-margin%average c50 across first 5 sessions is less than sample contrast minus a certain value
%                     cellCategory(i)=-1;
%                 elseif mean(c50(1:5))<sampleContrast+margin&&mean(c50(1:5))>30-margin%average c50 across first 5 sessions is less than sample contrast minus a certain value
%                     cellCategory(i)=0;
%                 end
%                 b=[sessions' c50'];
%                 [coefficientsC1 pC1]=corrcoef(b);%c50
%                 channels_p(i,1)=channels(i);
%                 channels_p(i,2)=coefficientsC1(2);
%                 channels_p(i,3)=pC1(2);
%             end
%         end
%     end
%     b=[all_sessions' all_c50'];
%     [allcoefficientsC0 allpC0]=corrcoef(b);%c50
%     for i=1:length(sessionNums)
%         ind=find(sessions==sessionNums(i));
%         meanc50(i)=mean(all_c50(ind));
%     end
%     figure
%     plot(sessionNums,meanc50,'xk');
%     aboveC=[];%c50 list
%     belowC=[];
%     aroundC=[];
%     aboveSe=[];%session number
%     belowSe=[];
%     aroundSe=[];
%     aboveSl=[];%slope
%     belowSl=[];
%     aroundSl=[];
%     for i=1:length(sessionNums)
%         ind=find(all_roc14(:,2)==sessionNums(i));
%         for indCount=1:length(ind)
%             cellNum=all_roc14(ind(indCount),1);
%             cellInd=find(channels==cellNum);
%             if cellCategory(cellInd)==1&&all_roc14(ind(indCount),18)>10&&all_roc14(ind(indCount),18)<80
%                 aboveC=[aboveC all_roc14(ind(indCount),18)];
%                 aboveSl=[aboveSl all_roc14(ind(indCount),17)];
%                 aboveSe=[aboveSe all_roc14(ind(indCount),2)];
%             elseif cellCategory(cellInd)==0&&all_roc14(ind(indCount),18)>10&&all_roc14(ind(indCount),18)<80
%                 aroundC=[aroundC all_roc14(ind(indCount),18)];
%                 aroundSl=[aroundSl all_roc14(ind(indCount),17)];
%                 aroundSe=[aroundSe all_roc14(ind(indCount),2)];
%             elseif cellCategory(cellInd)==-1&&all_roc14(ind(indCount),18)>10&&all_roc14(ind(indCount),18)<80
%                 belowC=[belowC all_roc14(ind(indCount),18)];
%                 belowSl=[belowSl all_roc14(ind(indCount),17)];
%                 belowSe=[belowSe all_roc14(ind(indCount),2)];
%             end
%         end
%     end
%     
%     if length(find(cellCategory==1))>1%high initial c50 average value
%         figure
%         for i=1:length(sessionNums)
%             index=find(aboveSe==sessionNums(i));
%             meanC(i)=mean(aboveC(index));
%         end
%         b=[aboveSe' aboveC'];
%         [allcoefficientsC1 allpC1]=corrcoef(b);%c50 above margin
%         b=[sessionNums' meanC'];
%         [allcoefficientsC1m allpC1m]=corrcoef(b);%c50 above margin
%         if allpC1m(2)<0.05
%             setFill=markerCol;
%         else
%             setFill='none';
%         end
%         xvals=1:length(sessionNums);
%         plot(xvals,meanC,markerSh,'MarkerFaceColor',setFill,'LineStyle','none','MarkerEdgeColor',markerCol);hold on
%         brob = robustfit(xvals,meanC);
%         plot(xvals,brob(1)+brob(2)*xvals,'k','LineWidth',2);hold on
%         titletext=['channels with high initial c50 (> ',num2str(sampleContrast+margin),'%), p=',num2str(allpC1m(2))];
%         title(titletext);
%         ylabel('mean c50');
%         xlabel('session number');
%         plot([1 length(sessionNums)],[sampleContrast sampleContrast],'--k')
%         yLimVals=get(gca,'YLim');
%         ind=find(yLimVals==sampleContrast);
%         if ind==1            
%             ylim([yLimVals(1)-5,yLimVals(2)]);
%         elseif ind==2
%             ylim([yLimVals(1),yLimVals(2)+5]);
%         end
%         if allpC1(2)<0.05
%             setFill=markerCol;
%         else
%             setFill='none';
%         end
%         figName=[folder,addText,'\roc_c50_ch_above_',num2str(sampleContrast+margin),'mean'];
%         printtext=['print -dtiff ',figName];
%         eval(printtext);
%         hgsave(figName);
%         figure
%         xvals=1:length(sessionNums);
%         for i=1:length(sessionNums)
%             ind=find(aboveSe==sessionNums(i));
%             aboveSe2(ind)=i;
%         end
%         % for i=1:length(sessionNums)
%         %     ind=find(aboveSe==sessionNums(i));
%         %     if ~isempty(ind)
%         %         yvals=aboveC(ind);
%         %         plot(i,yvals,'s',...
%         %             'MarkerFaceColor',[currcol(floor(i*2.5),1) currcol(floor(i*2.5),2) currcol(floor(i*2.5),3)],...
%         %             'MarkerEdgeColor',[currcol(floor(i*2.5),1) currcol(floor(i*2.5),2) currcol(floor(i*2.5),3)]);hold on
%         %     end
%         % end
%         plot(aboveSe2,aboveC,markerSh,'MarkerFaceColor',setFill,'LineStyle','none','MarkerEdgeColor',markerCol);hold on
%         brob = robustfit(aboveSe2,aboveC);
%         plot(xvals,brob(1)+brob(2)*xvals,'k','LineWidth',4);hold on
%         titletext=['channels with high initial c50 (> ',num2str(sampleContrast+margin),'%), p=',num2str(allpC1(2))];
%         title(titletext);
%         ylabel('c50');
%         xlabel('session number');
%         set(gca,'yLim',[20 60])
%         % plot([1 23],[30 30],'--k')
%         % ylim([25 60])
%         plot([1 length(sessionNums)],[sampleContrast sampleContrast],'--k')
%         yLimVals=get(gca,'YLim');
%         ind=find(yLimVals==sampleContrast);
%         if ind==1            
%             ylim([yLimVals(1)-5,yLimVals(2)]);
%         elseif ind==2
%             ylim([yLimVals(1),yLimVals(2)+5]);
%         end
%         % set(gca,'XTick',1:17)
%         % set(gca,'XTickLabel',24:40)
% %         set(gca,'XTick',1:5:17)
% %         set(gca,'XTickLabel',24:5:40)
% %         set(gca,'XLim',[1 13])
%         figName=[folder,addText,'\roc_c50_ch_above_',num2str(sampleContrast+margin)];
%         printtext=['print -dtiff ',figName];
%         eval(printtext);
%         hgsave(figName);
%     end
%     if length(find(cellCategory==0))>1%intermediate initial c50 average value
%         figure
%         for i=1:length(sessionNums)
%             index=find(aroundSe==sessionNums(i));
%             meanC(i)=mean(aroundC(index));
%         end
%         nanind=~isnan(meanC);
%         meanC=meanC(nanind);
%         xvals=sessionNums(nanind);
%         xvals=1:length(xvals);
%         b=[aroundSe' aroundC'];
%         [allcoefficientsC2 allpC2]=corrcoef(b);%c50 around margin
%         b=[xvals' meanC'];
%         [allcoefficientsC2m allpC2m]=corrcoef(b);%c50 around margin
%         if allpC2m(2)<0.05
%             setFill=markerCol;
%         else
%             setFill='none';
%         end
%         plot(xvals,meanC,markerSh,'MarkerFaceColor',setFill,'LineStyle','none','MarkerEdgeColor',markerCol);hold on
%         brob = robustfit(xvals,meanC);
%         plot(xvals,brob(1)+brob(2)*xvals,'k','LineWidth',2);hold on
%         titletext=['channels with intermediate initial c50 (>',num2str(sampleContrast-margin),'% & <',num2str(sampleContrast+margin),'%), p=',num2str(allpC2m(2))];
%         title(titletext);
%         ylabel('mean c50');
%         xlabel('session number');
%         plot([1 length(sessionNums)],[sampleContrast sampleContrast],'--k')
%         yLimVals=get(gca,'YLim');
%         ind=find(yLimVals==sampleContrast);
%         if ind==1            
%             ylim([yLimVals(1)-5,yLimVals(2)]);
%         elseif ind==2
%             ylim([yLimVals(1),yLimVals(2)+5]);
%         end
%         if allpC2(2)<0.05
%             setFill=markerCol;
%         else
%             setFill='none';
%         end
%         figName=[folder,addText,'\roc_c50_ch_',num2str(sampleContrast-margin),'_to_',num2str(sampleContrast+margin),'mean'];
%         printtext=['print -dtiff ',figName];
%         eval(printtext);
%         hgsave(figName);
%         figure
%         xvals=1:length(sessionNums);
%         for i=1:length(sessionNums)
%             ind=find(aroundSe==sessionNums(i));
%             aroundSe2(ind)=i;
%         end
%         %     for i=1:length(sessionNums)
%         %         ind=find(aroundSe==sessionNums(i));
%         %         yvals=aroundC(ind);
%         %         plot(i,yvals,'s',...
%         %             'MarkerFaceColor',[currcol(floor(i*2.5),1) currcol(floor(i*2.5),2) currcol(floor(i*2.5),3)],...
%         %             'MarkerEdgeColor',[currcol(floor(i*2.5),1) currcol(floor(i*2.5),2) currcol(floor(i*2.5),3)]);hold on
%         %     end
%         plot(aroundSe2,aroundC,markerSh,'MarkerFaceColor',setFill,'LineStyle','none','MarkerEdgeColor',markerCol);hold on
%         brob = robustfit(aroundSe2,aroundC);
%         plot(xvals,brob(1)+brob(2)*xvals,'k','LineWidth',2);hold on
%         titletext=['channels with intermediate initial c50 (>',num2str(sampleContrast-margin),'% & <',num2str(sampleContrast+margin),'%), p=',num2str(allpC2(2))];
%         title(titletext);
%         ylabel('c50');
%         xlabel('session number');
%         plot([1 length(sessionNums)],[sampleContrast sampleContrast],'--k')
%         yLimVals=get(gca,'YLim');
%         ind=find(yLimVals==sampleContrast);
%         if ind==1            
%             ylim([yLimVals(1)-5,yLimVals(2)]);
%         elseif ind==2
%             ylim([yLimVals(1),yLimVals(2)+5]);
%         end
%         figName=[folder,addText,'\roc_c50_ch_',num2str(sampleContrast-margin),'_to_',num2str(sampleContrast+margin)];
%         printtext=['print -dtiff ',figName];
%         eval(printtext);
%         hgsave(figName);
%     end
%     if length(find(cellCategory==-1))>1%low initial c50 average value
%         figure
%         for i=1:length(sessionNums)
%             index=find(belowSe==sessionNums(i));
%             meanC(i)=mean(belowC(index));
%         end
%         nanind=~isnan(meanC);
%         meanC=meanC(nanind);
%         xvals=sessionNums(nanind);
%         xvals=1:length(xvals);
%         b=[belowSe' belowC'];
%         [allcoefficientsC3 allpC3]=corrcoef(b);%c50 below margin
%         b=[xvals' meanC'];
%         [allcoefficientsC3m allpC3m]=corrcoef(b);%c50 below margin
%         if allpC3m(2)<0.05
%             setFill=markerCol;
%         else
%             setFill='none';
%         end
%         plot(xvals,meanC,markerSh,'MarkerFaceColor',setFill,'LineStyle','none','MarkerEdgeColor',markerCol);hold on
%         brob = robustfit(xvals,meanC);
%         plot(xvals,brob(1)+brob(2)*xvals,'k','LineWidth',2);hold on
%         titletext=['channels with low initial c50 (<',num2str(sampleContrast-margin),'%), p=',num2str(allpC3m(2))];
%         title(titletext);
%         ylabel('mean c50');
%         xlabel('session number');
%         plot([1 length(sessionNums)],[sampleContrast sampleContrast],'--k')
%         yLimVals=get(gca,'YLim');
%         yLimVals=get(gca,'YLim');
%         ind=find(yLimVals==sampleContrast);
%         if ind==1            
%             ylim([yLimVals(1)-5,yLimVals(2)]);
%         elseif ind==2
%             ylim([yLimVals(1),yLimVals(2)+5]);
%         end
%         if allpC2(2)<0.05
%             setFill=markerCol;
%         else
%             setFill='none';
%         end
%         figName=[folder,addText,'\roc_c50_ch_below_',num2str(sampleContrast-margin),'mean'];
%         printtext=['print -dtiff ',figName];
%         eval(printtext);
%         hgsave(figName);
%         figure
%         xvals=1:length(sessionNums);
%         for i=1:length(sessionNums)
%             ind=find(belowSe==sessionNums(i));
%             belowSe2(ind)=i;
%         end
%         for i=1:length(sessionNums)
%             ind=find(belowSe==sessionNums(i));
%             yvals=belowC(ind);
%             plot(i,yvals,'s',...
%                 'MarkerFaceColor',[currcol(floor(i*2.5),1) currcol(floor(i*2.5),2) currcol(floor(i*2.5),3)],...
%                 'MarkerEdgeColor',[currcol(floor(i*2.5),1) currcol(floor(i*2.5),2) currcol(floor(i*2.5),3)]);hold on
%         end
%         % plot(aroundSe2,aroundC,markerSh,'MarkerFaceColor',setFill,'LineStyle','none','MarkerEdgeColor',markerCol);hold on
%         brob = robustfit(belowSe2,belowC);
%         plot(xvals,brob(1)+brob(2)*xvals,'k','LineWidth',2);hold on
%         titletext=['channels with low initial c50 (<',num2str(sampleContrast-margin),'%), p=',num2str(allpC3m(2))];
%         title(titletext);
%         ylabel('c50');
%         xlabel('session number');
%         plot([1 length(sessionNums)],[sampleContrast sampleContrast],'--k')
%         yLimVals=get(gca,'YLim');
%         ind=find(yLimVals==sampleContrast);
%         if ind==1            
%             ylim([yLimVals(1)-5,yLimVals(2)]);
%         elseif ind==2
%             ylim([yLimVals(1),yLimVals(2)+5]);
%         end
%         figName=[folder,addText,'\roc_c50_ch_below_',num2str(sampleContrast-margin)];
%         printtext=['print -dtiff ',figName];
%         eval(printtext);
%         hgsave(figName);
%     end
% end
% 
%     %grand c50 figure:
%     figure
%     for i=1:length(sessionNums)
%         index=find(aboveSe==sessionNums(i));
%         meanC(i)=mean(aboveC(index));
%     end
%     xvals=1:length(sessionNums);
%     for i=1:length(sessionNums)
%         ind=find(aboveSe==sessionNums(i));
%         aboveSe2(ind)=i;
%     end
%     if allpC1(2)<0.05
%         setFill=markerCol;
%     else
%         setFill='none';
%     end
%     plot(aboveSe2,aboveC,markerSh,'MarkerFaceColor',setFill,'LineStyle','none','MarkerEdgeColor',markerCol);hold on
%     brob = robustfit(aboveSe2,aboveC);
%     plot(aboveSe2,brob(1)+brob(2)*aboveSe2,'k','LineWidth',2);hold on
%     titletext=['channels with high initial c50 (> ',num2str(sampleContrast+margin),'%), p=',num2str(allpC1(2))];
%     title(titletext);
%     for i=1:length(sessionNums)
%         index=find(aroundSe==sessionNums(i));
%         meanC(i)=mean(aroundSe(index));
%     end
%     xvals=1:length(sessionNums);
%     for i=1:length(sessionNums)
%         ind=find(aroundSe==sessionNums(i));
%         aroundSe2(ind)=i;
%     end
%     if allpC1m(2)<0.05
%         setFill=markerCol;
%     else
%         setFill='none';
%     end
%     plot(aroundSe2,aroundC,markerSh,'MarkerFaceColor',setFill,'LineStyle','none','MarkerEdgeColor',markerCol);hold on
%     brob = robustfit(aroundSe2,aroundC);
%     plot(aroundSe2,brob(1)+brob(2)*aroundSe2,'k','LineWidth',2);hold on
%     titletext=['channels with high  (> ',num2str(sampleContrast+margin),'%) and intermediate (>',num2str(sampleContrast-margin),'% & <',num2str(sampleContrast+margin),'%) initial c50, p=',num2str(allpC2m(2))];
%     title(titletext);
%     ylabel('c50');
%     xlabel('session number');
% 
% 
% % %rough check of correlation:
% b=[aboveSe' aboveC']
% [coefficientsC1 pC1]=corrcoef(b);%c50
% % b=[aboveSe' aboveSl']
% % [coefficientsS1 pS1]=corrcoef(b);%slope
% % b=[aroundSe' aroundC']
% % [coefficientsC2 pC2]=corrcoef(b);%c50
% % b=[aroundSe' aroundSl']
% % [coefficientsS2 pS2]=corrcoef(b);%slope
% % b=[belowSe' belowC']
% % [coefficientsC3 pC3]=corrcoef(b);%c50
% % b=[belowSe' belowSl']
% % [coefficientsS3 pS3]=corrcoef(b);%slope
% % coefficients(1)=coefficientsC1(2);%c50
% % coefficients(2)=coefficientsC2(2);
% % % coefficients(3)=coefficientsC3(2);
% % coefficients(4)=coefficientsS1(2);%slope
% % coefficients(5)=coefficientsS2(2);
% % % coefficients(6)=coefficientsS3(2);
% % p(1)=pC1(2);%c50
% % p(2)=pC2(2);
% % % p(3)=pC3(2);
% % p(4)=pS1(2);%slope
% % p(5)=pS2(2);
% % % p(6)=pS3(2);
% 
% %closer examination of correlation:
% for i=1:length(cellCategory)
%     if cellCategory(i)==1
%         figure(1)
%         ind=find(all_roc14(:,1)==channels(i));
%         brob = robustfit(all_roc14(ind,2),all_roc14(ind,18));
%         plot(all_roc14(ind,2),all_roc14(ind,18),'ok');hold on
%         plot(all_roc14(ind,2),brob(1)+brob(2)*all_roc14(ind,2),'k','LineWidth',2);hold on
%         aboveC1(i)=brob(1)+brob(2)*307;%first data point on fitted line
%         aboveC2(i)=brob(1)+brob(2)*342;%last data point on fitted line
%         figure(2)
%         brob = robustfit(all_roc14(ind,2),all_roc14(ind,17));
%         plot(all_roc14(ind,2),brob(1)+brob(2)*all_roc14(ind,2),'k','LineWidth',2);hold on
%         aboveSl1(i)=brob(1)+brob(2)*307;%first data point on fitted line
%         aboveSl2(i)=brob(1)+brob(2)*342;%last data point on fitted line
%     elseif cellCategory(cellInd)==0   
%         figure(1)     
%         ind=find(all_roc14(:,1)==channels(i));
%         brob = robustfit(all_roc14(ind,2),all_roc14(ind,18));
%         plot(all_roc14(ind,2),all_roc14(ind,18),'ok');hold on
%         plot(all_roc14(ind,2),brob(1)+brob(2)*all_roc14(ind,2),'k','LineWidth',2);hold on
%         aroundC1(i)=brob(1)+brob(2)*307;%first data point on fitted line
%         aroundC2(i)=brob(1)+brob(2)*342;%last data point on fitted line
%         figure(2)
%         brob = robustfit(all_roc14(ind,2),all_roc14(ind,17));
%         plot(all_roc14(ind,2),brob(1)+brob(2)*all_roc14(ind,2),'k','LineWidth',2);hold on
%         aroundSl1(i)=brob(1)+brob(2)*307;%first data point on fitted line
%         aroundSl2(i)=brob(1)+brob(2)*342;%last data point on fitted line
%     elseif cellCategory(cellInd)==-1
%         figure(1)
%         ind=find(all_roc14(:,1)==channels(i));
%         brob = robustfit(all_roc14(ind,2),all_roc14(ind,18));
%         plot(all_roc14(ind,2),all_roc14(ind,18),'ok');hold on
%         plot(all_roc14(ind,2),brob(1)+brob(2)*all_roc14(ind,2),'k','LineWidth',2);hold on
%         belowC1(i)=brob(1)+brob(2)*307;%first data point on fitted line
%         belowC2(i)=brob(1)+brob(2)*342;%last data point on fitted line
%         figure(2)
%         brob = robustfit(sessionNums,all_roc14(ind,17));
%         plot(all_roc14(ind,2),brob(1)+brob(2)*all_roc14(ind,2),'k','LineWidth',2);hold on
%         belowSl1(i)=brob(1)+brob(2)*307;%first data point on fitted line
%         belowSl2(i)=brob(1)+brob(2)*342;%last data point on fitted line
%     end
% end
% length(find(cellCategory==1))
% length(find(cellCategory==0))
% length(find(cellCategory==-1))
% 
% aboveC1=aboveC1(find(aboveC1~=0));
% aboveC2=aboveC2(find(aboveC2~=0));
% aboveSl1=aboveSl1(find(aboveSl1~=0));
% aboveSl2=aboveSl2(find(aboveSl2~=0));
% [h,p,ci,stats]=ttest(aboveC1,aboveC2);
% [h,p,ci,stats]=ttest(aboveSl1,aboveSl2);
% %3% margin:
% %coef p value: 0.8014
% %slope p value: 0.9714
% aroundC1=aroundC1(find(aroundC1~=0));
% aroundC2=aroundC2(find(aroundC2~=0));
% aroundSl1=aroundSl1(find(aroundSl1~=0));
% aroundSl2=aroundSl2(find(aroundSl2~=0));
% 
% %calculate a c50 for each session for each channel.
% %For each channel, calculate mean c50 across first 5 sessions.
% %Divide channels according to whether mean c50 falls above or below 30%.
% %Try categorising based on whether the value falls a certain small
% %percentage above or below 30%, e.g. 3% or 2%.
% %Check for correlation between c50s and session, using corrcoef.
% 
% loadMatText=['load ',folder,'\av_roc_chs14.mat av_roc14'];
% eval(loadMatText);
% X0=[30 2];
% options = optimset('Display','off','MaxFunEvals',10^6,'MaxIter',10^6,'TolFun',1.0E-6,'TolX',1.0E-6);
% figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
% for i=1:size(av_roc14,1)
%     perfvals=av_roc14(i,1:14);%neurometric performance for each session
%     perfvals2=VALUES(i,1:14);%psychometric performance for each session
%     [X,fval]=fminsearch('weib_sim_min_max',X0,options,testContrast,perfvals);
%     [X2,fval2]=fminsearch('weib_sim_min_max',X0,options,testContrast,perfvals2);
%     slopeNeuro(1,i)=X(2);
%     c50(1,i)=X(1);
%     slopePsycho(1,i)=X2(2);
%     c502(1,i)=X2(1);
%     subplot(ceil((size(av_roc14,1)+2)/5),5,i);
%     plot(testContrast,perfvals,'ok');
%     hold on
%     yvals=max(perfvals)-(max(perfvals)-min(perfvals))*exp(-((testContrast./X(1)).^X(2)));
%     plot(testContrast,yvals,'r');
%     line(sampleContrast,0:0.01:1);
%     ylim([min(perfvals),max(perfvals)]);
%     subplottitle=num2str(sessionNums(1,i));
%     title(subplottitle);
%     if i==1
%         ptext=sprintf('%s',psychoname);
%         orient landscape
%         yLimVals=get(gca,'YLim');
%         text('Position',[-10 yLimVals(2)+0.2],'FontSize',9,'String',ptext);
%     end
% end
% % subplot(ceil((size(av_roc14,1)+2)/5),5,i+1);
% % plot(sessionNums,slopeNeuro,'xk');
% % b=[sessionNums' slopeNeuro']
% % [coefficients1 p1]=corrcoef(b);
% % ptext=['coef: ',num2str(coefficients1(2)),' p: ',num2str(p1(2))];
% % title(ptext);
% % 
% % subplot(ceil((size(av_roc14,1)+2)/5),5,i+2);
% % plot(sessionNums,c50,'xk');
% % c=[sessionNums' c50']
% % [coefficients1 p1]=corrcoef(c);
% % ptext=['coef: ',num2str(coefficients1(2)),' p: ',num2str(p1(2))];
% % title(ptext);
% % printtext=['print -dpng ',folder,'\meanROCslopes14'];
% % eval(printtext);
% 
% plotCorrFigsOnly=1;
% if plotCorrFigsOnly==1
%     figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
%     subplot(1,3,1);
%     b=[sessionNums' slopeNeuro']
%     [coefficients1 p1]=corrcoef(b);
%     setFill='auto';
%     if p1(2)<0.05
%         plot(sessionNums,slopeNeuro,'Marker',markerSh,'MarkerFaceColor',markerCol,'LineStyle','none','MarkerEdgeColor',markerCol);
%         brob = robustfit(sessionNums,slopeNeuro);
%         plot(sessionNums,brob(1)+brob(2)*sessionNums,'k','LineWidth',2);
%     else
%         plot(sessionNums,slopeNeuro,markerSh,'MarkerFaceColor',setFill,'LineStyle','none','MarkerEdgeColor',markerCol);
%     end
%     ptext=['coef: ',num2str(coefficients1(2)),' p: ',num2str(p1(2))];
%     title(ptext);
%     
%     subplot(1,3,2);
%     c=[sessionNums' c50']
%     [coefficients1 p2]=corrcoef(c);
%     setFill='auto';
%     if p2(2)<0.05
%         plot(sessionNums,c50,'Marker',markerSh,'MarkerFaceColor',markerCol,'LineStyle','none','MarkerEdgeColor',markerCol);hold on
%         brob = robustfit(sessionNums,c50);
%         plot(sessionNums,brob(1)+brob(2)*sessionNums,'k','LineWidth',2);
%     else
%         plot(sessionNums,c50,markerSh,'MarkerFaceColor',setFill,'LineStyle','none','MarkerEdgeColor',markerCol);
%     end
%     ptext=['coef: ',num2str(coefficients1(2)),' p: ',num2str(p2(2))];
%     title(ptext);
%     
%     subplot(1,3,3);
%     d=[slopePsycho' slopeNeuro']
%     [coefficients1 p3]=corrcoef(d);
%     setFill='auto';
%     if p3(2)<0.05
%         plot(slopePsycho,slopeNeuro,'Marker',markerSh,'MarkerFaceColor',markerCol,'LineStyle','none','MarkerEdgeColor',markerCol);
%         brob = robustfit(slopePsycho,slopeNeuro);
%         plot(slopePsycho,brob(1)+brob(2)*sessionNums,'k','LineWidth',2);
%     else
%         plot(slopePsycho,slopeNeuro,'Marker',markerSh,'MarkerFaceColor',setFill,'LineStyle','none','MarkerEdgeColor',markerCol);
%     end
%     ptext=['coef: ',num2str(coefficients1(2)),' p: ',num2str(p3(2))];
%     title(ptext);
%     printtext=['print -dpng ',folder,'\meanROCslopes14'];
%     eval(printtext);
% end
% 
% % 
% % 
% % matText=['save ',folder,'\all_roc_chs.mat all_roc']
% % eval(matText);
% % av_roc=zeros(length(sessionNums),6);
% % for i=1:length(sessionNums)
% %     ind=find(all_roc(:,1)==sessionNums(i));
% %     av_roc(i,:)=mean(all_roc(ind,2:7));
% % end
% % matText=['save ',folder,'\av_roc_chs.mat av_roc']
% % eval(matText);
% % 
% % for cond=1:6
% %     a=[sessionNums' av_roc(:,cond)];
% %     [coefficients1 p1]=corrcoef(a);
% %     coefficients(cond)=coefficients1(2);
% %     ps(cond)=p1(2);
% % end
% 
% load F:\blanco\v4_1_roc_analysis\excludeSTNbest_sessions\all_roc_chs6.mat all_roc
% load F:\blanco\v1_roc_analysis\all_roc_chs6.mat all_roc
% for cond=1:6
%     pooled_roc=[];
%     pooled_sessionNums=[];
%     for i=1:length(sessionNums)
%         ind=find(all_roc(:,1)==sessionNums(i));
%         pooled_roc=[pooled_roc;all_roc(ind,cond+1)];
%         pooled_sessionNums=[pooled_sessionNums;all_roc(ind,1)];
%     end
%     a=[pooled_sessionNums pooled_roc];
%     [coefficients1 p1]=corrcoef(a);
%     coefficients(cond)=coefficients1(2);
%     ps(cond)=p1(2);
% end
% %V4_1:
% %coefs:
% %[-0.0439,0.0580,0.0266,0.1499,0.08612,0.06223;]
% %ps: 
% %[0.3001,0.1712,0.5309,0.00038057,0.04199,0.14204;]
% 
% %V1_1:
% %coefs:
% %[0.02529,0.1468,0.006772,0.10423,0.05085,0.11989;]
% %ps:
% %[0.6658,0.01176,0.90795,0.07436,0.38499,0.039939;]
% 
% fighandle1=  figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
% set(fighandle1, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
% plot(1:size(av_roc,1),av_roc(:,1),':xy');hold on
% plot(1:size(av_roc,1),av_roc(:,2),':xc');hold on
% plot(1:size(av_roc,1),av_roc(:,3),':xm');hold on
% plot(1:size(av_roc,1),av_roc(:,4),':xb');hold on
% plot(1:size(av_roc,1),av_roc(:,5),':xr');hold on
% plot(1:size(av_roc,1),av_roc(:,6),':xk');hold on
% lsline
% 
% yLimVals=get(gca,'ylim');
% xLimVals=get(gca,'xlim');
% unitSpace=(yLimVals(2)-yLimVals(1))/30;
% text('Position',[xLimVals(1)+(xLimVals(2)-xLimVals(1))/10 yLimVals(2)-unitSpace*6],'FontSize',9,'String',[num2str(testContrast(5)),'%'],'Color',[1 1 0]);
% text('Position',[xLimVals(1)+(xLimVals(2)-xLimVals(1))/10 yLimVals(2)-unitSpace*5],'FontSize',9,'String',[num2str(testContrast(6)),'%'],'Color',[0 1 1]);
% text('Position',[xLimVals(1)+(xLimVals(2)-xLimVals(1))/10 yLimVals(2)-unitSpace*4],'FontSize',9,'String',[num2str(testContrast(7)),'%'],'Color',[1 0 1]);
% text('Position',[xLimVals(1)+(xLimVals(2)-xLimVals(1))/10 yLimVals(2)-unitSpace*3],'FontSize',9,'String',[num2str(testContrast(8)),'%'],'Color',[0 0 1]);
% text('Position',[xLimVals(1)+(xLimVals(2)-xLimVals(1))/10 yLimVals(2)-unitSpace*2],'FontSize',9,'String',[num2str(testContrast(9)),'%'],'Color',[1 0 0]);
% text('Position',[xLimVals(1)+(xLimVals(2)-xLimVals(1))/10 yLimVals(2)-unitSpace],'FontSize',9,'String',[num2str(testContrast(10)),'%'],'Color',[0 0 0]);
% 
% set(gca,'XTick',1:size(av_roc,1))
% set(gca,'XTickLabel',sessionNums);
% 
% title('mean ROC values across channels');
% 
% orient landscape
% currdir=cd
% printdir=['cd ',folder];
% eval(printdir);
% rocFigName='mean_roc_all_chs';%where 'w' stands for 'time window'
% printtext=sprintf('print -dpng %s',rocFigName);
% eval(printtext);
% chdirtext=sprintf('cd ''%s''',currdir);
% eval(chdirtext);
% cd
% 
% %V4:
% % 33%
% % coefficients1 =0.3185
% % p1 =0.1386
% % 32%
% % coefficients1 =0.2779
% % p1 =0.1992
% % 31%
% % coefficients1 =0.6139
% % p1 =0.0018
% % 29%
% % coefficients1 =0.0473
% % p1 =0.8304
% % 28%
% % coefficients1 =0.1703
% % p1 =0.4373
% % 27%
% % coefficients1 =-0.2143
% % p1 =0.3261
% 
% 
% %V1:
% % coefficients from cond 5 to 10:
% % [0.1875,0.4605,-0.0086767,0.35086,0.14967,0.3957;]
% % p values from cond 5 to 10:
% % [0.5033,0.08406,0.9755,0.19976,0.5944,0.14425;]