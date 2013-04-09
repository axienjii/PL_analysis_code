function plot_all_bj_behav_figs(roving,areas)
%Written by Xing 27/10/12.
%Plot figures for both monkeys, for behavioural paper.

printFigs=1;
onExternalHD=0;
if onExternalHD==1
    rootFolder='G:\PL_backup_060413';
else
    rootFolder='F:';
end
animalTexts=[{'subject B'} {'subject J'}];
animals=[{'blanco'} {'jack'}];
if roving==0
    areaTexts=[{'V4'} {'V1'}];
    areas=[{'v4_1'} {'v4_2'} {'v1_1'}];
elseif roving==1
    areaTexts={'20' '30' '40'};
end
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        [sampleContrasts allTestContrasts]=area_metadata(area);
        for sampleInd=1:length(sampleContrasts)
            sampleContrast=sampleContrasts(sampleInd);
            testContrast=allTestContrasts(sampleInd,:);
            
            analysisFolderAppend=[];
            midInd=find(testContrast<sampleContrast);
            midInd=midInd(end);
            loadText=['load ',rootFolder,'\PL\psycho_data\',animal,'\allMeanPerf\allMeanPerf_',area,'_',num2str(sampleContrast),'.mat allMeanPerf'];
            eval(loadText)
            if strcmp(area,'v4_1')
                allMeanPerf=allMeanPerf(1:end-1,:);%exclude session with horizontal stimulus
                maxV4_1(animalInd)=size(allMeanPerf,1);
            end
            if roving==0
                if areaInd<3
                    plotAreaInd=1;%V4_1 & V4_2
                else
                    plotAreaInd=2;%V1
                end
            elseif roving==1
                plotAreaInd=sampleInd;
            end
            
            %values stored in allMeanPerf:
            %1: session
            %2: average perf across conditions
            %3 to 2+numconds: prop report higher
            %3+numconds: PSE
            %4+numconds: slope
            %5+numconds to 6+numconds: max and min
            %7+numconds to 8+numconds: average correct trial RT across conditions and std of all RT
            %9+numconds to 8+numconds*2: mean correct trial RT per condition
            %9+numconds*2 to 8+numconds*3: std of correct trial RT per condition
            %9+numconds*3 to 15+numconds*3: mean(all_ave_RT) mean(fix1_dur) mean(fix2_dur) mean(fix3_dur) mean(sample_dur) mean(test_dur) mean(trial_dur)
            %16+numconds*4 to 22+numconds*3: std of variables in previous line
            %23+numconds*3 to 24+numconds*3: average error trial RT across conditions and std of all error trial RT
            %25+numconds*3 to 24+numconds*4: mean error trial RT per condition
            %25+numconds*4 to 24+numconds*5: std of error trial RT per condition
            %25+numconds*5 to 26+numconds*5: mean and std of RTs combined across all conditions, and across correct and error trials
            X0=[2 30 0.2 0.1];
            options = optimset('Display','off','MaxFunEvals',Inf,'MaxIter',10^6,'TolFun',1.0E-6,'TolX',1.0E-5);
            xvals=0:1:testContrast(end)+10;
            numconds=length(testContrast);
            numsessions=size(allMeanPerf,1);
            for i=1:numsessions
                VALUES=allMeanPerf(:,3:2+numconds);
            end
            
            %psychometric curves for all sessions
            psychoCurves(animalInd+2*(areaInd-1))=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
            for i=1:numsessions
                figure(psychoCurves(animalInd+2*(areaInd-1)))
                prop_corr=VALUES(i,1:numconds);%psychometric performance for each session
                X=fminsearch(@fit_weibull,X0,[],testContrast,prop_corr,[],'least_square',[1 0 1 0],[20 0 max(prop_corr) 0],[1 0 0 0],[0 0 0 0]);
                subplot(6,6,i);%check
                plot(testContrast,prop_corr,'ok');
                hold on
                yvals=1-X(4)-X(3).*exp(-(xvals./X(2)).^X(1));
                plot(xvals,yvals,'k');
                line(sampleContrast,0:0.01:1);
                line(allMeanPerf(i,3+length(testContrast)),0:0.01:1,'Color','r');
                ylim([min(prop_corr),max(prop_corr)]);
                %     subplottitle=num2str(allMeanPerf(i,1));
                subplottitle=num2str(i);
                title(subplottitle);
                set(gca,'FontSize',[6],'YLim',[0,1.01],'XLim',[0,testContrast(end)+10],'YTickMode','manual');%'YTick',[0.1]
                if i==1
                    %         ptext=sprintf('%s',psychoname);
                    orient landscape
                    yLimVals=get(gca,'YLim');
                    %         text('Position',[-10 yLimVals(2)+0.4],'FontSize',9,'String',ptext);
                end
                ind=find(testContrast<sampleContrast);
                ind=ind(end);
                perfvals=[1-VALUES(i,1:ind) VALUES(i,1+ind:end)];
                stdPerf(i)=std(perfvals);
                if strcmp(animal,'jack')&&strcmp(area,'v4_1')
                    if i==1||i==2
                        %figure with 2 example sessions- 1 early and 1 late:
                        if i==1
                            psychoCurvesj2sess=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.5]);
                            subplot(1,2,1);
                            prop_corr=[0.64516129 0.609756098 0.632911392 0.552486188 0.462962963 0.595238095 0.534759358 0.515463918 0.505050505 0.523560209 0.421940928 0.469483568 0.746268657 0.925925926];
                            prop_corr=[1-prop_corr(1:7) prop_corr(8:end)];
                            X=fminsearch(@fit_weibull,X0,[],testContrast,prop_corr,[],'least_square',[0 0 0 0],[],[0 0 0 0],[]);
                            yvals=1-X(4)-X(3).*exp(-(xvals./X(2)).^X(1));
                            PSE=X(2).*(-log((0.5-X(4))/X(3))).^(1/X(1));
                        else
                            figure(psychoCurvesj2sess)
                            subplot(1,2,2);
                            prop_corr=[0.917431193000000,0.961538462000000,0.934579439000000,0.775193798000000,0.757575758000000,0.662251656000000,0.636942675000000,0.520833333000000,0.641025641000000,0.680272109000000,0.793650794000000,0.925925926000000,0.990099010000000,1];
                            prop_corr=[1-prop_corr(1:7) prop_corr(8:end)];
                            X=fminsearch(@fit_weibull,X0,[],testContrast,prop_corr,[],'least_square',[0 0 0 0],[],[0 0 0 0],[]);
                            yvals=1-X(4)-X(3).*exp(-(xvals./X(2)).^X(1));
                            PSE=X(2).*(-log((0.5-X(4))/X(3))).^(1/X(1));
                        end
                        plot(testContrast,prop_corr,'ok');
                        hold on
                        plot(xvals,yvals,'k');
                        line(sampleContrast,0:0.01:1);
                        plot([PSE PSE],[0 1],'Color','r','LineStyle','--');
                        %         plot([allMeanPerf(i,3+length(testContrast)) allMeanPerf(i,3+length(testContrast))],[0 1],'Color','r','LineStyle','--');
                        ylim([min(prop_corr),max(prop_corr)]);
                        %     subplottitle=num2str(allMeanPerf(i,1));
                        subplottitle=num2str(i);
                        title(subplottitle);
                        set(gca,'FontSize',[6],'YLim',[0,1.01],'XLim',[0,testContrast(end)+10],'YTickMode','manual');%'YTick',[0.1]
                    end
                end
            end
            % printtext=sprintf('print -dpng %s',psychoname);
            % eval(printtext);
            
            colmapText=colormap(jet(size(testContrast,2)));
            colmapText=[colmapText(1,:);132/255 22/255 216/255;202/255 65/255 223/255;colmapText(3:7,:);157/255 212/255 61/255;colmapText(10:12,:);178/255 111/255 12/255;colmapText(end,:)];
            
            % %simple proportion correct figure without error bars
            if animalInd==1&&areaInd==1&&sampleInd==1
                if roving==0
                    pcpseslFig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
                elseif roving==1
                    pcpseslFig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.5, 0.8]);
                end
            else
                figure(pcpseslFig);
            end
            subplot(length(areaTexts),3,1+(plotAreaInd-1)*3);
            animalMarkerFaceColors={'none' 'k'};
            animalMarkerEdgeColors={'auto' 'k'};
            if strcmp(area,'v4_2')
                plot(3+maxV4_1(animalInd):size(allMeanPerf,1)+maxV4_1(animalInd)+2,allMeanPerf(:,2)*100,'ok','LineStyle','none','MarkerFaceColor',animalMarkerFaceColors{animalInd},'MarkerEdgeColor',animalMarkerEdgeColors{animalInd});
            else
                plot(1:size(allMeanPerf,1),allMeanPerf(:,2)*100,'ok','LineStyle','none','MarkerFaceColor',animalMarkerFaceColors{animalInd},'MarkerEdgeColor',animalMarkerEdgeColors{animalInd});
            end
            if strcmp(animal,'jack')&&strcmp(area,'v1_1')
                ylim([65 95]);
                xlim([0 24]);
            end
            hold on
            if strcmp(area,'v4_1')%horizontal Gabor session
                animalMarkerFaceColors={'none' 'b'};
                plot(size(allMeanPerf,1),allMeanPerf(end,2)*100,'Marker','o','Color','b','LineStyle','none','MarkerFaceColor',animalMarkerFaceColors{animalInd});
            end
            if ~strcmp(area,'v4_2')
                [chiselinearperf(animalInd+2*(plotAreaInd-1)) chiseperf(animalInd+2*(plotAreaInd-1))]=bj_linearexpo_fitting(testContrast,allMeanPerf(:,2)*100,0,0,'ROC');
            end
            
            % %simple proportion correct without error bars,
            % for first and last 30% of trials within session
            if animalInd==1&&areaInd==1&&sampleInd==1
                if roving==0
                    pcpsesl30Fig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
                elseif roving==1
                    pcpsesl30Fig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.5, 0.8]);
                end
            else
                figure(pcpsesl30Fig);
            end
            subplot(length(areaTexts),3,1+(plotAreaInd-1)*3);
            animalMarkerFaceColors={'none' 'k'};
            if strcmp(area,'v4_2')
                plot(3+maxV4_1(animalInd):size(allMeanPerf,1)+maxV4_1(animalInd)+2,allMeanPerf(:,27+numconds*5),'Marker','o','Color','k','LineStyle','none','MarkerFaceColor',animalMarkerFaceColors{animalInd});
            else
                plot(1:size(allMeanPerf,1),allMeanPerf(:,27+numconds*5),'Marker','o','Color','k','LineStyle','none','MarkerFaceColor',animalMarkerFaceColors{animalInd});
            end
            hold on
            animalMarkerFaceColors={'none' 'r'};
            if strcmp(area,'v4_2')
                plot(3+maxV4_1(animalInd):size(allMeanPerf,1)+maxV4_1(animalInd)+2,allMeanPerf(:,34+numconds*6),'Marker','o','Color','r','LineStyle','none','MarkerFaceColor',animalMarkerFaceColors{animalInd});
            else
                plot(1:size(allMeanPerf,1),allMeanPerf(:,34+numconds*6),'Marker','o','Color','r','LineStyle','none','MarkerFaceColor',animalMarkerFaceColors{animalInd});
                [hperfIndpair(plotAreaInd+length(areaTexts)*(animalInd-1),1),pperfIndpair(plotAreaInd+length(areaTexts)*(animalInd-1),1)]=ttest(allMeanPerf(:,27+numconds*5),allMeanPerf(:,34+numconds*6));
            end
            
            % index of difference in simple proportion correct without error bars,
            % between first and last 30% of trials within session
            if animalInd==1&&areaInd==1&&sampleInd==1
                if roving==0
                    pcpsesl30indFig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
                elseif roving==1
                    pcpsesl30indFig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.5, 0.8]);
                end
            else
                figure(pcpsesl30indFig);
            end
            subplot(length(areaTexts),3,1+(plotAreaInd-1)*3);
            pcInd=(allMeanPerf(:,34+numconds*6)-allMeanPerf(:,27+numconds*5))./allMeanPerf(:,27+numconds*5);
            animalMarkerFaceColors={'none' 'k'};
            if strcmp(area,'v4_2')
                plot(3+maxV4_1(animalInd):size(allMeanPerf,1)+maxV4_1(animalInd)+2,pcInd,'Marker','o','Color','k','LineStyle','none','MarkerFaceColor',animalMarkerFaceColors{animalInd});
            else
                plot(1:size(allMeanPerf,1),pcInd,'Marker','o','Color','k','LineStyle','none','MarkerFaceColor',animalMarkerFaceColors{animalInd});
                [hperfInd(plotAreaInd+length(areaTexts)*(animalInd-1),1),pperfInd(plotAreaInd+length(areaTexts)*(animalInd-1),1)]=ttest(pcInd);
            end
            hold on
            
            %proportion correct with error bars
%             if animalInd==1&&areaInd==1&&sampleInd==1
%                 pcpseslFig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
%             else
%                 figure(pcpseslFig);
%             end
%             errorbar(1:size(allMeanPerf,1),allMeanPerf(:,2),stdPerf,'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','Color','k');
%             hold on
            
            %
%             figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
%             for i=1:numconds
%                 plot(allMeanPerf(:,1),allMeanPerf(:,2+i),'Color',[1/i 1/i 1/i],'LineStyle','none','Marker','o','MarkerFaceColor',[1/i 1/i 1/i],'MarkerEdgeColor',[1/i 1/i 1/i]);hold on
%             end
            
            %RT, divided by condition:
            if animalInd==1&&areaInd==1&&sampleInd==1
                rt_condFig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
            else
                figure(rt_condFig);
            end
            subplot(length(areaTexts),2,animalInd+2*(plotAreaInd-1));
            f = fittype('a+b*x');
            colTexts='ymcrgbk';
            markerTexts='+x';
            markerSizes=[5 6];
            for i=1:numconds
                if i<=ind
                    colText=colTexts(i);
                    markerText=markerTexts(1);%sample of lower contrast than test
                    markerS=markerSizes(1);
                else
                    colText=colTexts(i-ind);
                    markerText=markerTexts(2);%sample of higher contrast than test
                    markerS=markerSizes(2);
                end
                fit3 = fit(allMeanPerf(:,1),allMeanPerf(:,8+numconds+i),f,'StartPoint',[1 1],'Robust','on');
                if strcmp(area,'v4_2')
                    errorbar(3+maxV4_1(animalInd):size(allMeanPerf,1)+maxV4_1(animalInd)+2,allMeanPerf(:,8+numconds+i),allMeanPerf(:,8+numconds*2+i),'Color',colText,'LineStyle','none','Marker',markerText,'MarkerFaceColor',colText,'MarkerEdgeColor',colText,'MarkerSize',markerS);hold on%'MarkerFaceColor',[1/i 1/i 1/i],
                else
                    errorbar(1:size(allMeanPerf,1),allMeanPerf(:,8+numconds+i),allMeanPerf(:,8+numconds*2+i),'Color',colText,'LineStyle','none','Marker',markerText,'MarkerFaceColor',colText,'MarkerEdgeColor',colText,'MarkerSize',markerS);hold on%'MarkerFaceColor',[1/i 1/i 1/i],
                end
                %     plot(fit3,'b:');
                colText=['-',colText];
                %     plot(fit3,colText);hold on%,'LineWidth',1
                % lsline
                yLimVals=get(gca,'ylim');
                xLimVals=get(gca,'xlim');
                unitSpace=(yLimVals(2)-yLimVals(1))/30;
            end
            
            %RT error, divided by condition:
            if animalInd==1&&areaInd==1&&sampleInd==1
                rterror_condFig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
            else
                figure(rterror_condFig);
            end
            subplot(length(areaTexts),2,animalInd+2*(plotAreaInd-1));
            f = fittype('a+b*x');
            colTexts='ymcrgbk';
            markerTexts='+x';
            for i=1:numconds
                if i<=ind
                    colText=colTexts(i);
                    markerText=markerTexts(1);%sample of lower contrast than test
                    markerS=markerSizes(1);
                else
                    colText=colTexts(i-ind);
                    markerText=markerTexts(2);%sample of higher contrast than test
                    markerS=markerSizes(2);
                end
                RTerr_noNan=allMeanPerf(:,24+numconds*3+i);
                RTerr_noNan(isnan(RTerr_noNan))=0;
                %     fit3 = fit(allMeanPerf(:,1),RTerr_noNan,f,'StartPoint',[1 1],'Robust','on');
                if strcmp(area,'v4_2')
                    errorbar(3+maxV4_1(animalInd):size(allMeanPerf,1)+maxV4_1(animalInd)+2,allMeanPerf(:,24+numconds*3+i),allMeanPerf(:,24+numconds*4+i),'Color',colText,'LineStyle','none','Marker',markerText,'MarkerFaceColor',colText,'MarkerEdgeColor',colText,'MarkerSize',markerS);hold on%'MarkerFaceColor',[1/i 1/i 1/i],
                else
                    errorbar(1:size(allMeanPerf,1),allMeanPerf(:,24+numconds*3+i),allMeanPerf(:,24+numconds*4+i),'Color',colText,'LineStyle','none','Marker',markerText,'MarkerFaceColor',colText,'MarkerEdgeColor',colText);hold on%'MarkerFaceColor',[1/i 1/i 1/i],
                end
                %     plot(fit3,'b:');
                colText=['-',colText];
                %     plot(fit3,colText);hold on%,'LineWidth',1
                % lsline
                yLimVals=get(gca,'ylim');
                xLimVals=get(gca,'xlim');
                unitSpace=(yLimVals(2)-yLimVals(1))/30;
            end
            
            % text('Position',[xLimVals(1)+(xLimVals(2)-xLimVals(1))/10 yLimVals(2)-unitSpace*6],'FontSize',9,'String',[num2str(testContrast(midInd-2)),'%'],'Color',[1 1 0]);
            % text('Position',[xLimVals(1)+(xLimVals(2)-xLimVals(1))/10 yLimVals(2)-unitSpace*5],'FontSize',9,'String',[num2str(testContrast(midInd-1)),'%'],'Color',[0 1 1]);
            % text('Position',[xLimVals(1)+(xLimVals(2)-xLimVals(1))/10 yLimVals(2)-unitSpace*4],'FontSize',9,'String',[num2str(testContrast(midInd)),'%'],'Color',[1 0 1]);
            % text('Position',[xLimVals(1)+(xLimVals(2)-xLimVals(1))/10 yLimVals(2)-unitSpace*3],'FontSize',9,'String',[num2str(testContrast(midInd+1)),'%'],'Color',[0 0 1]);
            % text('Position',[xLimVals(1)+(xLimVals(2)-xLimVals(1))/10 yLimVals(2)-unitSpace*2],'FontSize',9,'String',[num2str(testContrast(midInd+2)),'%'],'Color',[1 0 0]);
            % text('Position',[xLimVals(1)+(xLimVals(2)-xLimVals(1))/10 yLimVals(2)-unitSpace],'FontSize',9,'String',[num2str(testContrast(midInd+3)),'%'],'Color',[0 0 0]);
            
            
            %mean RT without error bars
            % figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8,
            % 0.8]);
            % plot(allMeanPerf(:,1),allMeanPerf(:,7+numconds)/1000,'ok','LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k');
            
            %mean RT for error trials with error bars
            % rterrorFig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
            if animalInd==1&&areaInd==1&&sampleInd==1
                if roving==0
                    rtFig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
                elseif roving==1
                    rtFig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.5, 0.8]);
                end
            else
                figure(rtFig);
            end
            subplot(length(areaTexts),2,animalInd+2*(plotAreaInd-1));
            if strcmp(area,'v4_2')
                errorbar(3+maxV4_1(animalInd):size(allMeanPerf,1)+maxV4_1(animalInd)+2,allMeanPerf(:,23+numconds*3)/1000,allMeanPerf(:,24+numconds*3)/1000,'LineStyle','none','Marker','o','MarkerEdgeColor','r','Color','r','MarkerFaceColor','r');
            else
                errorbar(1:size(allMeanPerf,1),allMeanPerf(:,23+numconds*3)/1000,allMeanPerf(:,24+numconds*3)/1000,'LineStyle','none','Marker','o','MarkerEdgeColor','r','Color','r','MarkerFaceColor','r');
                sessions=1:size(allMeanPerf,1);sessions=sessions';
                [rmeanRT(plotAreaInd*2,animalInd),pmeanRT(plotAreaInd*2,animalInd)]=corr(sessions,allMeanPerf(:,23+numconds*3)/1000);
            end
            hold on
            
            %mean RT with error bars
            figure(rtFig)
            subplot(length(areaTexts),2,animalInd+2*(plotAreaInd-1));
            if strcmp(area,'v4_2')
                errorbar(3+maxV4_1(animalInd):size(allMeanPerf,1)+maxV4_1(animalInd)+2,allMeanPerf(:,7+numconds)/1000,allMeanPerf(:,8+numconds)/1000,'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','Color','k');
            else
                errorbar(1:size(allMeanPerf,1),allMeanPerf(:,7+numconds)/1000,allMeanPerf(:,8+numconds)/1000,'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','Color','k');hold on
                sessions=1:size(allMeanPerf,1);sessions=sessions';
                [rmeanRT(plotAreaInd*2-1,animalInd),pmeanRT(plotAreaInd*2-1,animalInd)]=corr(sessions,allMeanPerf(:,7+numconds)/1000);
            end
            hold on
            
            %difference between mean RTs for correct and error trials
            diffMeanRT=(allMeanPerf(:,23+numconds*3)-allMeanPerf(:,7+numconds))./allMeanPerf(:,7+numconds);
            diffMeanRTabs(areaInd,animalInd)=mean(allMeanPerf(:,23+numconds*3)-allMeanPerf(:,7+numconds))/1000;%calculate difference between correct and error trial RTs for each session, find mean difference across sessions.
            diffMeanRTabs(areaInd,animalInd+2)=std(allMeanPerf(:,23+numconds*3)-allMeanPerf(:,7+numconds))/1000;%calculate difference between correct and error trial RTs for each session, find mean difference across sessions.
            if animalInd==1&&areaInd==1&&sampleInd==1
                diffmeanRTFig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.4]);%[left, bottom, width, height]
            else
                figure(diffmeanRTFig);
            end
            subplot(1,length(areaTexts),plotAreaInd);
            animalMarkerFaceColors={'none' 'k'};
            if strcmp(area,'v4_2')
                plot(3+maxV4_1(animalInd):size(allMeanPerf,1)+maxV4_1(animalInd)+2,diffMeanRT,'LineStyle','none','Marker','o','MarkerFaceColor',animalMarkerFaceColors{animalInd},'MarkerEdgeColor','k');
            else
                plot(1:size(allMeanPerf,1),diffMeanRT,'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor',animalMarkerFaceColors{animalInd});
                [hRTInd(animalInd,areaInd),pRTInd(animalInd,areaInd)]=ttest(diffMeanRT);
            end
            hold on
            
            %calculate normalised RT for correct and error trials, for each condition
            %(subtract mean RT for session- across conditions and correct/error trials-
            %from the mean RT per condition, separately for correct and error trials):
            normcRT=[];
            normeRT=[];
            for i=1:size(allMeanPerf,1)
                normcRT(i,:)=allMeanPerf(i,9+numconds:8+numconds*2)-allMeanPerf(i,25+numconds*5)/1000;%normalised RT for correct trials, per condition
                normeRT(i,:)=allMeanPerf(i,25+numconds*3:24+numconds*4)-allMeanPerf(i,25+numconds*5)/1000;%normalised RT for correct trials, per condition
            end
            diffCERT=(normeRT-normcRT)./abs(normcRT);%calculate difference between error trial RT and correct trial RT for each condition. 'CERT': correct; error; RT
            %note that diffCERT is array with num of session in rows, and columns 1 to
            %number of test contrast conditions. Many cells are 'NaN' for easy
            %conditions, because there are no or few error trials.
            lowerdiffCERT=reshape(diffCERT(:,1:midInd),1,[]);
            if roving==0 
                [hdiffCERT(areaInd*2-1,animalInd),pdiffCERT(areaInd*2-1,animalInd)]=ttest(lowerdiffCERT);
            elseif roving==1
                [hdiffCERT(sampleInd*2-1,animalInd),pdiffCERT(sampleInd*2-1,animalInd)]=ttest(lowerdiffCERT);
            end
            higherdiffCERT=reshape(diffCERT(:,midInd+1:end),1,[]);
            if roving==0
                [hdiffCERT(areaInd*2,animalInd),pdiffCERT(areaInd*2,animalInd)]=ttest(higherdiffCERT);
            elseif roving==1
                [hdiffCERT(sampleInd*2,animalInd),pdiffCERT(sampleInd*2,animalInd)]=ttest(higherdiffCERT);
            end
            if animalInd==1&&areaInd==1&&sampleInd==1
                rtcediffFig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
            else
                figure(rtcediffFig);
            end
            subplot(length(areaTexts),2,animalInd+2*(plotAreaInd-1));
            for i=1:numconds
                if i<=ind
                    colText=colTexts(i);
                    markerText=markerTexts(1);%sample of lower contrast than test
                    markerS=markerSizes(1);
                else
                    colText=colTexts(i-ind);
                    markerText=markerTexts(2);%sample of higher contrast than test
                    markerS=markerSizes(2);
                end
                plot(1:size(allMeanPerf,1),diffCERT(:,i),'Color',colText,'LineStyle','none','Marker',markerText,'MarkerFaceColor',colText,'MarkerEdgeColor',colText,'MarkerSize',markerS);hold on%'MarkerFaceColor',[1/i 1/i 1/i],
            end
            if plotAreaInd==2
                [xAxis]=get(gca,'xLim');
                line([xAxis(1) xAxis(2)],[0 0],'Color','k','LineStyle',':');
                %     ylim([-50 50]);
            end
            
            %ANOVA statistical analysis with condition number and correct/incorrect as
            %factors:
            if ~strcmp(area,'v4_2')
                RTs{animalInd+2*(plotAreaInd-1)}=[normcRT;normeRT];
                ceInd=[zeros(size(allMeanPerf,1),length(testContrast))+1;zeros(size(allMeanPerf,1),length(testContrast))+2];
                %             condInd=zeros(size(allMeanPerf,1)*2,length(testContrast));%individual condition numbers
                %             for i=1:length(testContrast)
                %                 condInd(:,i)=i;
                %             end
                lowerContrastConds=sum(testContrast<sampleContrast);
                condInd=[zeros(size(allMeanPerf,1),lowerContrastConds)+1 zeros(size(allMeanPerf,1),length(testContrast)-lowerContrastConds)+2];%conditions divided by lower and higher test contrast than sample contrast
                condInd=[condInd;condInd];%first half of rows: correct, second half: error. earlier columns: conditions with lower test contrast, later columns: higher test cnotrast
                sessInd=zeros(size(allMeanPerf,1)*2,length(testContrast));
                for j=1:2
                    for i=1:size(allMeanPerf,1)
                        sessInd(i+(size(allMeanPerf,1))*(j-1),:)=i;
                    end
                end
                RTs{animalInd+2*(plotAreaInd-1)}=reshape(RTs{animalInd+2*(plotAreaInd-1)},1,size(allMeanPerf,1)*2*length(testContrast));
                ceInd=reshape(ceInd,1,size(allMeanPerf,1)*2*length(testContrast));
                condInd=reshape(condInd,1,size(allMeanPerf,1)*2*length(testContrast));
                sessInd=reshape(sessInd,1,size(allMeanPerf,1)*2*length(testContrast));
                figure
%                 [p,table,stats] = anovan(RTs{animalInd+2*(plotAreaInd-1)},{ceInd,condInd,sessInd})
%                 [p,table,stats] = anovan(RTs{animalInd+2*(plotAreaInd-1)},{ceInd,condInd,sessInd},'model','interaction')
%                 [p,table,stats] = anovan(RTs{animalInd+2*(plotAreaInd-1)},{ceInd,condInd,sessInd},'model','full')
%                 [c,m,h] = multcompare(stats)
%                 [p,table,stats] = anovan(RTs{animalInd+2*(plotAreaInd-1)},{ceInd,condInd})
                [pRTceANOVA(:,animalInd+2*(plotAreaInd-1)),table,stats] = anovan(RTs{animalInd+2*(plotAreaInd-1)},{ceInd,condInd},'model','interaction')
%                 [p,table,stats] = anovan(RTs{animalInd+2*(plotAreaInd-1)},{ceInd,condInd},'model','full')
                [cRTceANOVA(animalInd+2*(plotAreaInd-1),:),m,h] = multcompare(stats)
            end
            
            %compare RT between first and last 30% of trials within each session, for
            %correct trials only, RT averaged across mean RTs for test contrast
            %conditions:
            if ~strcmp(area,'v4_2')
                if animalInd==1&&areaInd==1&&sampleInd==1
                    if roving==0
                        RT30=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
                    elseif roving==1
                        RT30=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.5]);
                    end
                else
                    figure(RT30);
                end
                subplot(1,length(areaTexts),plotAreaInd);
                animalMarkerFaceColors={'none' 'k'};
                plot(1:size(allMeanPerf,1),allMeanPerf(:,32+6*length(testContrast))/1000,'ko','MarkerEdgeColor','k','Color','k','MarkerFaceColor',animalMarkerFaceColors{animalInd});hold on
                plot(1:size(allMeanPerf,1),allMeanPerf(:,39+7*length(testContrast))/1000,'ro','MarkerEdgeColor','r','Color','r','MarkerFaceColor',animalMarkerFaceColors{animalInd});hold on
                
                if animalInd==1&&areaInd==1&&sampleInd==1
                    if roving==0
                    RT30ind=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
                    elseif roving==1
                    RT30ind=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.5]);
                    end
                else
                    figure(RT30ind);
                end
                subplot(1,length(areaTexts),plotAreaInd);
                plot(1:size(allMeanPerf,1),(allMeanPerf(:,32+6*length(testContrast))-allMeanPerf(:,39+7*length(testContrast)))./allMeanPerf(:,39+7*length(testContrast)),'ko','MarkerEdgeColor','k','MarkerFaceColor',animalMarkerFaceColors{animalInd});hold on
                if roving==0
                    line([0 30],[0 0],'LineStyle',':');
                elseif roving==1
                    line([0 35],[0 0],'LineStyle',':');
                end
               [h,pRTind(animalInd+2*(plotAreaInd-1)),ciRTind(animalInd+2*(plotAreaInd-1),:),stats] = ttest(allMeanPerf(:,32+6*length(testContrast)),allMeanPerf(:,39+7*length(testContrast)))
            end
            
            %proportion of trials monkey reported that test was higher contrast,
            %divided by condition:
            if animalInd==1&&areaInd==1&&sampleInd==1
                if roving==0
                    pc_condFig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
                elseif roving==1
                    pc_condFig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.5, 0.8]);
                end
            else
                figure(pc_condFig);
            end
            subplot(length(areaTexts),2,animalInd+2*(plotAreaInd-1));
            markerText=markerTexts(2);markerS=8;
            if strcmp(animal,'blanco')&&strcmp(area,'v1_1')
                xlim([0 17]);
            end
            %         if strcmp(animal,'jack')&&strcmp(area,'v4_2')%might have to shift into following i loop
            %             xlim([0 38]);
            %         end
            if strcmp(animal,'jack')&&strcmp(area,'v1_1')
                xlim([0 22]);
            end
            for i=1:numconds
                if strcmp(area,'v4_2')
                    plot(3+maxV4_1(animalInd):2+size(allMeanPerf,1)+maxV4_1(animalInd),allMeanPerf(:,2+i)*100,'Color',colmapText(i,:),'LineStyle','none','Marker',markerText,'MarkerFaceColor',colmapText(i,:),'MarkerEdgeColor',colmapText(i,:),'MarkerSize',markerS);hold on%'MarkerFaceColor',[1/i 1/i 1/i],
                else
                    if strcmp(animal,'blanco')&&strcmp(area,'v4_1')
                        xlim([0 size(allMeanPerf,1)+maxV4_1(animalInd)+4]);
                    elseif strcmp(animal,'jack')&&strcmp(area,'v4_2')
                        yLimVals=get(gca,'ylim');
                        xLimVals=get(gca,'xlim');
                        unitSpace=(yLimVals(2)-yLimVals(1))/30;
                        text('Position',[xLimVals(2)+(xLimVals(2)-xLimVals(1))/25 yLimVals(1)+unitSpace*i*2],'FontSize',9,'String',[markerText,'  ',num2str(testContrast(i)),'%'],'Color',colmapText(i,:));
                    end
                    plot(1:size(allMeanPerf,1),allMeanPerf(:,2+i)*100,'Color',colmapText(i,:),'LineStyle','none','Marker',markerText,'MarkerFaceColor',colmapText(i,:),'MarkerEdgeColor',colmapText(i,:),'MarkerSize',markerS);hold on%'MarkerFaceColor',[1/i 1/i 1/i],
%                     if strcmp(area,'v1_2_1')%try a linear fit as little learning occurs during roving stage 4
%                         plotLinear=1;
%                         startpoint5=0;
%                     else
                        plotLinear=0;
                        startpoint5=1;
%                     end
                    [chiselinearTemp(i,:) chiseTemp(i,:) coefperflinearTemp(i,:) coefperfTemp(i,:) aRSlinearTemp(i,:) aRSTemp(i,:)]=bj_linearexpo_fitting(testContrast,allMeanPerf(:,2+i)*100,i,startpoint5,'ROC',plotLinear);
                end
                %     yLimVals=get(gca,'ylim');
                %     xLimVals=get(gca,'xlim');
                %     unitSpace=(yLimVals(2)-yLimVals(1))/30;
            end
            if ~strcmp(area,'v4_2')
                chiselinear{animalInd+2*(plotAreaInd-1)}=chiselinearTemp;
                chise{animalInd+2*(plotAreaInd-1)}=chiseTemp;
                coefperflinear{animalInd+2*(plotAreaInd-1)}=coefperflinearTemp;
                coefperf{animalInd+2*(plotAreaInd-1)}=coefperfTemp;
                aRSlinear{animalInd+2*(plotAreaInd-1)}=aRSlinearTemp;
                aRS{animalInd+2*(plotAreaInd-1)}=aRSTemp;
            end
            ylim([0 100]);
            legend('hide');
            xlabel('');
            ylabel('');
            if strcmp(animal,'blanco')&&strcmp(area,'v1_1')
                xlim([0 38]);
            end
            %         if strcmp(animal,'jack')&&strcmp(area,'v1_1')
            %             xlim([0 23]);
            %         end
            if strcmp(animal,'blanco')&&strcmp(area,'v1_1')
                xlim([0 18]);
            end
            if strcmp(animal,'jack')&&strcmp(area,'v1_1')
                xlim([0 23]);
            end
            
            if ~strcmp(area,'v4_2')
                if animalInd==1&&areaInd==1&&sampleInd==1
                    if roving==0
                    pc_condcoefFig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
                    elseif roving==1
                    pc_condcoefFig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.5, 0.8]);
                    end
                else
                    figure(pc_condcoefFig);
                end
                subplot(length(areaTexts),2,animalInd+2*(plotAreaInd-1));
                [rperfcoeff(animalInd+2*(plotAreaInd-1),1:2) pperfcpef(animalInd+2*(plotAreaInd-1),1:2)]=bj_plot_coef_expo_perf(coefperf{animalInd+2*(plotAreaInd-1)}(:,1),sampleContrast,testContrast,ind);
            end
            
            %to check values of coefficients:
            %         figure;plot(1:14,coefperf{animalInd+2*(plotAreaInd-1)}(:,1));
            %         figure;plot(1:14,coefperf{animalInd+2*(plotAreaInd-1)}(:,2));
            %         figure;plot(1:14,coefperf{animalInd+2*(plotAreaInd-1)}(:,3));
            
            %to create test data for exp curve:
            %         xtest=1:length(testContrast);
            %         figure
            %         for testind=1:length(testContrast)
            %             ytest(testind)=exp(-testind);
            %         end
            %         plot(xtest,ytest)
            %         figure
            %         for testind=1:length(testContrast)
            %             ytest(testind)=-exp(-testind);
            %         end
            %         plot(xtest,ytest)
            
            % %
            % figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
            % for i=1:numconds
            %     errorbar(allMeanPerf(:,1),allMeanPerf(:,8+i+numconds),allMeanPerf(:,8+i+numconds*2),'Color',[1/i 1/i 1/i],'LineStyle','none','Marker','o','MarkerFaceColor',[1/i 1/i 1/i],'MarkerEdgeColor',[1/i 1/i 1/i]);hold on
            % end
            % for i=1:numconds
            %     plot(allMeanPerf(:,1),allMeanPerf(:,8+i+numconds),'Color',[1/i 1/i 1/i],'LineStyle','none','Marker','o','MarkerFaceColor',[1/i 1/i 1/i],'MarkerEdgeColor',[1/i 1/i 1/i]);hold on
            % end
            
            %PSE
            figure(pcpseslFig)
            subplot(length(areaTexts),3,3+(plotAreaInd-1)*3);
            animalMarkerFaceColors={'none' 'k'};
            if strcmp(area,'v4_2')
                plot(3+maxV4_1(animalInd):size(allMeanPerf,1)+maxV4_1(animalInd)+2,allMeanPerf(:,17),'ok','MarkerFaceColor',animalMarkerFaceColors{animalInd});
                if strcmp(animal,'jack')
                    xlim([0 38]);
                    ylim([26 40]);
                end
            else
                plot(1:size(allMeanPerf,1),allMeanPerf(:,3+length(testContrast)),'ok','MarkerFaceColor',animalMarkerFaceColors{animalInd});
                hold on
                if strcmp(area,'V4_2')%horizontal Gabor
                    animalMarkerFaceColors={'none' 'b'};
                    plot(size(allMeanPerf,1),allMeanPerf(end,3+length(testContrast)),'Marker','o','Color','b','LineStyle','none','MarkerFaceColor',animalMarkerFaceColors{animalInd});
                end
                [chiselinearpse(animalInd+2*(plotAreaInd-1)) chisepse(animalInd+2*(plotAreaInd-1))]=bj_linearexpo_fitting(testContrast,allMeanPerf(:,3+length(testContrast)),0,0,'ROC');
                if strcmp(area,'v4_1')&&strcmp(animal,'jack')
                    xlim([0 31]);
                    ylim([25 55]);
                end
                if strcmp(area,'v1_1')&&strcmp(animal,'jack')
                    xlim([0 24]);
                    ylim([20 34]);
                end
            end
            
            %slope
            figure(pcpseslFig)
            subplot(length(areaTexts),3,2+(plotAreaInd-1)*3);
            animalMarkerFaceColors={'none' 'k'};
            if strcmp(area,'v4_2')
                plot(3+maxV4_1(animalInd):size(allMeanPerf,1)+maxV4_1(animalInd)+2,allMeanPerf(:,4+length(testContrast)),'ok','MarkerFaceColor',animalMarkerFaceColors{animalInd});
                if strcmp(animal,'jack')
                    xlim([0 38]);
                    ylim([1 11]);
                end
            else
                plot(1:size(allMeanPerf,1),allMeanPerf(:,4+length(testContrast)),'ok','MarkerFaceColor',animalMarkerFaceColors{animalInd});
                hold on
                if strcmp(area,'V4_2')%horizontal Gabor
                    animalMarkerFaceColors={'none' 'b'};
                    plot(size(allMeanPerf,1),allMeanPerf(end,4+length(testContrast)),'Marker','o','Color','b','LineStyle','none','MarkerFaceColor',animalMarkerFaceColors{animalInd});
                end
                [chiselinearsl(animalInd+2*(plotAreaInd-1)) chisesl(animalInd+2*(plotAreaInd-1))]=bj_linearexpo_fitting(testContrast,allMeanPerf(:,4+length(testContrast)),0,0,'ROC');
                if strcmp(area,'v4_1')&&strcmp(animal,'jack')
                    xlim([0 31]);
                    ylim([2 14]);
                end
                if strcmp(area,'v1_1')&&strcmp(animal,'jack')
                    xlim([0 24]);
                    ylim([0 7]);
                end
            end
            
            % slope and PSE, for first and last 30% of trials within session
            figure(pcpsesl30Fig)
            subplot(length(areaTexts),3,2+(plotAreaInd-1)*3);
            animalMarkerFaceColors={'none' 'k'};
            if strcmp(area,'v4_2')
                plot(3+maxV4_1(animalInd):size(allMeanPerf,1)+maxV4_1(animalInd)+2,allMeanPerf(:,29+numconds*6),'Marker','o','Color','k','LineStyle','none','MarkerFaceColor',animalMarkerFaceColors{animalInd});
            else
                plot(1:size(allMeanPerf,1),allMeanPerf(:,29+numconds*6),'Marker','o','Color','k','LineStyle','none','MarkerFaceColor',animalMarkerFaceColors{animalInd});
            end
            hold on
            animalMarkerFaceColors={'none' 'r'};
            if strcmp(area,'v4_2')
                plot(3+maxV4_1(animalInd):size(allMeanPerf,1)+maxV4_1(animalInd)+2,allMeanPerf(:,36+numconds*7),'Marker','o','Color','r','LineStyle','none','MarkerFaceColor',animalMarkerFaceColors{animalInd});
            else
                plot(1:size(allMeanPerf,1),allMeanPerf(:,36+numconds*7),'Marker','o','Color','r','LineStyle','none','MarkerFaceColor',animalMarkerFaceColors{animalInd});
                [hperfIndpair(plotAreaInd+length(areaTexts)*(animalInd-1),2),pperfIndpair(plotAreaInd+length(areaTexts)*(animalInd-1),2)]=ttest(allMeanPerf(:,29+numconds*6),allMeanPerf(:,36+numconds*7));
            end
            subplot(length(areaTexts),3,3+(plotAreaInd-1)*3);
            animalMarkerFaceColors={'none' 'k'};
            if strcmp(area,'v4_2')
                plot(3+maxV4_1(animalInd):size(allMeanPerf,1)+maxV4_1(animalInd)+2,allMeanPerf(:,28+numconds*6),'Marker','o','Color','k','LineStyle','none','MarkerFaceColor',animalMarkerFaceColors{animalInd});
            else
                plot(1:size(allMeanPerf,1),allMeanPerf(:,28+numconds*6),'Marker','o','Color','k','LineStyle','none','MarkerFaceColor',animalMarkerFaceColors{animalInd});
            end
            hold on
            animalMarkerFaceColors={'none' 'r'};
            if strcmp(area,'v4_2')
                plot(3+maxV4_1(animalInd):size(allMeanPerf,1)+maxV4_1(animalInd)+2,allMeanPerf(:,35+numconds*7),'Marker','o','Color','r','LineStyle','none','MarkerFaceColor',animalMarkerFaceColors{animalInd});
            else
                plot(1:size(allMeanPerf,1),allMeanPerf(:,35+numconds*7),'Marker','o','Color','r','LineStyle','none','MarkerFaceColor',animalMarkerFaceColors{animalInd});
                [hperfIndpair(plotAreaInd+length(areaTexts)*(animalInd-1),3),pperfIndpair(plotAreaInd+length(areaTexts)*(animalInd-1),3)]=ttest(allMeanPerf(:,28+numconds*6),allMeanPerf(:,35+numconds*7));
            end
            
            % index of difference in slope and PSE,
            % between first and last 30% of trials within session
            figure(pcpsesl30indFig)
            subplot(length(areaTexts),3,2+(plotAreaInd-1)*3);
            animalMarkerFaceColors={'none' 'k'};
            slInd=(allMeanPerf(:,36+numconds*7)-allMeanPerf(:,29+numconds*6))./allMeanPerf(:,29+numconds*6);
            if strcmp(area,'v4_2')
                plot(3+maxV4_1(animalInd):size(allMeanPerf,1)+maxV4_1(animalInd)+2,slInd,'Marker','o','Color','k','LineStyle','none','MarkerFaceColor',animalMarkerFaceColors{animalInd});
            else
                plot(1:size(allMeanPerf,1),slInd,'Marker','o','Color','k','LineStyle','none','MarkerFaceColor',animalMarkerFaceColors{animalInd});
                [hperfInd(plotAreaInd+2*(animalInd-1),2),pperfInd(plotAreaInd+length(areaTexts)*(animalInd-1),2)]=ttest(slInd);
            end
            hold on
            subplot(length(areaTexts),3,3+(plotAreaInd-1)*3);
            animalMarkerFaceColors={'none' 'k'};
            pseInd=(allMeanPerf(:,35+numconds*7)-allMeanPerf(:,28+numconds*6))./allMeanPerf(:,28+numconds*6);
            if strcmp(area,'v4_2')
                plot(3+maxV4_1(animalInd):size(allMeanPerf,1)+maxV4_1(animalInd)+2,pseInd,'Marker','o','Color','k','LineStyle','none','MarkerFaceColor',animalMarkerFaceColors{animalInd});
            else
                plot(1:size(allMeanPerf,1),pseInd,'Marker','o','Color','k','LineStyle','none','MarkerFaceColor',animalMarkerFaceColors{animalInd});
                [hperfInd(plotAreaInd+2*(animalInd-1),3),pperfInd(plotAreaInd+length(areaTexts)*(animalInd-1),3)]=ttest(pseInd);
            end
            hold on
        end
    end
end

%compare sizes of adjusted R-square values, between linear and polynomial fits, for PC divided by test contrast
aRSlinexpo1=[];
aRSlinexpo2=[];
for animalInd=1:length(animals)
    for areaInd=1:2%hard-coded for V4_1 and V1_1
        aRSlinexpo1=[aRSlinexpo1;aRSlinear{animalInd+2*(plotAreaInd-1)}(:,:)];
        aRSlinexpo2=[aRSlinexpo2;aRS{animalInd+2*(plotAreaInd-1)}(:,:)];
    end
end
aRSlinexpo=[aRSlinexpo1 aRSlinexpo2];
[h,p,ci,stats]=ttest(aRSlinexpo(:,1),aRSlinexpo(:,2))
mean(aRSlinexpo(:,1)-aRSlinexpo(:,2));
figure
diff=aRSlinexpo(:,1)-aRSlinexpo(:,2);hold on
plot(sort(diff));
plot([0 56],[0 0],'k:');
sum(diff<0)/length(diff)

if roving==0
    % %proportion correct
    figure(pcpseslFig)
    subplot(length(areaTexts),3,1);
    ylim([50 90]);
    % ylim([0.55 0.9]);
    subplot(length(areaTexts),3,3);
    ylim([27 40]);
    subplot(length(areaTexts),3,4);
    ylim([60 95]);
    % ylim([0.65 0.95]);
    xlim([0 24]);
    
    % %proportion correct, slope, PSE: first and last 30% of trials within each
    % session
    figure(pcpsesl30Fig)
    subplot(length(areaTexts),3,1);
    ylim([0.5 0.9]);
    xlim([0 38]);
    subplot(length(areaTexts),3,2);
    ylim([0 20]);
    xlim([0 38]);
    subplot(length(areaTexts),3,3);
    ylim([26 38]);
    xlim([0 38]);
    subplot(length(areaTexts),3,4);
    ylim([0.65 0.95]);
    xlim([0 24]);
    subplot(length(areaTexts),3,5);
    ylim([0 10]);
    xlim([0 24]);
    subplot(length(areaTexts),3,6);
    ylim([16 36]);
    xlim([0 24]);
    
    % %proportion correct, slope, PSE: first and last 30% of trials within each
    % session, index of difference
    figure(pcpsesl30indFig)
    subplot(length(areaTexts),3,1);
    line([0 38],[0 0],'LineStyle',':','Color','k');
    ylim([-0.2 0.2]);
    xlim([0 38]);
    subplot(length(areaTexts),3,2);
    line([0 38],[0 0],'LineStyle',':','Color','k');
    ylim([-4 4]);
    xlim([0 38]);
    subplot(length(areaTexts),3,3);
    line([0 38],[0 0],'LineStyle',':','Color','k');
    ylim([-0.1 0.1]);
    xlim([0 38]);
    subplot(length(areaTexts),3,4);
    line([0 38],[0 0],'LineStyle',':','Color','k');
    ylim([-0.25 0.25]);
    xlim([0 24]);
    subplot(length(areaTexts),3,5);
    line([0 38],[0 0],'LineStyle',':','Color','k');
    ylim([-4 4]);
    xlim([0 24]);
    subplot(length(areaTexts),3,6);
    line([0 38],[0 0],'LineStyle',':','Color','k');
    ylim([-0.8 0.8]);
    xlim([0 24]);
    
    % figure(pcFig)
    % subplot(length(areaTexts),2,1);
    % xlim([0 38]);
    % ylim([0.55 1.05]);
    % subplot(length(areaTexts),2,2);
    % xlim([0 35]);
    % ylim([0.4 1.05]);
    % subplot(length(areaTexts),2,3);
    % xlim([0 18]);
    % ylim([0.5 1.05]);
    % subplot(length(areaTexts),2,4);
    % xlim([0 23]);
    % ylim([0.5 1.05]);
    
    figure(pc_condFig)
    subplot(length(areaTexts),2,1);
    % ylim([-0.1 1.1]);
    xlim([0 30]);
    ylim([-5 105]);
    subplot(length(areaTexts),2,2);
    % ylim([-0.1 1.1]);
    xlim([0 34]);
    ylim([-5 105]);
    subplot(length(areaTexts),2,3);
    % ylim([-0.1 1.1]);
    ylim([-5 105]);
    subplot(length(areaTexts),2,4);
    % ylim([-0.1 1.1]);
    ylim([-5 105]);
    
    % figure(pc_condFig)
    % subplot(length(areaTexts),2,1);
    % ylim([-20 105]);
    % subplot(length(areaTexts),2,2);
    % ylim([-20 105]);
    % subplot(length(areaTexts),2,3);
    % ylim([-20 105]);
    % subplot(length(areaTexts),2,4);
    % ylim([-20 105]);
    
    figure(pc_condcoefFig)
    subplot(length(areaTexts),2,1);
    ylim([0 60]);
    subplot(length(areaTexts),2,2);
    ylim([0 60]);
    subplot(length(areaTexts),2,3);
    ylim([0 60]);
    subplot(length(areaTexts),2,4);
    ylim([0 60]);
    
    %mean RT with error bars
    figure(rtFig)
    subplot(length(areaTexts),2,1);
    xlim([0 30]);
    ylim([80 300]);
    subplot(length(areaTexts),2,2);
    xlim([0 35]);
    ylim([110 230]);
    subplot(length(areaTexts),2,3);
    xlim([0 18]);
    ylim([100 230]);
    subplot(length(areaTexts),2,4);
    xlim([0 23]);
    ylim([100 240]);
    
    %difference between mean RTs for correct and error trials
    figure(diffmeanRTFig)
    subplot(1,length(areaTexts),1);
    xlim([0 38]);
    ylim([-0.05 0.2]);
    line([0 38],[0 0],'LineStyle',':','Color','k');
    subplot(1,length(areaTexts),2);
    xlim([0 23]);
    ylim([-0.05 0.2]);
    line([0 23],[0 0],'LineStyle',':','Color','k');
    
    
    %mean RT with error bars per cond
    figure(rt_condFig)
    subplot(length(areaTexts),2,1);
    xlim([0 30]);
    ylim([20 320]);
    subplot(length(areaTexts),2,2);
    xlim([0 35]);
    ylim([100 250]);
    subplot(length(areaTexts),2,3);
    xlim([0 18]);
    ylim([80 240]);
    subplot(length(areaTexts),2,4);
    xlim([0 23]);
    ylim([80 250]);
    
    %mean error trial RT with error bars per cond
    figure(rterror_condFig)
    subplot(length(areaTexts),2,1);
    xlim([0 30]);
    ylim([10 470]);
    subplot(length(areaTexts),2,2);
    xlim([0 35]);
    ylim([70 330]);
    subplot(length(areaTexts),2,3);
    xlim([0 18]);
    ylim([0 360]);
    subplot(length(areaTexts),2,4);
    xlim([0 23]);
    ylim([100 360]);
    
    %difference between correct and error trial RT
    figure(rtcediffFig)
    subplot(length(areaTexts),2,1);
    xlim([0 30]);
    ylim([-10 10]);
    subplot(length(areaTexts),2,2);
    xlim([0 35]);
    ylim([-10 10]);
    subplot(length(areaTexts),2,3);
    xlim([0 18]);
    ylim([-10 10]);
    subplot(length(areaTexts),2,4);
    xlim([0 23]);
    ylim([-10 10]);
elseif roving==1    
    if strcmp(area,'v1_2_1')
    % %proportion correct
    figure(pcpseslFig)
%     subplot(length(areaTexts),3,1);
%     ylim([50 90]);
    subplot(length(areaTexts),3,2);
    ylim([1 5]);
%     % ylim([0.55 0.9]);
%     subplot(length(areaTexts),3,3);
%     ylim([27 40]);
%     subplot(length(areaTexts),3,4);
%     ylim([60 95]);
%     % ylim([0.65 0.95]);
%     xlim([0 24]);
    
    % %proportion correct, slope, PSE: first and last 30% of trials within each
    % session
    figure(pcpsesl30Fig)
%     subplot(length(areaTexts),3,1);
%     ylim([0.5 0.9]);
%     xlim([0 38]);
    subplot(length(areaTexts),3,2);
    ylim([0 8]);
%     xlim([0 38]);
    subplot(length(areaTexts),3,3);
    ylim([20 35]);
%     xlim([0 38]);
%     subplot(length(areaTexts),3,4);
%     ylim([0.65 0.95]);
%     xlim([0 24]);
    subplot(length(areaTexts),3,5);
    ylim([0 20]);
%     xlim([0 24]);
%     subplot(length(areaTexts),3,6);
%     ylim([16 36]);
%     xlim([0 24]);
    subplot(length(areaTexts),3,8);
    ylim([0 10]);
    subplot(length(areaTexts),3,9);
    ylim([25 40]);
    
    % %proportion correct, slope, PSE: first and last 30% of trials within each
    % session, index of difference
    figure(pcpsesl30indFig)
    subplot(length(areaTexts),3,1);
    line([0 38],[0 0],'LineStyle',':','Color','k');
%     ylim([-0.2 0.2]);
    xlim([0 37]);
    subplot(length(areaTexts),3,2);
    line([0 38],[0 0],'LineStyle',':','Color','k');
    ylim([-2 5]);
    xlim([0 37]);
    subplot(length(areaTexts),3,3);
    line([0 38],[0 0],'LineStyle',':','Color','k');
%     ylim([-0.1 0.1]);
    xlim([0 37]);
    subplot(length(areaTexts),3,4);
    line([0 38],[0 0],'LineStyle',':','Color','k');
%     ylim([-0.25 0.25]);
    xlim([0 37]);
    subplot(length(areaTexts),3,5);
    line([0 38],[0 0],'LineStyle',':','Color','k');
    ylim([-2 4]);
    xlim([0 37]);
    subplot(length(areaTexts),3,6);
    line([0 38],[0 0],'LineStyle',':','Color','k');
%     ylim([-0.8 0.8]);
    xlim([0 37]);
    subplot(length(areaTexts),3,7);
    line([0 38],[0 0],'LineStyle',':','Color','k');
%     ylim([-0.8 0.8]);
    xlim([0 37]);
    subplot(length(areaTexts),3,8);
    line([0 38],[0 0],'LineStyle',':','Color','k');
    ylim([-2 2]);
    xlim([0 37]);
    subplot(length(areaTexts),3,9);
    line([0 38],[0 0],'LineStyle',':','Color','k');
    ylim([-0.2 0.2]);
    xlim([0 37]);
    
    figure(pc_condFig)
    subplot(length(areaTexts),2,1);
    % ylim([-0.1 1.1]);
%     xlim([0 37]);
    ylim([-5 105]);
    subplot(length(areaTexts),2,2);
    % ylim([-0.1 1.1]);
    xlim([0 17]);
    ylim([-5 105]);
    subplot(length(areaTexts),2,3);
    % ylim([-0.1 1.1]);
    ylim([-5 105]);
    subplot(length(areaTexts),2,4);
    % ylim([-0.1 1.1]);
    xlim([0 17]);
    ylim([-5 105]);
    subplot(length(areaTexts),2,5);
    % ylim([-0.1 1.1]);
    ylim([-5 105]);
    subplot(length(areaTexts),2,6);
    % ylim([-0.1 1.1]);
    xlim([0 17]);
    ylim([-5 105]);
    
   
    figure(pc_condcoefFig)
%     subplot(length(areaTexts),2,1);
%     ylim([0 60]);
%     subplot(length(areaTexts),2,2);
%     ylim([0 60]);
%     subplot(length(areaTexts),2,3);
%     ylim([0 60]);
%     subplot(length(areaTexts),2,4);
%     ylim([0 60]);
    
    %mean RT with error bars
    figure(rtFig)
    subplot(length(areaTexts),2,1);
%     xlim([0 40]);
    ylim([40 250]);
    subplot(length(areaTexts),2,2);
    xlim([0 17]);
    ylim([130 210]);
    subplot(length(areaTexts),2,3);
%     xlim([0 18]);
    ylim([40 250]);
    subplot(length(areaTexts),2,4);
    xlim([0 17]);
    ylim([130 220]);
    subplot(length(areaTexts),2,6);
    xlim([0 17]);
    ylim([130 210]);
    
    %difference between mean RTs for correct and error trials
    figure(diffmeanRTFig)
    subplot(1,length(areaTexts),1);
    xlim([0 38]);
%     ylim([-0.05 0.2]);
    line([0 38],[0 0],'LineStyle',':','Color','k');
    subplot(1,length(areaTexts),2);
%     xlim([0 23]);
    ylim([-0.05 0.25]);
    line([0 35],[0 0],'LineStyle',':','Color','k');
    subplot(1,length(areaTexts),3);
%     xlim([0 23]);
%     ylim([-0.05 0.2]);
    line([0 35],[0 0],'LineStyle',':','Color','k');
    
    
    %mean RT with error bars per cond
    figure(rt_condFig)
%     subplot(length(areaTexts),2,1);
%     xlim([0 38]);
%     ylim([20 320]);
    subplot(length(areaTexts),2,2);
    xlim([0 17]);
    ylim([100 250]);
    subplot(length(areaTexts),2,3);
%     xlim([0 18]);
    ylim([0 250]);
    subplot(length(areaTexts),2,4);
    xlim([0 17]);
    ylim([100 250]);
    subplot(length(areaTexts),2,5);
%     xlim([0 17]);
    ylim([0 300]);
    subplot(length(areaTexts),2,6);
    xlim([0 17]);
    ylim([100 250]);
    
    %mean error trial RT with error bars per cond
    figure(rterror_condFig)
    subplot(length(areaTexts),2,1);
%     xlim([0 38]);
    ylim([0 300]);
    subplot(length(areaTexts),2,2);
    xlim([0 17]);
    ylim([130 240]);
%     subplot(length(areaTexts),2,3);
%     xlim([0 18]);
%     ylim([0 360]);
%     subplot(length(areaTexts),2,4);
%     xlim([0 23]);
%     ylim([100 360]);
    subplot(length(areaTexts),2,5);
%     xlim([0 18]);
    ylim([0 250]);
    subplot(length(areaTexts),2,6);
    xlim([0 17]);
    ylim([110 280]);
    
    %difference between correct and error trial RT
    figure(rtcediffFig)
    subplot(length(areaTexts),2,1);
%     xlim([0 38]);
    ylim([-10 10]);
    subplot(length(areaTexts),2,2);
%     xlim([0 35]);
    ylim([-10 10]);
    subplot(length(areaTexts),2,3);
%     xlim([0 18]);
    ylim([-10 10]);
    subplot(length(areaTexts),2,4);
%     xlim([0 23]);
    ylim([-10 10]);
    subplot(length(areaTexts),2,5);
%     xlim([0 18]);
    ylim([-10 10]);
    subplot(length(areaTexts),2,6);
%     xlim([0 23]);
    ylim([-10 10]);
    elseif strcmp(area,'v1_2_2')
    % %proportion correct
    figure(pcpseslFig)
    subplot(length(areaTexts),3,1);
    xlim([0 24]);
%     ylim([50 90]);
%     % ylim([0.55 0.9]);
    subplot(length(areaTexts),3,2);
    xlim([0 24]);
    subplot(length(areaTexts),3,3);
    xlim([0 24]);
%     ylim([27 40]);
    subplot(length(areaTexts),3,4);
    xlim([0 24]);
%     ylim([60 95]);
%     % ylim([0.65 0.95]);
    subplot(length(areaTexts),3,5);
    xlim([0 24]);
    subplot(length(areaTexts),3,6);
    xlim([0 24]);
    ylim([20 50]);
    subplot(length(areaTexts),3,7);
    xlim([0 24]);
    subplot(length(areaTexts),3,8);
    xlim([0 24]);
    subplot(length(areaTexts),3,9);
    xlim([0 24]);
    
    % %proportion correct, slope, PSE: first and last 30% of trials within each
    % session
    figure(pcpsesl30Fig)
    subplot(length(areaTexts),3,1);
%     ylim([0.5 0.9]);
    xlim([0 24]);
    subplot(length(areaTexts),3,2);
    ylim([0 15]);
    xlim([0 24]);
    subplot(length(areaTexts),3,3);
%     ylim([20 35]);
    xlim([0 24]);
    subplot(length(areaTexts),3,4);
%     ylim([0.65 0.95]);
    xlim([0 24]);
    subplot(length(areaTexts),3,5);
    ylim([0 10]);
    xlim([0 24]);
    subplot(length(areaTexts),3,6);
    ylim([20 50]);
    xlim([0 24]);
    subplot(length(areaTexts),3,7);
    xlim([0 24]);
    subplot(length(areaTexts),3,8);
    ylim([0 17]);
    xlim([0 24]);
    subplot(length(areaTexts),3,9);
    ylim([25 65]);
    xlim([0 24]);
    
    % %proportion correct, slope, PSE: first and last 30% of trials within each
    % session, index of difference
    figure(pcpsesl30indFig)
    subplot(length(areaTexts),3,1);
    line([0 24],[0 0],'LineStyle',':','Color','k');
%     ylim([-0.2 0.2]);
    xlim([0 24]);
    subplot(length(areaTexts),3,2);
    line([0 24],[0 0],'LineStyle',':','Color','k');
    ylim([-2 4]);
    xlim([0 24]);
    subplot(length(areaTexts),3,3);
    line([0 24],[0 0],'LineStyle',':','Color','k');
%     ylim([-0.1 0.1]);
    xlim([0 24]);
    subplot(length(areaTexts),3,4);
    line([0 24],[0 0],'LineStyle',':','Color','k');
%     ylim([-0.25 0.25]);
    xlim([0 24]);
    subplot(length(areaTexts),3,5);
    line([0 24],[0 0],'LineStyle',':','Color','k');
    ylim([-2 2]);
    xlim([0 24]);
    subplot(length(areaTexts),3,6);
    line([0 24],[0 0],'LineStyle',':','Color','k');
%     ylim([-0.8 0.8]);
    xlim([0 24]);
    subplot(length(areaTexts),3,7);
    line([0 24],[0 0],'LineStyle',':','Color','k');
%     ylim([-0.8 0.8]);
%     xlim([0 37]);
    subplot(length(areaTexts),3,8);
    line([0 24],[0 0],'LineStyle',':','Color','k');
    ylim([-2 4]);
    xlim([0 24]);
    subplot(length(areaTexts),3,9);
    line([0 24],[0 0],'LineStyle',':','Color','k');
    ylim([-0.2 0.3]);
    xlim([0 24]);
    
    figure(pc_condFig)
    subplot(length(areaTexts),2,1);
    % ylim([-0.1 1.1]);
    xlim([0 16]);
    ylim([-5 105]);
    subplot(length(areaTexts),2,2);
    % ylim([-0.1 1.1]);
    xlim([0 23]);
    ylim([-5 105]);
    subplot(length(areaTexts),2,3);
    xlim([0 16]);
    % ylim([-0.1 1.1]);
    ylim([-5 105]);
    subplot(length(areaTexts),2,4);
    % ylim([-0.1 1.1]);
    xlim([0 23]);
    ylim([-5 105]);
    subplot(length(areaTexts),2,5);
    xlim([0 16]);
    % ylim([-0.1 1.1]);
    ylim([-5 105]);
    subplot(length(areaTexts),2,6);
    % ylim([-0.1 1.1]);
    xlim([0 23]);
    ylim([-5 105]);
    
   
    figure(pc_condcoefFig)
%     subplot(length(areaTexts),2,1);
%     ylim([0 60]);
    subplot(length(areaTexts),2,2);
    ylim([-50 100]);
%     subplot(length(areaTexts),2,3);
%     ylim([0 60]);
    subplot(length(areaTexts),2,4);
    ylim([-10 60]);
    
    %mean RT with error bars
    figure(rtFig)
    subplot(length(areaTexts),2,1);
    xlim([0 16]);
%     ylim([40 250]);
    subplot(length(areaTexts),2,2);
    xlim([0 23]);
%     ylim([130 210]);
    subplot(length(areaTexts),2,3);
    xlim([0 16]);
%     ylim([40 250]);
    subplot(length(areaTexts),2,4);
    xlim([0 23]);
%     ylim([130 220]);
    subplot(length(areaTexts),2,5);
    xlim([0 16]);
    subplot(length(areaTexts),2,6);
    xlim([0 23]);
%     ylim([130 210]);
    
    %difference between mean RTs for correct and error trials
    figure(diffmeanRTFig)
    subplot(1,length(areaTexts),1);
    xlim([0 23]);
%     ylim([-0.05 0.2]);
    line([0 38],[0 0],'LineStyle',':','Color','k');
    subplot(1,length(areaTexts),2);
    xlim([0 23]);
    ylim([-0.05 0.25]);
    line([0 35],[0 0],'LineStyle',':','Color','k');
    subplot(1,length(areaTexts),3);
    xlim([0 23]);
%     ylim([-0.05 0.2]);
    line([0 35],[0 0],'LineStyle',':','Color','k');
    
    
    %mean RT with error bars per cond
    figure(rt_condFig)
    subplot(length(areaTexts),2,1);
%     xlim([0 38]);
    ylim([0 270]);
%     subplot(length(areaTexts),2,2);
%     xlim([0 17]);
%     ylim([100 250]);
%     subplot(length(areaTexts),2,3);
%     xlim([0 18]);
%     ylim([0 250]);
%     subplot(length(areaTexts),2,4);
%     xlim([0 17]);
%     ylim([100 250]);
%     subplot(length(areaTexts),2,5);
%     xlim([0 17]);
%     ylim([0 300]);
%     subplot(length(areaTexts),2,6);
%     xlim([0 17]);
%     ylim([100 250]);
    
    %mean error trial RT with error bars per cond
    figure(rterror_condFig)
    subplot(length(areaTexts),2,1);
%     xlim([0 38]);
%     ylim([0 300]);
    subplot(length(areaTexts),2,2);
    xlim([0 23]);
    ylim([50 250]);
%     subplot(length(areaTexts),2,3);
%     xlim([0 18]);
%     ylim([0 360]);
    subplot(length(areaTexts),2,4);
    xlim([0 23]);
    ylim([70 250]);
    subplot(length(areaTexts),2,5);
%     xlim([0 18]);
%     ylim([0 250]);
    subplot(length(areaTexts),2,6);
    xlim([0 23]);
    ylim([50 300]);
    
    %difference between correct and error trial RT
    figure(rtcediffFig)
    subplot(length(areaTexts),2,1);
%     xlim([0 38]);
    ylim([-10 10]);
    subplot(length(areaTexts),2,2);
%     xlim([0 35]);
    ylim([-10 10]);
    subplot(length(areaTexts),2,3);
%     xlim([0 18]);
    ylim([-10 10]);
    subplot(length(areaTexts),2,4);
%     xlim([0 23]);
    ylim([-10 10]);
    subplot(length(areaTexts),2,5);
%     xlim([0 18]);
    ylim([-10 10]);
    subplot(length(areaTexts),2,6);
%     xlim([0 23]);
    ylim([-10 10]);
    elseif strcmp(area,'v1_2_3')
        figure(pc_condcoefFig)
        ylim([-20 70]);
    end
end

if printFigs==1
    if roving==0
        rovingText='non_roving';
    elseif roving==1
        rovingText=['roving_',area];
    end
    folderPrint=fullfile(rootFolder,'PL','psycho_data','behavioural_figures',rovingText);
    formats=[{'epsc'} {'png'}];
    for i=1:length(formats)
        format=formats{i};
        if roving==0
            figure(psychoCurvesj2sess)
            figName=[folderPrint,'\_jack_psychoj2sess'];
            printtext=sprintf('print -d%s %s -r600',format,figName);
            set(gcf, 'PaperPositionMode', 'auto');
            eval(printtext);
        end
        for animalInd=1:length(animals)
            for areaInd=1:length(areas)
                for sampleInd=1:length(sampleContrasts)
                figure(psychoCurves(animalInd+2*(areaInd-1)))
                figName=[folderPrint,'\_',animals{animalInd},'_psycho_',areas{areaInd},'_',sampleContrasts(sampleInd)];
                printtext=sprintf('print -d%s %s -r600',format,figName);
                set(gcf, 'PaperPositionMode', 'auto');
                eval(printtext);
                end
            end
        end
        figure(pc_condFig)
        figName=[folderPrint,'\_pc_cond'];
        printtext=sprintf('print -d%s %s -r600',format,figName);
        set(gcf, 'PaperPositionMode', 'auto');
        eval(printtext);
        figure(pc_condcoefFig)
        figName=[folderPrint,'\_pc_cond_coeff'];
        printtext=sprintf('print -d%s %s -r600',format,figName);
        set(gcf, 'PaperPositionMode', 'auto');
        eval(printtext);
        figure(rtFig)
        figName=[folderPrint,'\_rt'];
        printtext=sprintf('print -d%s %s -r600',format,figName);
        set(gcf, 'PaperPositionMode', 'auto');
        eval(printtext);
        figure(diffmeanRTFig)
        figName=[folderPrint,'\_meanrt_diff'];
        printtext=sprintf('print -d%s %s -r600',format,figName);
        set(gcf, 'PaperPositionMode', 'auto');
        eval(printtext);
        %     figure(rterrorFig)
        %     figName=[folderPrint,'\_rterr'];
        %     printtext=sprintf('print -d%s %s -r600',format,figName);
        %     set(gcf, 'PaperPositionMode', 'auto');
        %eval(printtext);
        figure(rt_condFig)
        figName=[folderPrint,'\_rt_cond'];
        printtext=sprintf('print -d%s %s -r600',format,figName);
        set(gcf, 'PaperPositionMode', 'auto');
        eval(printtext);
        figure(rterror_condFig)
        figName=[folderPrint,'\_rterror_cond'];
        printtext=sprintf('print -d%s %s -r600',format,figName);
        set(gcf, 'PaperPositionMode', 'auto');
        eval(printtext);
        figure(rtcediffFig)
        %         figName=[folderPrint,'\_rterrordiff_cond'];
        %         printtext=sprintf('print -d%s %s -r600',format,figName);
        %         set(gcf, 'PaperPositionMode', 'auto');
        %         eval(printtext);
        figure(RT30)
        figName=[folderPrint,'\_rt30'];
        printtext=sprintf('print -d%s %s -r600',format,figName);
        set(gcf, 'PaperPositionMode', 'auto');
        eval(printtext);
        figure(RT30ind)
        figName=[folderPrint,'\_rt30_diff'];
        printtext=sprintf('print -d%s %s -r600',format,figName);
        set(gcf, 'PaperPositionMode', 'auto');
        eval(printtext);
        figure(pcpseslFig)
        figName=[folderPrint,'\_pcpsesl'];
        printtext=sprintf('print -d%s %s -r600',format,figName);
        set(gcf, 'PaperPositionMode', 'auto');
        eval(printtext);
        figure(pcpsesl30Fig)
        figName=[folderPrint,'\_pcpsesl30'];
        printtext=sprintf('print -d%s %s -r600',format,figName);
        set(gcf, 'PaperPositionMode', 'auto');
        eval(printtext);
        figure(pcpsesl30indFig)
        figName=[folderPrint,'\_pcpsesl30ind'];
        printtext=sprintf('print -d%s %s -r600',format,figName);
        set(gcf, 'PaperPositionMode', 'auto');
        eval(printtext);
    end
end
%refer to tables in
%V:\thielelab\Groups\ThieleGroup\monkey_data\blanco\_grid\blanco jack psychometric analysis.xlsx
pRTceANOVA%Two-way ANOVA on RTs- X1: correct vs error trial; X2: lower test contrast vs higher test contrast
pperfIndpair%paired t-test between early and late 30% of trials
pperfInd%t-test to see if mean of difference between early and late values is significant different from 0
pRTind%t-test for differences between early vs late trials mean RT for correct trials
rmeanRT%correlation between RT and time (conducted separately for correct and error trials), p-vals
pmeanRT%correlation between RT and time (conducted separately for correct and error trials), R-vals
close all