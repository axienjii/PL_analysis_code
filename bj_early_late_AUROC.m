function bj_early_late_AUROC
if nargin<3||isempty(animals)
    animals=[{'blanco'} {'jack'}];
    % animals={'blanco'};
end
roving=0;
if nargin<4||isempty(areas)
    areas=[{'v4_1'} {'v4_2'} {'v1_1'} {'v1_2'}];
    areas=[{'v4_1'} {'v1_1'} {'v1_2_1'} {'v1_2_2'} {'v1_2_3'}];
    if roving==0
        areas=[{'v4_1'} {'v1_1'}];
    elseif roving==1
        areas=[{'v1_2_1'} {'v1_2_2'}];
        areas=[{'v1_2_2'}];
    end
end
calcStats=0;
if calcStats==1
    fig1=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.3, 0.4]); %
    set(fig1, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
    fig2=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.3, 0.8]); %
    set(fig2, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
    allParams={[]};
    sessArr={[]};
    animalArr={[]};
    areaArr={[]};
    statsTable=[];
    for animalInd=1:length(animals)
        animal=animals{animalInd};
        for areaInd=1:length(areas)
            area=areas{areaInd};
            [sampleContrasts testContrasts]=area_metadata(area);
            channels=main_channels(animal,area);
            sessionNums=main_raw_sessions_final(animal,area,[],0);
            for sampleContrastsInd=1:length(sampleContrasts)
                sampleContrast=sampleContrasts(sampleContrastsInd);
                testContrast=testContrasts(sampleContrastsInd,:);
                loadText=['load F:\PL\ROC_zero_one\',animal,'\new_vs_old_sglroc3acrosschannels\cumulative_ROCs_old_new_',area,'_',num2str(sampleContrast),'_cutoff10.mat'];
                eval(loadText)
                %compare AUROC values directly, between earlyand late sessions
                %(first and last 30% of sessions)
                [h p stats]=ttest(threshold82higher(1:floor(length(threshold82higher)*0.3)),threshold82higher(end-floor(length(threshold82higher)*0.3)+1:end));
                [h p stats]=ttest(threshold82lower(1:floor(length(threshold82higher)*0.3)),threshold82lower(end-floor(length(threshold82higher)*0.3)+1:end));
                rocEarly=all_rocvals(1:floor(length(threshold82higher)*0.3),2:end);
                rocLate=all_rocvals(end-floor(length(threshold82higher)*0.3)+1:end,2:end);
                figure(fig1)
                subplot(2,2,animalInd+2*(areaInd-1));
                for condInd=1:size(rocEarly,2)
                    [h p stats]=ttest2(rocEarly(:,condInd),rocLate(:,condInd));%unpaired t-test
                    if p<.05
                        MarkerType='o';
                    else
                        MarkerType='.';
                    end
                    plot(condInd,mean(rocEarly(:,condInd)),'LineStyle','none','Marker',MarkerType,'Color','k');hold on
                    plot(condInd,mean(rocLate(:,condInd)),'LineStyle','none','Marker',MarkerType,'Color','r');hold on
                    errorbar(condInd,mean(rocEarly(:,condInd)),std(rocEarly(:,condInd)),'LineStyle','none','Marker',MarkerType,'Color','k');hold on
                    errorbar(condInd,mean(rocLate(:,condInd)),std(rocLate(:,condInd)),'LineStyle','none','Marker',MarkerType,'Color','r');hold on
                end
                ylim([0 1]);
                if animalInd+2*(areaInd-1)==1
                    xlabel('contrast (%)');
                    ylabel('AUROC');
                    title('Monkey 1');
                elseif animalInd+2*(areaInd-1)==2
                    title('Monkey 2');
                end
                %compare slope, PNE, min, max, & lower & upper thresholds
                figure(fig2)
                params=[slopeNeuroNew;PNENew;minRateNew;maxRateNew;threshold82lower;threshold82higher];
                allParams{sampleContrastsInd}=[allParams{sampleContrastsInd} params(:,1:floor(size(params,2)*0.3)) params(:,end-floor(size(params,2)*0.3)+1:end)];
                sessArr{sampleContrastsInd}=[sessArr{sampleContrastsInd} zeros(size(params,1),floor(size(params,2)*0.3))+1 zeros(size(params,1),floor(size(params,2)*0.3))+2];%early and late sessions
                animalArr{sampleContrastsInd}=[animalArr{sampleContrastsInd} zeros(size(params,1),2*floor(size(params,2)*0.3))+animalInd];
                areaArr{sampleContrastsInd}=[areaArr{sampleContrastsInd} zeros(size(params,1),2*floor(size(params,2)*0.3))+areaInd];
                paramText=[{'slope'} {'PNE'} {'min'} {'max'} {'upper threshold'} {'lower threshold'}];
                for paramInd=1:size(params,1)
                    [h(paramInd) p(paramInd) stats]=ttest2(params(paramInd,1:floor(size(params,2)*0.3)),params(paramInd,end-floor(size(params,2)*0.3)+1:end));%unpaired t-test
                    subplot(size(params,1),4,animalInd+2*(areaInd-1)+4*(paramInd-1));
                    if p(paramInd)<.05
                        MarkerType='o';
                    else
                        MarkerType='.';
                    end
                    %                 plot(1,params(paramInd,1:floor(size(params,2)*0.3)),'LineStyle','none','Marker',MarkerType,'Color','k');hold on%early
                    %                 plot(2,params(paramInd,end-floor(size(params,2)*0.3)+1:end),'LineStyle','none','Marker',MarkerType,'Color','r');hold on%late
                    plot(1,mean(params(paramInd,1:floor(size(params,2)*0.3))),'LineStyle','none','Marker',MarkerType,'Color','k');hold on
                    plot(2,mean(params(paramInd,end-floor(size(params,2)*0.3)+1:end)),'LineStyle','none','Marker',MarkerType,'Color','r');hold on
                    errorbar(1,mean(params(paramInd,1:floor(size(params,2)*0.3))),std(params(paramInd,1:floor(size(params,2)*0.3))),'LineStyle','none','Marker',MarkerType,'Color','k');hold on
                    errorbar(2,mean(params(paramInd,end-floor(size(params,2)*0.3)+1:end)),std(params(paramInd,end-floor(size(params,2)*0.3)+1:end)),'LineStyle','none','Marker',MarkerType,'Color','r');hold on
                    statsTable=[statsTable;mean(params(paramInd,1:floor(size(params,2)*0.3))) std(params(paramInd,1:floor(size(params,2)*0.3))) mean(params(paramInd,end-floor(size(params,2)*0.3)+1:end)) std(params(paramInd,end-floor(size(params,2)*0.3)+1:end))];
                    if paramInd==1
                        title(['Monkey ',num2str(animalInd)]);
                    end
                    xlim([0 3]);
                    set(gca,'XTick',[1 2],'XTickLabel',[{' '} {' '}]);
                    %                 yLimVals=get(gca,'YLim');
                    %                 set(gca,'YLim',[0 yLimVals(2)]);
                    if animalInd+2*(areaInd-1)+4*(paramInd-1)==1
                        xlabel('contrast (%)');
                        set(gca,'XTick',[1 2],'XTickLabel',[{'early'} {'late'}]);
                    end
                    if animalInd==1&&areaInd==1
                        ylabel(paramText{paramInd});
                    end
                end
            end
        end
    end
    statsTable=[statsTable(1:6,:);statsTable(13:18,:);statsTable(7:12,:);statsTable(19:24,:)];%rows 1 to 6: monkey 1. rows 7 to 12, monkey 2. rows 1 to 12, V4. rows 13 to 24, V1. within each group of 6 rows, slope, PNE, min, max, lower threshold, upper threshold
    for sampleContrastsInd=1:length(sampleContrasts)
        for paramInd=1:size(params,1)
            [p,t,stats]=anovan(allParams{sampleContrastsInd}(paramInd,:),{sessArr{sampleContrastsInd}(paramInd,:),animalArr{sampleContrastsInd}(paramInd,:),areaArr{sampleContrastsInd}(paramInd,:)},'model','full');
            p1{paramInd}=p;
            t1{paramInd}=t;
            stats1{paramInd}=stats;
            figure;
            [cRTceANOVA{paramInd},m,h]=multcompare(stats,'dimension',[1 2 3])
            Fs(paramInd,1)=t{2,6};
        end
    end
end
if roving==0
    pcpseslFig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.4, 0.5]);
end
markerTexts='+x';
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        [sampleContrasts testContrasts]=area_metadata(area);
        channels=main_channels(animal,area);
        sessionNums=main_raw_sessions_final(animal,area,[],0);
        if roving==1
            if animalInd==1
                pcpseslFig(areaInd)=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.4, 0.5]);
            else
                figure(pcpseslFig(areaInd))
            end
        end
        for sampleContrastsInd=1:length(sampleContrasts)
            sampleContrast=sampleContrasts(sampleContrastsInd);
            testContrast=testContrasts(sampleContrastsInd,:); 
            colmapText=colormap(jet(size(testContrast,2)));
            colmapText=[colmapText(1,:);132/255 22/255 216/255;202/255 65/255 223/255;colmapText(3:7,:);157/255 212/255 61/255;colmapText(10:12,:);178/255 111/255 12/255;colmapText(end,:)];
            loadText=['load F:\PL\ROC_zero_one\',animal,'\new_vs_old_sglroc3acrosschannels\cumulative_ROCs_old_new_',area,'_',num2str(sampleContrast),'_cutoff10.mat'];
            eval(loadText)
            if roving==1
                figure(pcpseslFig(areaInd))
            elseif roving==0
                figure(pcpseslFig)
            end
            allMeanPerf=all_rocvals;
            if roving==0
                subplot(length(areas),2,animalInd+2*(areaInd-1));
            elseif roving==1
                subplot(length(sampleContrasts),2,animalInd+2*(sampleContrastsInd-1));
            end
            markerText=markerTexts(2);markerS=8;
            for i=1:length(testContrast)
                plot(1:size(allMeanPerf,1),allMeanPerf(:,1+i)*100,'Color',colmapText(i,:),'LineStyle','none','Marker',markerText,'MarkerFaceColor',colmapText(i,:),'MarkerEdgeColor',colmapText(i,:),'MarkerSize',markerS);hold on%'MarkerFaceColor',[1/i 1/i 1/i],
                if animalInd==2&&areaInd==1
                    plotLinear=0;
                    startpoint5=1;
                else
                    plotLinear=0;
                    startpoint5=1;
                end
                [chiselinearTemp(i,:) chiseTemp(i,:) coefperflinearTemp(i,:) coefperfTemp(i,:) aRSlinearTemp(i,:) aRSTemp(i,:)]=bj_linearexpo_fitting(testContrast,allMeanPerf(:,1+i)*100,i,startpoint5,'ROC',plotLinear,[],[]);
            end
            xlim([0 size(allMeanPerf,1)+1]);
            ylim([0 100]);
            for i=1:length(testContrast)
                if strcmp(animal,'jack')
                    yLimVals=get(gca,'ylim');
                    xLimVals=get(gca,'xlim');
                    unitSpace=(yLimVals(2)-yLimVals(1))/30;
                    text('Position',[xLimVals(2)+(xLimVals(2)-xLimVals(1))/25 yLimVals(1)+unitSpace*i*2],'FontSize',9,'String',[markerText,'  ',num2str(testContrast(i)),'%'],'Color',colmapText(i,:));
                end
            end
            if roving==0
                chiselinear{animalInd+2*(areaInd-1)}=chiselinearTemp;
                chise{animalInd+2*(areaInd-1)}=chiseTemp;
                coefperflinear{animalInd+2*(areaInd-1)}=coefperflinearTemp;
                coefperf{animalInd+2*(areaInd-1)}=coefperfTemp;
                aRSlinear{animalInd+2*(areaInd-1)}=aRSlinearTemp;
                aRS{animalInd+2*(areaInd-1)}=aRSTemp;
            elseif roving==1
                chiselinear{animalInd+2*(sampleContrastsInd-1)}=chiselinearTemp;
                chise{animalInd+2*(sampleContrastsInd-1)}=chiseTemp;
                coefperflinear{animalInd+2*(sampleContrastsInd-1)}=coefperflinearTemp;
                coefperf{animalInd+2*(sampleContrastsInd-1)}=coefperfTemp;
                aRSlinear{animalInd+2*(sampleContrastsInd-1)}=aRSlinearTemp;
                aRS{animalInd+2*(sampleContrastsInd-1)}=aRSTemp;
            end
            legend('hide');
            xlabel('');
            ylabel('');
            set(gca,'YTick',0:20:100,'YTickLabel',0:0.2:1);
            if areaInd==1
                title(['Monkey ',num2str(animalInd)]);
                if animalInd==1
                    xlabel('session number');
                    ylabel('PROBMAT');
                end
            end
            if roving==0
                if animalInd==1&&areaInd==1&&sampleContrastsInd==1
                    pc_condcoefFig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
                else
                    figure(pc_condcoefFig);
                end
            elseif roving==1
                if animalInd==1&&sampleContrastsInd==1
                    pc_condcoefFig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.5, 0.8]);
                else
                    figure(pc_condcoefFig);
                end
            end
            if roving==0
                subplot(length(areas),2,animalInd+2*(areaInd-1));
            elseif roving==1
                subplot(length(sampleContrasts),2,animalInd+2*(sampleContrastsInd-1));
            end
            ind=find(testContrast<sampleContrast);
            ind=ind(end);
            [rperfcoeff(animalInd+2*(sampleContrastsInd-1),1:2) pperfcpef(animalInd+2*(sampleContrastsInd-1),1:2)]=bj_plot_coef_expo_perf(coefperf{animalInd+2*(sampleContrastsInd-1)}(:,1),sampleContrast,testContrast,ind);
        end
    end
end
if roving==1&&strcmp(area,'v1_2_2')
    subplot(3,2,1);
    xlabel('difference in contrast (%)');
    ylabel('coefficient a');
    subplot(3,2,2);
    ylim([-50 70]);
    xlim([0 70]);
    subplot(3,2,1);
    xlim([0 70]);
    subplot(3,2,5);
    xlim([0 50]);
    subplot(3,2,6);
    xlim([0 50]);    
end
%compare sizes of adjusted R-square values, between linear and polynomial fits, for PC divided by test contrast
aRSlinexpo1=[];
aRSlinexpo2=[];
for animalInd=1:length(animals)
    for areaInd=1:2%hard-coded for V4_1 and V1_1
        aRSlinexpo1=[aRSlinexpo1;aRSlinear{animalInd+2*(areaInd-1)}(:,:)];
        aRSlinexpo2=[aRSlinexpo2;aRS{animalInd+2*(areaInd-1)}(:,:)];
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
