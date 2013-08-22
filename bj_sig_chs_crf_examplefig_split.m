function bj_sig_chs_crf_examplefig_split
%Written by Xing 20/08/13
%Modified from bj_sig_chs_crf_examplefig, places V4 and V1 channels in separate figures.
%Draw CRFs for channels with significant changes in slope and/or PNE
excludeSessHighSSE=0;
useColMap=1;
analysisType='CRF';
combineSlopeC50=1;
onExternalHD=0;
if onExternalHD==1
    rootFolder='G:\PL_backup_060413';
else
    rootFolder='F:';
end
areaTexts=[{'V4'} {'V1'}];
calculateTangent=1;
plotSlopeFig=1;
if plotSlopeFig==1
    for areaInd=1:2
        figExamples(areaInd)=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.6, 0.8]); %
        set(figExamples(areaInd), 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
        %example figures for channels with significant changes in slope at 30%:
        if excludeSessHighSSE==0
            animals=[{'blanco'} {'jack'} {'blanco'} {'jack'}];
            allChannels=[{[2,24,42,53,59,60]} {[1 2 3 6 8 10 24 35 40 41 52 53]} {[12]} {[9,12,17,18,20,21,22,28,32,55]}];%Spearman's correlation results, Bonferroni corrected
            areas=[{'v4_1'} {'v4_1'} {'v4_1'} {'v1_1'}];
            if combineSlopeC50==1
                if areaInd==1
                    animals=[{'blanco'} {'jack'} {'blanco'} {'jack'}];
                    allChannels=[{[2,24,42,59,60]} {[1 3 6 8 10 24 35 40 41 52 53]} {[12]} {[54]}];%Spearman's correlation results, Bonferroni corrected
                    areas=[{'v4_1'} {'v4_1'} {'v4_1'} {'v4_1'}];
                    sigC50=[1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];%flag channels that have sig C50 changes
                    sigSlope=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];%flag channels that have sig slope changes
                elseif areaInd==2
                    animals={'jack'};
                    allChannels={[9,17,18,20,21,22,28,32,55,7,19,25]};%Spearman's correlation results, Bonferroni corrected
                    areas={'v1_1'};
                    sigC50=[1 1 1 1 1 1 1 1 1 1 1 1];%flag channels that have sig C50 changes
                    sigSlope=[1 1 1 1 1 1 1 1 1 0 1 0];%flag channels that have sig slope changes
                end
            end
        end
        sampleContrast=30;
        allChInd=0;
        for animalInd=1:length(animals)
            animal=animals{animalInd};
            channels=allChannels{animalInd};
            area=areas{animalInd};
            for chInd=1:length(channels)
                allChInd=allChInd+1;
                %         figROCnew=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
                %         set(figROCnew, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                [sampleContrasts testContrasts]=area_metadata(area);
                testContrast=testContrasts;
                matname=['CRF_Ch',num2str(channels(chInd)),'_',num2str(sampleContrast),'_1024_to_1536'];
                pathname=fullfile(rootFolder,'PL',analysisType,animal,area,matname);
                loadText=['load ',pathname,'.mat'];
                eval(loadText)
                %         subplot(ceil(29/5),5,allChInd);
                if combineSlopeC50==1
                    subplot(ceil(length(sigC50)/5),5,allChInd);
                else
                    subplot(ceil(length(sigC50)/5),5,allChInd);
                end
                xvals=testContrast(1):1:testContrast(end);
                %             copperCols=colormap(copper(size(CRFmat,1)));
                %             copperCols=colormap(cool(size(CRFmat,1)+ceil(size(CRFmat,1)*0.25)));
                copperCols1=[];
                copperCols2=[];
                for colMapInd=1:ceil(size(CRFmat,1)/2)
                    copperCols1(colMapInd,:)=[1 0 (colMapInd-1)/ceil((size(CRFmat,1))/2)];
                end
                for colMapInd=1:size(CRFmat,1)-size(CRFmat,1)/2
                    copperCols2(colMapInd,:)=[1-colMapInd/floor((size(CRFmat,1))/2) 0 1];
                end
                copperCols=[copperCols1;copperCols2];
                chSSE=[];
                slopeNeuro=[];
                c50=[];
                diffc50=[];
                minRate=[];
                maxRate=[];
                yLimData=[];
                threshold82lower=[];
                threshold82higher=[];
                for sessionInd=1:size(CRFmat,1)
                    datavals=CRFmat{sessionInd,3};
                    [slopeNeuroNew,c50,diffPNENew,minRateNew,maxRateNew,chSSENew,xvals,yvals]=nr_fitting(datavals,sampleContrast,testContrast,sessionInd,slopeNeuro,chSSE,c50,minRate,maxRate,diffc50,0,calculateTangent,[],animal,area);
                    if useColMap==1
                        plot(xvals,yvals,'Color',copperCols(sessionInd,:));
                    else
                        plot(xvals,yvals,'Color',[1-sessionInd/size(CRFmat,1) 0 sessionInd/size(CRFmat,1)]);
                    end
                    hold on
                    if strcmp(area,'v4_1')
                        xlim([10 60]);
                    elseif strcmp(area,'v1_1')
                        xlim([5 90]);
                    end
                    title(num2str(channels(chInd)));
                    title(allChInd);
                    if allChInd==1
                        xlabel('contrast (%)');
                        ylabel('firing rate (spikes/s)');
                    end
                end
                ylimvals=get(gca,'YLim');
                for sessionInd=1:size(CRFmat,1)
                    if useColMap==1
                        if combineSlopeC50==1&&sigC50(chInd)==1
                            plot([c50(sessionInd) c50(sessionInd)],[0 ylimvals(2)],'Color',copperCols(sessionInd,:));
                        end
                    else
                        if combineSlopeC50==1&&sigC50(chInd)==1
                            plot([c50(sessionInd) c50(sessionInd)],[1-sessionInd/size(CRFmat,1) 0 sessionInd/size(CRFmat,1)]);
                        end
                    end
                end
                if strcmp(area,'v4_1')
                    if sigSlope(chInd)
                        text(12,ylimvals(2)-ylimvals(2)/10,'S','Color',[0.1 0.1 0.1]);
                    end
                    if sigC50(chInd)
                        text(16,ylimvals(2)-ylimvals(2)/10,'C','Color',[0.1 0.1 0.1]);
                    end
                elseif strcmp(area,'v1_1')
                    if sigSlope(chInd)
                        text(7,ylimvals(2)-ylimvals(2)/10,'S','Color',[0.1 0.1 0.1]);
                    end
                    if sigC50(chInd)
                        text(13,ylimvals(2)-ylimvals(2)/10,'C','Color',[0.1 0.1 0.1]);
                    end
                end
                ylim([0 ylimvals(2)]);
            end
        end
        if combineSlopeC50==1
            for newsubplotInd=1:length(sigC50)
                subplot(ceil(length(sigC50)/5),5,newsubplotInd);
                xlimvals=get(gca,'XLim');set(gca,'XTick',[xlimvals(1),xlimvals(2)]);
                ylimvals=get(gca,'YLim');set(gca,'YTick',[ylimvals(1),ylimvals(2)]);
            end
        end
        imagename=['example_sig_chs_change_slope',areaTexts{areaInd}];
        if excludeSessHighSSE==0
            imagename=[imagename,'_allSess'];
        end
        if useColMap==1
            imagename=[imagename,'_copper'];
            imagename=[imagename,'_cool'];
        end
        pathname=fullfile(rootFolder,'PL',analysisType,imagename);
        printtext=sprintf('print -dpng %s.png',pathname);
        set(gcf,'PaperPositionMode','auto')
        eval(printtext);
    end
end
