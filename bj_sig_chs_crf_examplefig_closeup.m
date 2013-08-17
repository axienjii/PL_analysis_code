function bj_sig_chs_crf_examplefig_closeup
%Written by Xing 10/8/13
%Modified from bj_sig_chs_crf_examplefig to
%draw CRFs in detail for four example channels with significant changes in slope and/or PNE
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
calculateTangent=1;
plotSlopeFig=1;
if plotSlopeFig==1
    %example figures for channels with significant changes in slope at 30%:
    if excludeSessHighSSE==0
        if combineSlopeC50==1
            animals=[{'jack'} {'blanco'} {'jack'} {'jack'}];
            allChannels=[{[53]} {[42]} {[10]} {[20]}];%Spearman's correlation results, Bonferroni corrected
            areas=[{'v4_1'} {'v4_1'} {'v4_1'} {'v1_1'}];
        end
    end
    figSlope=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.3, 0.8]); %
    set(figSlope, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
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
            subplot(4,2,allChInd*2-1);
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
                [slopeNeuro,c50,diffPNENew,minRateNew,maxRateNew,chSSENew,xvals,yvals]=nr_fitting(datavals,sampleContrast,testContrast,sessionInd,slopeNeuro,chSSE,c50,minRate,maxRate,diffc50,0,calculateTangent,[],animal,area);
                if useColMap==1
                    plot(xvals,yvals,'Color',copperCols(sessionInd,:));
                else
                    plot(xvals,yvals,'Color',[1-sessionInd/size(CRFmat,1) 0 sessionInd/size(CRFmat,1)]);
                end
                if allChInd<=2
%                     text(xvals(end)+3,yvals(end)*1.3,num2str(sessionInd),'Color',copperCols(sessionInd,:));
                end
                hold on
                if strcmp(area,'v4_1')
                    xlim([10 60+5]);
                elseif strcmp(area,'v1_1')
                    xlim([5 90+10]);
                end
                title(num2str(channels(chInd)));
                title(allChInd);
                if allChInd==1
                    xlabel('contrast (%)');
                    ylabel('firing rate (spikes/s)');
                end
            end
            ylimvals=get(gca,'YLim');
            ylim([0 ylimvals(2)]);  
            for sessionInd=1:size(CRFmat,1)
                if allChInd>=3
                    subplot(4,2,allChInd*2-1);
                    plot([c50(sessionInd) c50(sessionInd)],[0 ylimvals(2)],'Color',copperCols(sessionInd,:));
                    %text(c50(sessionInd),ylimvals(2)+3,num2str(sessionInd),'Color',copperCols(sessionInd,:));
                    subplot(4,2,allChInd*2);
                    plot(sessionInd,c50(sessionInd),'Color',copperCols(sessionInd,:),'MarkerFaceColor',copperCols(sessionInd,:),'Marker','o');hold on
                else
                    subplot(4,2,allChInd*2);
                    plot(sessionInd,slopeNeuro(sessionInd),'Color',copperCols(sessionInd,:),'MarkerFaceColor',copperCols(sessionInd,:),'Marker','o');hold on
                end
            end
            if allChInd<=2
                xlabel('session number');
                ylabel('slope');
            else
                xlabel('session number');
                ylabel('C_5_0');
            end     
            xlim([0 size(CRFmat,1)+1]);
        end
    end
    imagename='example_sig_chs_change_slope_detailed';
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
if combineSlopeC50==1
    for newsubplotInd=1:length(sigC50)
        subplot(ceil(length(sigC50)/5),5,newsubplotInd);
        xlimvals=get(gca,'XLim');set(gca,'XTick',[xlimvals(1),xlimvals(2)]);
        ylimvals=get(gca,'YLim');set(gca,'YTick',[ylimvals(1),ylimvals(2)]);
    end
end

