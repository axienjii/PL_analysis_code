function jack_psycho_batch_V1_roving_paper
%Similar to blanco_SE_roc_batch, contains list of sessions and data
%directories, however, as this analyses psychometric performance, does not
%take spike activity into account and thus the function
%'blanco_2target_psycho_EV_v2' is only run once per session.
%(blanco_SE_roc_batch runs the function blanco_SE_roc 64 times per session- once for each channel).
%If intention is to create data file containing list of Weibull constants
%from each session, must first activate lines of code in
%blanco_2target_psycho_EV_v2 so that values are written to file.

% allMeanPerf=[];save allMeanPerf_20.mat allMeanPerf
% allMeanPerf=[];save allMeanPerf_30.mat allMeanPerf
% allMeanPerf=[];save allMeanPerf_40.mat allMeanPerf

%set accordingly:
animals=[{'blanco'} {'jack'}];
area='v1_4';
psychoOnly=1;
plotLinear=0;
roving=1;
readData=0;
onExternalHD=0;
if onExternalHD==1
    rootFolder='G:\PL_backup_060413';
else
    rootFolder='F:';
end
for animalInd=1:2
    animal=animals{animalInd};
    
    taskCols=[204/256 153/256 256/256;255/256 222/256 173/256;154/256 205/256 50/256];%purple;orange;green
    taskColsDark=[153/256 50/256 204/256;210/256 105/256 30/256;34/256 139/256 34/256];%purple;orange;green
    taskColsVDark=[84/256 3/256 163/256;166/256 121/256 3/256;79/256 117/256 2/256];%purple;orange;green, for fitted line
    
    %RGB version:
    taskCols=[1 0.7 0.7;154/256 205/256 50/256;0.7 0.7 1];%red;green;blue
    taskColsDark=[0.8 0 0;34/256 139/256 34/256;0 0 0.7];
    taskColsVDark=[0.5 0 0;79/256 117/256 2/256;0 0 0.4];
    
    plotRoving=1;
    preflankerOnly=0;%1: flankers absent, 2: flankers present 3: control task, flankers absent 4: control task, flankers present 0: all roving 5: pre-flanker, flanker & post-flanker V1_2
    if plotRoving==1
        if strcmp(animal,'jack')
            if preflankerOnly==1||preflankerOnly==2%show only part of the data
                areas={'v1_2'};
                areaText=[{'V1 RFs location (-0.7,-1.3)'}];
            elseif preflankerOnly==3||preflankerOnly==4%show only part of the control data
                areas={'v1_4'};
                areaText=[{'control location (-3.5,-3)'}];
            elseif preflankerOnly==0
                areas=[{'v1_2'} {'v1_4'}];
                areaText=[{'V1 RFs location (-0.7,-1.3)'} {'control location (-3.5,-3)'}];
            elseif preflankerOnly==5
                areas=[{'v1_2'}];
                areaText=[{'V1 RFs location (-0.7,-1.3)'} {'control location (-3.5,-3)'}];
            end
        elseif strcmp(animal,'blanco')
            areas={'v1_2'};
            areaText={'V1 RFs location (-3.5,-3)'};
        end
        if animalInd==1
            fig=figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.2, 0.04, 0.6, 0.9]);
            set(fig,'PaperUnits','centimeters','PaperType','A4','PaperOrientation', 'portrait', 'PaperPosition', [0.63452 0.63452 31.5/(3-length(areas)) 28.41]);
        end
        for areaInd=1:length(areas)
            area=areas{areaInd};
            [flankerType flankerSessions]=getFlankerType(animal,area);
            sampleContrasts=20:10:40;
            count=0;
            if preflankerOnly==1%show only the pre-flanker data
                flankerType=flankerType(1);
                flankerSessions=flankerSessions(1);
            elseif preflankerOnly==2%show only the pre-flanker and flanker data
                flankerType=flankerType(1:2);
                flankerSessions=flankerSessions(1:2);
            elseif preflankerOnly==3%show only the control pre-flanker data
                flankerType=flankerType(1);
                flankerSessions=flankerSessions(1);
            elseif preflankerOnly==4%show only the control pre-flanker and flanker data
                flankerType=flankerType(1:2);
                flankerSessions=flankerSessions(1:2);
            end
            for flankerInd=1:length(flankerType)
                taskCol=taskCols(flankerInd,:);
                taskColDark=taskColsDark(flankerInd,:);
                markerFaceTypes=[{'none'} {taskCol} {taskColDark}];
                for i=1:length(sampleContrasts)
                    sampleContrast=sampleContrasts(i);
                    if sampleContrast==30
                        markerFaceType=markerFaceTypes{2};
                        fittedLineCol=taskColsDark(flankerInd,:);
                    elseif sampleContrast==20
                        markerFaceType=markerFaceTypes{1};
                        fittedLineCol=taskCols(flankerInd,:);
                    elseif sampleContrast==40
                        markerFaceType=markerFaceTypes{3};
                        fittedLineCol=taskColsVDark(flankerInd,:);
                    end
                    loadText=['load ',rootFolder,'\PL\psycho_data\',animal,'\allMeanPerf\allMeanPerf_',area,'_',num2str(flankerInd),'_',num2str(sampleContrast),'.mat allMeanPerf'];
                    eval(loadText)
                    perf=allMeanPerf(:,2);
                    PSE=allMeanPerf(:,15);
                    slope=allMeanPerf(:,16);
                    subplot(3,2,animalInd);
                    plot(count+1:count+length(flankerSessions{flankerInd}),perf,'LineStyle','none','Marker','o','MarkerEdgeColor',taskColDark,'MarkerFaceColor',markerFaceType,'MarkerSize',5);hold on
                    bj_linearexpo_fitting([count+1:count+length(flankerSessions{flankerInd})]',perf,0,0,'ROC',plotLinear,count,fittedLineCol);
                    if length(areas)==2
                        ylim([0.5 0.95]);
                    else
                        ylim([0.6 0.95]);
                    end
                    title(areaText{areaInd});
                    if areaInd==1
                        ylabel('proportion correct');
                    end
                    set(gca, 'box', 'off');
                    subplot(3,2,animalInd+4);
                    plot(count+1:count+length(flankerSessions{flankerInd}),PSE,'LineStyle','none','Marker','o','MarkerEdgeColor',taskColDark,'MarkerFaceColor',markerFaceType,'MarkerSize',5);hold on
                    bj_linearexpo_fitting([count+1:count+length(flankerSessions{flankerInd})]',PSE,0,0,'ROC',plotLinear,count,fittedLineCol);
                    if length(areas)==2
                        ylim([20 50]);
                    else
                        ylim([15 50]);
                    end
                    if areaInd==1
                        ylabel('PSE');
                    end
                    set(gca, 'box', 'off');
                    subplot(3,2,animalInd+2);
                    plot(count+1:count+length(flankerSessions{flankerInd}),slope,'LineStyle','none','Marker','o','MarkerEdgeColor',taskColDark,'MarkerFaceColor',markerFaceType,'MarkerSize',5);hold on
                    bj_linearexpo_fitting([count+1:count+length(flankerSessions{flankerInd})]',slope,0,0,'ROC',plotLinear,count,fittedLineCol);
                    if length(areas)==2
                        ylim([1 9]);
                    else
                        ylim([1 11]);
                    end
                    if areaInd==1
                        ylabel('slope');
                    end
                    set(gca, 'box', 'off');
                end
                count=count+length(flankerSessions{flankerInd});
            end
        end
        if strcmp(animal,'blanco')
            if preflankerOnly==1
                subplot(3,1,1);
                ylim([0.7 0.9]);
                xlim([0 35]);
                subplot(3,1,2);
                ylim([0 8]);
                xlim([0 35]);
                subplot(3,1,3);
                ylim([25 40]);
                xlim([0 35]);
            end
            if preflankerOnly==0
                subplot(3,1,1);
                xlim([0 57]);
                subplot(3,1,2);
                xlim([0 57]);
                subplot(3,1,3);
                xlim([0 57]);
            end
            if preflankerOnly==5
                subplot(3,2,1);
                xlabel('session number');
                xlim([0 57]);
                ylim([0.7 0.95]);
                subplot(3,2,3);
                xlim([0 57]);
                ylim([0 11]);
                subplot(3,2,5);
                ylim([15 55]);
                xlim([0 57]);
            end
        end
        if strcmp(animal,'jack')
            if strcmp(area,'v1_4')&&preflankerOnly~=0
                subplot(3,1,1);
                xlabel('session number');
                subplot(3,1,2);
                ylim([0 9]);
                subplot(3,1,3);
                ylim([20 40]);
            elseif strcmp(area,'v1_2')&&preflankerOnly==2
                subplot(3,1,1);
                xlabel('session number');
                ylim([0.5 0.9]);
                subplot(3,1,2);
                ylim([0 11]);
                subplot(3,1,3);
                ylim([20 45]);
            elseif strcmp(area,'v1_4')&&preflankerOnly==0
                subplot(3,2,1);
                xlabel('session number');
                xlim([0 43]);
                subplot(3,2,2);
                xlim([0 51]);
                ylim([0.6 0.9]);
                subplot(3,2,3);
                ylim([0 9]);
                xlim([0 43]);
                subplot(3,2,4);
                ylim([0 9]);
                xlim([0 51]);
                subplot(3,2,5);
                ylim([20 45]);
                xlim([0 43]);
                subplot(3,2,6);
                ylim([20 40]);
                xlim([0 51]);
            elseif strcmp(area,'v1_4')&&preflankerOnly==4
                subplot(3,1,1);
                xlabel('session number');
                xlim([0 43]);
                subplot(3,1,2);
                xlim([0 43]);
                subplot(3,1,3);
                ylim([20 40]);
                xlim([0 43]);
            elseif strcmp(area,'v1_2')&&preflankerOnly==5
                subplot(3,2,2);
                xlabel('session number');
                xlim([0 43]);
                ylim([0.5 0.95]);
                subplot(3,2,4);
                xlim([0 43]);
                ylim([0 11]);
                subplot(3,2,6);
                ylim([20 45]);
                xlim([0 43]);
            end
            if preflankerOnly==1
                subplot(3,1,1);
                ylim([0.7 0.9]);
                subplot(3,1,2);
                ylim([0 11]);
                subplot(3,1,3);
                ylim([20 40]);
            end
        end
        plotLegend=0;
        if plotLegend==1
            subplot(3,length(areas),length(areas));
            ylims=get(gca,'YLim');
            xlims=get(gca,'XLim');
            sampleTexts=[{'20% sample'} {'30% sample'} {'40% sample'}];
            flankerTexts=[{'no flankers'} {'flankers'} {'no flankers'}];
            markerFaceTypes=[{'none'} {[180/256 180/256 180/256]} {'k'}];
            markerFaceTypesPreflankerOnly=[{'none'} {[204/256 153/256 256/256]} {[153/256 50/256 204/256]}];%purple
            markerFaceTypesFlankerOnly=[{'none'} {[255/256 222/256 173/256]} {[210/256 105/256 30/256]}];%green
            taskColsMedium=[147/256 112/256 219/256;255/256 165/256 0/256;154/256 205/256 50/256];%purple;orange;green
            for i=1:3
                if preflankerOnly==0
                    if length(areas)==2
                        keyText=[sampleTexts{i},'         ',flankerTexts{i}];
                    else
                        keyText=[sampleTexts{i},'            ',flankerTexts{i}];
                    end
                    plot(xlims(1)+(xlims(2)-xlims(1))/20,ylims(1)+(ylims(2)-ylims(1))/12*(i),'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor',markerFaceTypes{i});
                    plot(xlims(1)+(xlims(2)-xlims(1))/3.6,ylims(1)+(ylims(2)-ylims(1))/12*(i),'LineStyle','none','Marker','o','MarkerEdgeColor',taskColsMedium(i,:),'MarkerFaceColor',taskColsMedium(i,:));
                    text('Position',[xlims(1)+(xlims(2)-xlims(1))/15,ylims(1)+(ylims(2)-ylims(1))/12*(i)],'LineStyle','none','FontSize',9,'String',keyText,'Color','k');
                elseif preflankerOnly==1||preflankerOnly==2||preflankerOnly==3||preflankerOnly==4
                    if length(areas)==2
                        keyText=[sampleTexts{i}];
                    else
                        keyText=[sampleTexts{i}];
                    end
                    plot(xlims(1)+(xlims(2)-xlims(1))/20,ylims(1)+(ylims(2)-ylims(1))/12*(i),'LineStyle','none','Marker','o','MarkerEdgeColor',taskColsDark(1,:),'MarkerFaceColor',markerFaceTypesPreflankerOnly{i});
                    text('Position',[xlims(1)+(xlims(2)-xlims(1))/15,ylims(1)+(ylims(2)-ylims(1))/12*(i)],'LineStyle','none','FontSize',9,'String',keyText,'Color','k');
                    if preflankerOnly==2||preflankerOnly==4
                        plot(xlims(1)+xlims(2)-(xlims(2)-xlims(1))/4,ylims(1)+(ylims(2)-ylims(1))/12*(i),'LineStyle','none','Marker','o','MarkerEdgeColor',taskColsDark(2,:),'MarkerFaceColor',markerFaceTypesFlankerOnly{i});
                        text('Position',[xlims(1)+xlims(2)-(xlims(2)-xlims(1))/5,ylims(1)+(ylims(2)-ylims(1))/12*(i)],'LineStyle','none','FontSize',9,'String',keyText,'Color','k');
                    end
                end
            end
        end
    end
end
formats=[{'eps'} {'png'}];
formatsCommand=[{'epsc'} {'png'}];
for j=1:2
    if length(areas)==2
        printFileName=[rootFolder,'\PL\psycho_data\',animal,'\PC_PSE_slope_',areas{1},'_',areas{2},'.',formats{j}];
    elseif length(areas)==1
        if preflankerOnly==0
            printFileName=[rootFolder,'\PL\psycho_data\',animal,'\PC_PSE_slope_',areas{1},'.',formats{j}];
        elseif preflankerOnly==1
            printFileName=[rootFolder,'\PL\psycho_data\',animal,'\PC_PSE_slope_',areas{1},'preflankers.',formats{j}];
        elseif preflankerOnly==2
            printFileName=[rootFolder,'\PL\psycho_data\',animal,'\PC_PSE_slope_v1_2_1_v1_2_2.',formats{j}];
        elseif preflankerOnly==3
            printFileName=[rootFolder,'\PL\psycho_data\',animal,'\PC_PSE_slope_',areas{1},'preflankers.',formats{j}];
        elseif preflankerOnly==4
            printFileName=[rootFolder,'\PL\psycho_data\',animal,'\PC_PSE_slope_v1_4_1_v1_4_2.',formats{j}];
        end
        print(sprintf('-d%s',formatsCommand{j}),'-r300',printFileName)
    end
end
close all

%analyse performance between pairs of sample contrasts
plotPairSample=0;
preflankerOnly=1;%1: flankers absent, 2: flankers present 
if plotPairSample==1
    animals=[{'blanco'} {'jack'}];
    animalText=[{'Monkey 1'} {'Monkey 2'}];
    if preflankerOnly==1
        fig=figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.04, 0.04, 0.4, 0.35]);
    else
        fig=figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.04, 0.04, 0.4 0.35]);
    end
    set(fig,'PaperUnits','centimeters','PaperType','A4','PaperOrientation', 'portrait');
    for animalInd=1:length(animals)
        animal=animals{animalInd};
        markerCols='kcm';%20 vs 30; 30 vs 40; 20 vs 40
        comparisons=[1 2;2 3;1 3];
        [flankerType flankerSessions]=getFlankerType(animal,area);
        sampleContrasts=20:10:40;
        if preflankerOnly==1%show only the pre-flanker data
            flankerType=flankerType(1);
            flankerSessions=flankerSessions(1);
        elseif preflankerOnly==2%show only the pre-flanker data
            flankerType=flankerType(2);
            flankerSessions=flankerSessions(2);
        end
%         markerFaceTypes=[{'none'} {taskCol} {taskColDark}];
        perf=[];PSE=[];slope=[];
        for i=1:length(sampleContrasts)
            sampleContrast=sampleContrasts(i);
            loadText=['load ',rootFolder,'\PL\psycho_data\',animal,'\allMeanPerf\allMeanPerf_v1_2_',num2str(flankerType),'_',num2str(sampleContrast),'.mat allMeanPerf'];
            eval(loadText)
            perf(i,:)=allMeanPerf(:,2);
            PSE(i,:)=allMeanPerf(:,15);
            slope(i,:)=allMeanPerf(:,16);
        end
        for comparisonInd=1:3
            subplot(1,2,animalInd);
%             if sampleContrast==30
%                 markerFaceType=markerFaceTypes{2};
%             elseif sampleContrast==20
%                 markerFaceType=markerFaceTypes{1};
%             elseif sampleContrast==40
%                 markerFaceType=markerFaceTypes{3};
%             end
            plot(perf(comparisons(comparisonInd,2),:),perf(comparisons(comparisonInd,1),:),'LineStyle','none','Marker','o','MarkerEdgeColor',markerCols(comparisonInd),'MarkerFaceColor',markerCols(comparisonInd));hold on
            hold on
            [rPairSamp(animalInd,comparisonInd),pPairSamp(animalInd,comparisonInd)]=corr(perf(comparisons(comparisonInd,1),:)',perf(comparisons(comparisonInd,2),:)');
            dfs(animalInd,comparisonInd)=length(perf(comparisons(comparisonInd,1),:))-2;
            title(animalText{animalInd});
            if animalInd==1
                ylabel('performance');
                xlabel('performance');
            end
        end
        if preflankerOnly==1
            subplot(1,2,1)
            xlim([0.75 0.9]);
            ylim([0.7 0.9]);
            subplot(1,2,2)
            ylim([0.7 0.9]);
        elseif preflankerOnly==2
            subplot(1,2,1)
            xlim([0.7 0.95]);
            ylim([0.6 0.92]);
            subplot(1,2,2)
            ylim([0.5 0.85]);
        end
        subplot(1,2,1);
        ylims=get(gca,'YLim');
        xlims=get(gca,'XLim');
        sampleTexts=[{'20% vs 30%'} {'30% vs 40%'} {'20% vs 40%'}];
        flankerTexts=[{'no flankers'} {'flankers'} {'no flankers'}];
        markerFaceTypes=[{'none'} {[180/256 180/256 180/256]} {'k'}];
        markerFaceTypesPreflankerOnly=[{'none'} {[204/256 153/256 256/256]} {[153/256 50/256 204/256]}];
        taskColsMedium=[147/256 112/256 219/256;255/256 165/256 0/256;154/256 205/256 50/256];%purple;orange;green
        for i=1:3
            if preflankerOnly==0
                if length(areas)==2
                    keyText=[sampleTexts{i}];
                else
                    keyText=[sampleTexts{i}];
                end
                plot(xlims(1)+(xlims(2)-xlims(1))/20,ylims(1)+(ylims(2)-ylims(1))/10*(4-i),'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor',markerFaceTypes{i});
                plot(xlims(1)+(xlims(2)-xlims(1))/3.6,ylims(1)+(ylims(2)-ylims(1))/10*(4-i),'LineStyle','none','Marker','o','MarkerEdgeColor',taskColsMedium(i,:),'MarkerFaceColor',taskColsMedium(i,:));
                text('Position',[xlims(1)+(xlims(2)-xlims(1))/15,ylims(1)+(ylims(2)-ylims(1))/10*(4-i)],'LineStyle','none','FontSize',9,'String',keyText,'Color','k');
            elseif preflankerOnly==1||preflankerOnly==2
                keyText=[sampleTexts{i}];
                plot(xlims(1)+(xlims(2)-xlims(1))/20,ylims(1)+(ylims(2)-ylims(1))/15*(4-i),'LineStyle','none','Marker','o','MarkerEdgeColor',markerCols(i),'MarkerFaceColor',markerCols(i));
                text('Position',[xlims(1)+(xlims(2)-xlims(1))/15,ylims(1)+(ylims(2)-ylims(1))/15*(4-i)],'LineStyle','none','FontSize',9,'String',keyText,'Color','k');
            end
        end
    end
    formats=[{'eps'} {'png'}];
    formatsCommand=[{'epsc'} {'png'}];
    for j=1:2
        if preflankerOnly==0
            printFileName=[rootFolder,'\PL\psycho_data\roving_sample_PC_',areas{1},'.',formats{j}];
        elseif preflankerOnly==1
            printFileName=[rootFolder,'\PL\psycho_data\roving_sample_PC_v1_2_1.',formats{j}];
        elseif preflankerOnly==2
            printFileName=[rootFolder,'\PL\psycho_data\roving_sample_PC_v1_2_2.',formats{j}];
        end
        set(gcf, 'PaperPositionMode', 'auto');
        print(sprintf('-d%s',formatsCommand{j}),'-r300',printFileName)
    end
end
close all