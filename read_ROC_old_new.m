function read_ROC_old_new(cutoff,excludeSuppressed,normalised,roving)
analysisType='ROC';
analysisType='ROC_zero_one';
newColCode=1;
fitCurves=1;
%Written by Xing 19/05/13
%Compare neurometric function parameters between new and old methods of
%ROC/AUROC value generation
onExternalHD=0;
if onExternalHD==1
    rootFolder='K:\PL_backup_190713';
else
    rootFolder='F:';
end
figGauss=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.45, 0.65]); %
set(figGauss, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 3.305 3.325/0.4*0.45 3.305/5*65]);
animals=[{'blanco'} {'jack'}];
if roving==0
    areas=[{'v4_1'} {'v1_1'}];
    allTableStats=[];
    allParTableStats=[];
    mpallTableStats=[];
elseif roving==1
    areas=[{'v1_2_1'} {'v1_2_2'} {'v1_2_3'}];
    areas=[{'v1_2_1'} {'v1_2_2'}];
    rovingCol=[0.5 0 0;0.75 0 0;1 0 0];
    rovingCol=[1 0 0;10/255 170/255 60/255;0 0 1];
    if newColCode==1
        taskCols=[204/256 153/256 256/256;255/256 222/256 173/256;154/256 205/256 50/256];%purple;orange;green
        taskColsDark=[153/256 50/256 204/256;210/256 105/256 30/256;34/256 139/256 34/256];%purple;orange;green
        taskColsVDark=[84/256 3/256 163/256;166/256 121/256 3/256;79/256 117/256 2/256];%purple;orange;green, for fitted line
    end
    allTableStats=cell(1,3);
    allParTableStats=cell(1,3);
    mpallTableStats=cell(1,3);
    slopeThird=cell(1,length(animals));%compare last third of pre-flanker and first third of flanker sessions
    C50Third=cell(1,length(animals));
    minThird=cell(1,length(animals));
    maxThird=cell(1,length(animals));
    sampleThird=cell(1,length(animals));
    areaThird=cell(1,length(animals));
    area2=cell(1,length(animals));%to compare all flanker sessions with last third of pre-flanker and all of post-flanker
    area1=cell(1,length(animals));
    slope1=cell(1,length(animals));
    C501=cell(1,length(animals));
    min1=cell(1,length(animals));
    max1=cell(1,length(animals));
    slope2=cell(1,length(animals));
    C502=cell(1,length(animals));
    min2=cell(1,length(animals));
    max2=cell(1,length(animals));
    area3=cell(1,length(animals));
    slope3=cell(1,length(animals));
    C503=cell(1,length(animals));
    min3=cell(1,length(animals));
    max3=cell(1,length(animals));    
    allSlope=cell(1,length(animals));
    allC50=cell(1,length(animals));
    allMin=cell(1,length(animals));
    allMax=cell(1,length(animals));   
    allArea=cell(1,length(animals));   
    sample1=cell(1,length(animals));
    sample2=cell(1,length(animals));   
    sample3=cell(1,length(animals));   
    allSample=cell(1,length(animals));
end
for areaInd=1:length(areas)
    area=areas{areaInd};
    [sampleContrasts testContrasts]=area_metadata(area);
    if roving==1&&newColCode==1
        taskCol=taskCols(areaInd,:);
        taskColDark=taskColsDark(areaInd,:);
        markerFaceTypes=[{'none'} {taskCol} {taskColDark}];
    end
    for animalInd=1:length(animals)
        animal=animals{animalInd};
        for sampleContrastsInd=1:length(sampleContrasts)
            sampleContrast=sampleContrasts(sampleContrastsInd);
            if roving==1&&newColCode==1
                if sampleContrast==30
                    markerFaceType=markerFaceTypes{2};
                    fittedLineCol=taskColsDark(areaInd,:);
                elseif sampleContrast==20
                    markerFaceType=markerFaceTypes{1};
                    fittedLineCol=taskCols(areaInd,:);
                elseif sampleContrast==40
                    markerFaceType=markerFaceTypes{3};
                    fittedLineCol=taskColsVDark(areaInd,:);
                end
            end
%         if cutoff==1
%             loadText=['load ',rootFolder,'\PL\',analysisType,'\',animal,'\new_vs_old_sglrocmeanchannels\cumulative_ROCs_old_new_',area,'_30.mat'];
%             loadText=['load ',rootFolder,'\PL\',analysisType,'\',animal,'\new_vs_old_sglroc3acrosschannels\cumulative_ROCs_old_new_',area,'_30.mat'];
%         else
            subFolder='new_vs_old_sglrocmeanchannels';
            subFolder='new_vs_old_sglroc3acrosschannels';
            if excludeSuppressed==1
                subFolder=[subFolder,'_excludeSuppressed'];
            end
            if normalised
                subFolder=[subFolder,'_normalised'];
            end
            loadText=['load ',rootFolder,'\PL\',analysisType,'\',animal,'\',subFolder,'\cumulative_ROCs_old_new_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10),'.mat'];  
            %         end
            eval(loadText)
            newValues=[slopeNeuroNew;PNENew;minRateNew;maxRateNew];
            oldValues=[slopeNeuroOld;PNEOld;minRateOld;maxRateOld];
            dfs=[statsS.df statsP.df statsmin.df statsmax.df];
            tvals=[statsS.tstat statsP.tstat statsmin.tstat statsmax.tstat];
            pvals=[pS pP pmin pmax];
            if roving==1%compare activity levels from last third of flankerless sessions against first third of flanker sessions
                numThirdSess=floor(size(newValues,2)*0.3);
                if strcmp(area,'v1_2_1')%last third of pre-flanker sessions
                    sampleThird{animalInd}=[sampleThird{animalInd} zeros(1,numThirdSess)+sampleContrastsInd];
                    slopeThird{animalInd}=[slopeThird{animalInd} slopeNeuroNew(end-numThirdSess+1:end)];
                    C50Third{animalInd}=[C50Third{animalInd} PNENew(end-numThirdSess+1:end)];
                    minThird{animalInd}=[minThird{animalInd} minRateNew(end-numThirdSess+1:end)];
                    maxThird{animalInd}=[maxThird{animalInd} maxRateNew(end-numThirdSess+1:end)];
                    areaThird{animalInd}=[areaThird{animalInd} zeros(1,numThirdSess)+1];
                    slope1{animalInd}=[slope1{animalInd} slopeNeuroNew(end-numThirdSess+1:end)];
                    C501{animalInd}=[C501{animalInd} PNENew(end-numThirdSess+1:end)];
                    min1{animalInd}=[min1{animalInd} minRateNew(end-numThirdSess+1:end)];
                    max1{animalInd}=[max1{animalInd} maxRateNew(end-numThirdSess+1:end)];
                    area1{animalInd}=[area1{animalInd} zeros(1,numThirdSess)+1];
                    sample1{animalInd}=[sample1{animalInd} zeros(1,numThirdSess)+sampleContrastsInd];
                elseif strcmp(area,'v1_2_2')%first third of flanker sessions
                    sampleThird{animalInd}=[sampleThird{animalInd} zeros(1,numThirdSess)+sampleContrastsInd];
                    slopeThird{animalInd}=[slopeThird{animalInd} slopeNeuroNew(1:numThirdSess)];
                    C50Third{animalInd}=[C50Third{animalInd} PNENew(1:numThirdSess)];
                    minThird{animalInd}=[minThird{animalInd} minRateNew(1:numThirdSess)];
                    maxThird{animalInd}=[maxThird{animalInd} maxRateNew(1:numThirdSess)];
                    areaThird{animalInd}=[areaThird{animalInd} zeros(1,numThirdSess)+2];
                    slope2{animalInd}=[slope2{animalInd} slopeNeuroNew];%store all of flanker data
                    C502{animalInd}=[C502{animalInd} PNENew];
                    min2{animalInd}=[min2{animalInd} minRateNew];
                    max2{animalInd}=[max2{animalInd} maxRateNew];
                    area2{animalInd}=[area2{animalInd} zeros(1,length(slopeNeuroNew))+2];
                    sample2{animalInd}=[sample2{animalInd} zeros(1,length(slopeNeuroNew))+sampleContrastsInd];
                elseif strcmp(area,'v1_2_3')
                    slope3{animalInd}=[slope3{animalInd} slopeNeuroNew];%store all of post-flanker data
                    C503{animalInd}=[C503{animalInd} PNENew];
                    min3{animalInd}=[min3{animalInd} minRateNew];
                    max3{animalInd}=[max3{animalInd} maxRateNew]; 
                    area3{animalInd}=[area3{animalInd} zeros(1,length(slopeNeuroNew))+3];   
                    sample3{animalInd}=[sample3{animalInd} zeros(1,length(slopeNeuroNew))+sampleContrastsInd];
                end                                
            end
            for paramInd=1:4
                if roving==0
                    subplot(4,4,(paramInd-1)*4+animalInd+(areaInd-1)*2);
                    %calculate mean, SD for each parameter, for new & old methods:
                    tableStats(paramInd,1)=mean(oldValues(paramInd,:));
                    tableStats(paramInd,2)=std(oldValues(paramInd,:));
                    tableStats(paramInd,3)=mean(newValues(paramInd,:));
                    tableStats(paramInd,4)=std(newValues(paramInd,:));
                    %ct-test stats for difference between old & new values:
                    tableStats(paramInd,5)=dfs(paramInd);
                    tableStats(paramInd,6)=tvals(paramInd);
                    tableStats(paramInd,7)=pvals(paramInd);
                    %             plot(1:length(slopeNeuroNew),newValues(paramInd,:),'ro','LineStyle','none','MarkerFaceColor','r');hold on
                    %             plot(1:length(slopeNeuroNew),oldValues(paramInd,:),'bo','LineStyle','none','MarkerFaceColor','b');
                    %             plot(1:length(slopeNeuroNew),newValues(paramInd,:),'rx','LineStyle','none');hold on
                    %             plot(1:length(slopeNeuroNew),oldValues(paramInd,:),'bx','LineStyle','none','MarkerSize',7);
                    %             plot(1:length(slopeNeuroNew),oldValues(paramInd,:),'bx','LineStyle','none','MarkerSize',7);hold on
                    %             plot(1:length(slopeNeuroNew),newValues(paramInd,:),'r+','LineStyle','none');
                    plot(1:length(slopeNeuroNew),newValues(paramInd,:),'ro','LineStyle','none');hold on
                    plot(1:length(slopeNeuroNew),oldValues(paramInd,:),'b+','LineStyle','none','MarkerSize',7);
                    startPoint=0;
                elseif roving==1
                    if newColCode==0
                        subplot(4,4,(paramInd-1)*4+animalInd+(areaInd-1)*2);
                    elseif newColCode==1
                        subplot(4,2,(paramInd-1)*2+animalInd);%plot values for stages without and with flankers on same subplot
                        if areaInd==1
                            startPointSave(animalInd)=length(slopeNeuroNew)+1;%keep track of flankerless sessions
                            startPoint=0;
                        elseif areaInd==2
                            startPoint=startPointSave(animalInd);
                            startPointSave2(animalInd)=length(slopeNeuroNew);%keep track of flankerless+flanker sessions
                        elseif areaInd==3
                            startPoint=startPointSave(animalInd)+startPointSave2(animalInd);
                        end
                    end
                    if newColCode==0
                        plot(1:length(slopeNeuroNew),newValues(paramInd,:),'Marker','o','Color',rovingCol(sampleContrastsInd,:),'LineStyle','none');hold on
                        if fitCurves==1
                            bj_linearexpo_fitting([1:length(slopeNeuroNew)]',[newValues(paramInd,:)]',0,0,'ROC');
                        end
                    elseif newColCode==1
                        plot(1+startPoint:length(slopeNeuroNew)+startPoint,newValues(paramInd,:),'Marker','o','MarkerEdgeColor',taskColDark,'LineStyle','none','MarkerFaceColor',markerFaceType);hold on
                        if fitCurves==1
                            bj_linearexpo_fitting([1+startPoint:length(slopeNeuroNew)+startPoint]',[newValues(paramInd,:)]',0,0,'ROC',[],startPoint,fittedLineCol);
                        end
                    end
                end
                %plot(1:length(slopeNeuroNew),newValues(paramInd,:),'ro','LineStyle','none');hold on
                %plot(1:length(slopeNeuroNew),oldValues(paramInd,:),'bx','LineStyle','none','MarkerSize',7);
                if newColCode==0
                    xlim([0 length(slopeNeuroNew)+1]);
                    axis square
                elseif newColCode==1
                    xlim([0 length(slopeNeuroNew)+startPoint+1]);
                end
            end
            %calculate partial correlations for each parameter (based on new data only):
            [rho p]=partialcorr(slopeNeuroNew',(1:length(slopeNeuroNew))',[PNENew' minRateNew' maxRateNew'])
            parTableStats(1,8:10)=[length(slopeNeuroNew)-2 rho p];
            [rho p]=partialcorr(PNENew',(1:length(PNENew))',[slopeNeuroNew' minRateNew' maxRateNew'])
            parTableStats(2,8:10)=[length(PNENew)-2 rho p];
            [rho p]=partialcorr(minRateNew',(1:length(minRateNew))',[slopeNeuroNew' PNENew' maxRateNew'])
            parTableStats(3,8:10)=[length(minRateNew)-2 rho p];
            [rho p]=partialcorr(maxRateNew',(1:length(maxRateNew))',[slopeNeuroNew' PNENew' minRateNew'])
            parTableStats(4,8:10)=[length(maxRateNew)-2 rho p];
            allParTableStats{sampleContrastsInd}=[allParTableStats{sampleContrastsInd};parTableStats];
            %calculate regular correlations
            [rho p]=corr(slopeNeuroNew',(1:length(slopeNeuroNew))','type','Spearman')
            tableStats(1,8:10)=[length(slopeNeuroNew)-2 rho p];
            [rho p]=corr(PNENew',(1:length(PNENew))','type','Spearman')
            tableStats(2,8:10)=[length(PNENew)-2 rho p];
            [rho p]=corr(minRateNew',(1:length(minRateNew))','type','Spearman')
            tableStats(3,8:10)=[length(minRateNew)-2 rho p];
            [rho p]=corr(maxRateNew',(1:length(maxRateNew))','type','Spearman')
            tableStats(4,8:10)=[length(maxRateNew)-2 rho p];
            allTableStats{sampleContrastsInd}=[allTableStats{sampleContrastsInd};tableStats];
            %calculate official version of partial correlations
            [rho p]=partialcorr(slopeNeuroNew',(1:length(slopeNeuroNew))',[minRateNew' maxRateNew'])
            mptableStats(1,8:10)=[length(slopeNeuroNew)-2 rho p];
            [rho p]=partialcorr(PNENew',(1:length(PNENew))',[minRateNew' maxRateNew'])
            mptableStats(2,8:10)=[length(PNENew)-2 rho p];
            [rho p]=partialcorr(minRateNew',(1:length(minRateNew))',[slopeNeuroNew' maxRateNew'])
            mptableStats(3,8:10)=[length(minRateNew)-2 rho p];
            [rho p]=partialcorr(maxRateNew',(1:length(maxRateNew))',[slopeNeuroNew' minRateNew'])
            mptableStats(4,8:10)=[length(maxRateNew)-2 rho p];
            mpallTableStats{sampleContrastsInd}=[mpallTableStats{sampleContrastsInd};mptableStats];
            %calculate most meaningful partial correlations (IMO)
%         [rho p]=partialcorr(slopeNeuroNew',(1:length(slopeNeuroNew))',[minRateNew' maxRateNew'])
%         mptableStats(1,8:10)=[length(slopeNeuroNew)-2 rho p];
%         [rho p]=partialcorr(PNENew',(1:length(PNENew))',[minRateNew' maxRateNew'])
%         mptableStats(2,8:10)=[length(PNENew)-2 rho p];
%         [rho p]=partialcorr(minRateNew',(1:length(minRateNew))',maxRateNew')
%         mptableStats(3,8:10)=[length(minRateNew)-2 rho p];
%         [rho p]=partialcorr(maxRateNew',(1:length(maxRateNew))',minRateNew')
%         mptableStats(4,8:10)=[length(maxRateNew)-2 rho p];
%         mpallTableStats=[mpallTableStats;mptableStats];      
        end
    end
end
% if normalised==0
if roving==0||newColCode==0
    subplot(4,4,1);
    xlabel('Session number');
    ylabel('Slope');
    title('Monkey 1');
    subplot(4,4,2);
    title('Monkey 2');
    subplot(4,4,3);
    title('Monkey 1');
    subplot(4,4,4);
    title('Monkey 2');
    subplot(4,4,5);
    ylabel('PNE');
    subplot(4,4,9);
    ylabel('Minimum');
    subplot(4,4,13);
    ylabel('Maximum');
elseif newColCode==1
    subplot(4,2,1);
    xlabel('Session number');
    ylabel('Slope');
    title('Monkey 1');
    subplot(4,2,2);
    title('Monkey 2');
    subplot(4,2,3);
    ylabel('PNE');
    subplot(4,2,5);
    ylabel('Minimum');
    subplot(4,2,7);
    ylabel('Maximum');
end
if roving==0
    subplot(4,4,1);
    ylim([0.015 0.05]);
    subplot(4,4,2);
    ylim([0 0.035]);
    subplot(4,4,3);
    ylim([0.007 0.025]);
    subplot(4,4,5);
    ylim([30 39]);
    subplot(4,4,7);
    ylim([28 38]);
    subplot(4,4,12);
    ylim([0 0.4]);
    subplot(4,4,13);
    ylim([0.7 1]);
    subplot(4,4,14);
    ylim([0.7 1]);
    subplot(4,4,15);
    ylim([0.7 1]);
    subplot(4,4,16);
    ylim([0.7 1]);
elseif roving==1
    if newColCode==1
        subplot(4,2,1);
        ylim([0 0.03]);
        subplot(4,2,3);
        ylim([15 50]);
        subplot(4,2,4);
        ylim([25 50]);
        subplot(4,2,5);
        ylim([0 0.6]);
        subplot(4,2,6);
        ylim([0 0.4]);
        text(2,0.55,'no flankers')
        plot(2,0.5,'Marker','o','MarkerEdgeColor',taskColsDark(1,:),'MarkerFaceColor',taskColsDark(1,:));
        plot(2,0.45,'Marker','o','MarkerEdgeColor',taskColsDark(1,:),'MarkerFaceColor',taskCols(1,:));
        plot(2,0.4,'Marker','o','MarkerEdgeColor',taskColsDark(1,:),'MarkerFaceColor','none');
        text(3,0.5,'40%')
        text(3,0.45,'30%')
        text(3,0.4,'20%')
        text(30,0.55,'flankers')
        plot(30,0.5,'Marker','o','MarkerEdgeColor',taskColsDark(2,:),'MarkerFaceColor',taskColsDark(2,:));
        plot(30,0.45,'Marker','o','MarkerEdgeColor',taskColsDark(2,:),'MarkerFaceColor',taskCols(2,:));
        plot(30,0.4,'Marker','o','MarkerEdgeColor',taskColsDark(2,:),'MarkerFaceColor','none');
        text(31,0.5,'40%')
        text(31,0.45,'30%')
        text(31,0.4,'20%')
        subplot(4,2,7);
        ylim([0.75 1]);
        subplot(4,2,8);
        ylim([0.75 1]);
        pause%then show just flankerless data
        subplot(4,2,1);
        xlim([0 35]);
        subplot(4,2,3);
        xlim([0 35]);
        subplot(4,2,5);
        xlim([0 35]);
        subplot(4,2,7);
        xlim([0 35]);
        subplot(4,2,2);
        xlim([0 17]);
        subplot(4,2,4);
        xlim([0 17]);
        subplot(4,2,6);
        xlim([0 17]);
        plot(13,0.55,'Marker','o','MarkerEdgeColor',taskColsDark(1,:),'MarkerFaceColor',taskColsDark(1,:));
        plot(13,0.5,'Marker','o','MarkerEdgeColor',taskColsDark(1,:),'MarkerFaceColor',taskCols(1,:));
        plot(13,0.45,'Marker','o','MarkerEdgeColor',taskColsDark(1,:),'MarkerFaceColor','none');
        text(14,0.55,'40%')
        text(14,0.5,'30%')
        text(14,0.45,'20%')
        subplot(4,2,8);
        xlim([0 17]);
    elseif newColCode==0
        subplot(4,4,3);
        ylim([0 0.025]);
        subplot(4,4,4);
        ylim([0 0.07]);
        subplot(4,4,5);
        ylim([15 50]);
        subplot(4,4,6);
        ylim([25 50]);
        subplot(4,4,7);
        ylim([20 45]);
        subplot(4,4,8);
        ylim([25 50]);
        subplot(4,4,9);
        ylim([0 0.6]);
        subplot(4,4,10);
        ylim([0 0.2]);
        subplot(4,4,11);
        ylim([0 0.6]);
        subplot(4,4,12);
        ylim([0 0.4]);
        subplot(4,4,13);
        ylim([0.5 1]);
    end
%     subplot(4,4,14);
%     ylim([0.7 1]);
%     subplot(4,4,15);
%     ylim([0.7 1]);
%     subplot(4,4,16);
%     ylim([0.7 1]);
end
% end
% subplot(4,4,3);
% ylim([0.0025 0.007]);
% subplot(4,4,4);
% ylim([0.011 0.019]);
% subplot(4,4,8);
% ylim([33 40]);
% subplot(4,4,9);
% ylim([-0.4 0.4]);
% subplot(4,4,10);
% ylim([-0.2 0.5]);
% subplot(4,4,11);
% ylim([-0.2 0.4]);
% subplot(4,4,13);
% ylim([0.2 1]);
% subplot(4,4,14);
% ylim([0.2 0.8]);
% subplot(4,4,16);
% ylim([0.78 1]);
% %ignore outliers:
% subplot(4,4,9);
% ylim([0.05 0.4]);
% subplot(4,4,10);
% ylim([0.25 0.45]);
% subplot(4,4,13);
% ylim([0.2 0.6]);
% subplot(4,4,14);
% ylim([0.2 0.5]);
imagename=['bj_old_new_ROC_4_params_sessions_cutoff',num2str(cutoff*10)];
if excludeSuppressed==1
    imagename=['bj_old_new_ROC_4_params_sessions_cutoff',num2str(cutoff*10),'_excludeSuppressed'];
end
if roving==1
    imagename=[imagename,'_roving'];
end
if newColCode==1
    imagename=[imagename,'_newCol'];
end
pathname=fullfile(rootFolder,'PL','ROC',imagename);
printtext=sprintf('print -dpng %s.png',pathname);
set(gcf,'PaperPositionMode','auto')
eval(printtext);
cutoff
allTableStats

allStatsTable=[];
figure;
for animalInd=1:2    
    [p,t,stats]=anovan(slopeThird{animalInd},{areaThird{animalInd},sampleThird{animalInd},},'model','full');
    [c,m] = multcompare(stats)
    [p,t,stats]=anovan(C50Third{animalInd},{areaThird{animalInd},sampleThird{animalInd},},'model','full');
    [c,m] = multcompare(stats)
    [p,t,stats]=anovan(minThird{animalInd},{areaThird{animalInd},sampleThird{animalInd},},'model','full');
    [c,m] = multcompare(stats)
    [p,t,stats]=anovan(maxThird{animalInd},{areaThird{animalInd},sampleThird{animalInd},},'model','full');
    [c,m] = multcompare(stats)
    %to compare all flanker sessions with last third of pre-flanker and all
    %of post-flanker:
    allSlope{animalInd}=[slope1{animalInd} slope2{animalInd} slope3{animalInd}];
    allC50{animalInd}=[C501{animalInd} C502{animalInd} C503{animalInd}];
    allMin{animalInd}=[min1{animalInd} min2{animalInd} min3{animalInd}];
    allMax{animalInd}=[max1{animalInd} max2{animalInd} max3{animalInd}];
    allArea{animalInd}=[area1{animalInd} area2{animalInd} area3{animalInd}];
    allSample{animalInd}=[sample1{animalInd} sample2{animalInd} sample3{animalInd}];
    [p,t,stats]=anovan(allSlope{animalInd},{allArea{animalInd},allSample{animalInd},},'model','full');
    [c,m] = multcompare(stats)
    allStatsTable=[allStatsTable;t(2,3) t(5,3) t(2,6) p(1)];
    [p,t,stats]=anovan(allC50{animalInd},{allArea{animalInd},allSample{animalInd},},'model','full');
    [c,m] = multcompare(stats)
    allStatsTable=[allStatsTable;t(2,3) t(5,3) t(2,6) p(1)];
    [p,t,stats]=anovan(allMin{animalInd},{allArea{animalInd},allSample{animalInd},},'model','full');
    [c,m] = multcompare(stats)
    allStatsTable=[allStatsTable;t(2,3) t(5,3) t(2,6) p(1)];
    [p,t,stats]=anovan(allMax{animalInd},{allArea{animalInd},allSample{animalInd},},'model','full');
    [c,m] = multcompare(stats)
    allStatsTable=[allStatsTable;t(2,3) t(5,3) t(2,6) p(1)];
end
allStatsTable=[allStatsTable(1:4,:) allStatsTable(5:8,:)];%monkey 1 in first 4 columns, monkey 2 in next four oclumns; 1st row: slope, 2nd: C50, 3rd: min, 4th: max
