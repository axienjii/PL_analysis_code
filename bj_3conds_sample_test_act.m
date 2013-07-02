function bj_3conds_sample_test_act(animals,readData)
%written by Xing 12/04/13
%Calculates correlations in activity to sample and test

saveSampleTestAct=1;
comparisonType=1;
analysisType='adaptation';
onExternalHD=0;
if onExternalHD==1
    rootFolder='G:\PL_backup_060413';
else
    rootFolder='F:';
end
plotDual=0;
plotFigs=0;
animals=[{'blanco'} {'jack'}];
areas=[{'v4_1'} {'v1_1'} {'v1_2_1'} {'v1_2_2'} {'v1_2_3'}];
test_epochs={0 512 512*2 512*3};durSpon=150;
durSpon=150;%length of period prior to sample onset from which spontaneous rates are calculated. Can take on a value of up to 512 ms.
minTrials=10;%set value of minumum number of trials for inclusion of session
subPeriod=1;
% readData=0;
if readData==1
    for animalInd=1:length(animals)
        animal=animals{animalInd};
        for areaInd=1:length(areas)
            area=areas{areaInd};
            channels=main_channels(animal,area);
            sessionNums = main_raw_sessions_final(animal,area,[],0);
            cellEpochTimes={0 512 512*2 512*3};%{[0 40 300] 529 [529*2 529*2+40 529*2+300] 529*3}
            [sampleContrasts testContrasts]=area_metadata(area);
            if strcmp(area,'v4_1')
                conds=[8 9 10];%test contrast of 31, 32, 33%
            elseif strcmp(area,'v1_1')
                conds=8;%test contrast of 32%
            elseif strcmp(area,'v1_2_1')||strcmp(area,'v1_2_2')||strcmp(area,'v1_2_3')
                conds=[4;9;8];%test contrast of 22, 35, 42% for respective sample contrast conditions
            end
            for sampleContrastsInd=1:length(sampleContrasts)
                sampleContrast=sampleContrasts(sampleContrastsInd);
                testContrast=testContrasts(sampleContrastsInd,:);
                %             colmapText=colormap(jet(size(testContrast,2)));
                %             colmapText=[colmapText(1,:);132/255 22/255 216/255;202/255 65/255 223/255;colmapText(3:7,:);157/255 212/255 61/255;colmapText(10:12,:);178/255 111/255 12/255;colmapText(end,:)];
                colmapText=[1 0 0;0.05 0.73 0.24;0 0 1];
                hs=[];
                ps=[];
                if plotFigs==1
                    figSess=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
                    set(figSess, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                end
                if length(conds)==1
                    allAUROC={[]};
                    cis={[]};
                    stats={[]};
                elseif length(conds)==3
                    allAUROC=[{[]};{[]};{[]}];
                    cis=[{[]};{[]};{[]}];
                    stats=[{[]};{[]};{[]}];
                end
                for h=1:length(channels)
                    channel=channels(h);
                    AUROC=[];
                    if plotFigs==1
                        figCh=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
                        set(figCh, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                    end
                    for i=1:length(sessionNums)
                        if comparisonType==1
                            matFolder=['F:\PL\spikeData\',animal];
                        elseif comparisonType==2
                            matFolder=['F:\PL\spikeData\',animal,'\plotRedArtifacts\correct_trials_only'];
                        end
                        chStr=[num2str(channel),'_',num2str(sessionNums(i)),'_',num2str(sampleContrast),'.mat'];
                        matPath=fullfile(matFolder,chStr);
                        matExists=0;
                        if exist(matPath,'file')
                            matExists=1;
                        end
                        if matExists==1
                            valsText=['load ',matPath,' matarray'];
                            eval(valsText);
                            for subplotMultiple=1:length(conds(sampleContrastsInd,:))%1:3
                                cond=conds(sampleContrastsInd,subplotMultiple);
                                trials2=length(matarray{cond,2});
                                trials4=length(matarray{cond,4});
                                for epoch=[2 4]
                                    actList=[];
                                    for n=1:min([trials2 trials4])
                                        temp=length(matarray{cond,epoch}{n});
                                        actList(n)=temp*1000/512;
                                    end
                                    if epoch==2
                                        epoch2{cond,subPeriod}=actList;%store activity levels from each subperiod for ROC calculation
                                    elseif epoch==4
                                        epoch4{cond,subPeriod}=actList;%store activity levels from each subperiod for ROC calculation
                                    end
                                end
                                if length(epoch2{cond,subPeriod})~=length(epoch4{cond,subPeriod})
                                    pauseHere=1;
                                end
                                higherTest=0;
                                lowerTest=0;
                                for rowInd=1:length(epoch4{cond,subPeriod})
                                    if epoch4{cond,subPeriod}(1,rowInd)>epoch2{cond,subPeriod}(1,rowInd)
                                        higherTest=higherTest+1;
                                    elseif epoch4{cond,subPeriod}(1,rowInd)<epoch2{cond,subPeriod}(1,rowInd)
                                        lowerTest=lowerTest+1;
                                    end
                                end
                                AUROC(cond,i)=higherTest/(higherTest+lowerTest);
                                subplot(ceil(length(sessionNums)/5),5*length(conds(sampleContrastsInd,:)),(i-1)*length(conds(sampleContrastsInd,:))+subplotMultiple);
                                plot(epoch2{cond,subPeriod},epoch4{cond,subPeriod},'x','Color',colmapText(subplotMultiple,:));hold on%test act against sample act
                                if (i-1)*3+subplotMultiple==1
                                    xlabel('sample');
                                    ylabel('test');
                                end
                                xlims=get(gca,'XLim');
                                ylims=get(gca,'YLim');
                                newLims=[min([xlims ylims]) max([xlims ylims])];
                                xlim([newLims(1) newLims(2)]);
                                ylim([newLims(1) newLims(2)]);
                                plot([newLims(1) newLims(2)],[newLims(1) newLims(2)],'k--');
                                axis square
                                if subplotMultiple==2
                                    title(sessionNums(i));
                                end
                            end
                        end
                    end
                    if saveSampleTestAct==1
                        if plotDual==1
                            imagename=['ch_',num2str(channel),'_sample_vs_test_act_',area,'_',num2str(sampleContrast),'_splitPlots'];
                        else
                            imagename=['ch_',num2str(channel),'_sample_vs_test_act_',area,'_',num2str(sampleContrast)];
                        end
                        pathname=fullfile(rootFolder,'PL',analysisType,animal,area,imagename);
                        foldername=fullfile(rootFolder,'PL',analysisType,animal,area);
                        if ~exist(foldername,'dir')
                            mkdir(foldername)
                        end
                        printtext=sprintf('print -dpng %s.png',pathname);
                        set(gcf,'PaperPositionMode','auto')
                        eval(printtext);
                    end
                    for subplotMultiple=1:length(conds(sampleContrastsInd,:))%1:3
                        cond=conds(sampleContrastsInd,subplotMultiple);
                        [coefficients(subplotMultiple,1,h) coefficients(subplotMultiple,2,h)]=corr(AUROC(cond,:)',[1:length(sessionNums)]','type','Spearman');%r and h values
                    end
                    if plotFigs==1
                        close(figCh);
                        figure(figSess);
                        subplot(ceil(length(channels)/5),5,h);
                        for subplotMultiple=1:length(conds(sampleContrastsInd,:))%1:3
                            if coefficients(subplotMultiple,2,h)<0.05
                                markerType='o';
                            else
                                markerType='x';
                            end
                            plot(1:length(sessionNums),AUROC(cond,:),markerType,'Color',colmapText(subplotMultiple,:));
                            hold on
                            plot([0 length(sessionNums)+1],[0.5 0.5],'k--');
                            xlim([0 length(sessionNums)+1]);
                        end
                        yLimVals=get(gca,'YLim');
                        xLimVals=get(gca,'XLim');
                        allPlotLims(h,:)=[xLimVals yLimVals];
                        title(num2str(channel));
                        for subplotMultiple=1:length(conds(sampleContrastsInd,:))
                            if coefficients(subplotMultiple,2,h)<0.05
                                pCol=colmapText(subplotMultiple,:);
                            else
                                pCol='k';
                            end
                            ptext=sprintf('   r= %f',round(10000*coefficients(subplotMultiple,1,h))/10000);
                            text('Position',[xLimVals(1) yLimVals(2)+0.1*subplotMultiple*(yLimVals(2)-yLimVals(1))],'FontSize',8,'String',ptext,'Color',pCol);
                            ptext=sprintf('p= %f',round(10000*coefficients(subplotMultiple,2,h))/10000);
                            text('Position',[xLimVals(2)-(xLimVals(2)-xLimVals(1))/3 yLimVals(2)+0.1*subplotMultiple*(yLimVals(2)-yLimVals(1))],'FontSize',8,'String',ptext,'Color',pCol);
                        end
                    end
                    allCondAUROC=[];
                    for subplotMultiple=1:length(conds(sampleContrastsInd,:))%1:3
                        allCondAUROC=[allCondAUROC AUROC(conds(subplotMultiple),:)];
                        allAUROC{subplotMultiple}=[allAUROC{subplotMultiple};AUROC(conds(subplotMultiple),:)];%store AUROC values across channels
                        [hs(subplotMultiple,h) ps(subplotMultiple,h) ci stat]=ttest(allCondAUROC,0.5);
                        cis{subplotMultiple}=[cis{subplotMultiple};ci];
                        stats{subplotMultiple}=[stats{subplotMultiple};stat];
                    end
                end
                matName=['all_ch_sample_test_3conds_ROCs_',area];
                matPathName=fullfile(rootFolder,'PL',analysisType,animal,'coef_mat',matName);
                matFolderName=fullfile(rootFolder,'PL',analysisType,animal,'coef_mat');
                if ~exist(matFolderName,'dir')
                    mkdir(matFolderName);
                end
                saveText=['save ',matPathName,' hs ps cis stats coefficients allAUROC'];
                eval(saveText);
                
                if plotFigs==1
                    for h=1:length(channels)
                        subplot(ceil(length(channels)/5),5,h);
                        ylim([allPlotLims(h,3) allPlotLims(h,4)]);
                    end
                    figure(figSess);
                    set(gcf,'PaperPositionMode','auto')
                    subFolderName=[analysisType,'_coef_images'];
                    coefImagename=[analysisType,'_coefs_',area,'_',num2str(sampleContrast)];
                    coefFolder=fullfile('F:','PL',analysisType,animal,subFolderName);
                    if ~exist(coefFolder,'dir')
                        mkdir(coefFolder)
                    end
                    coefPathname=fullfile(coefFolder,coefImagename);
                    printtext=sprintf('print -dpng %s.png',coefPathname);
                    eval(printtext);
                end
            end
        end
    end
end

allStats=[];
alltStats=[];
allChtStats=[];
areas=[{'v4_1'} {'v1_1'}];% {'v1_2_1'} {'v1_2_2'} {'v1_2_3'}
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        channels=main_channels(animal,area);
        sessionNums = main_raw_sessions_final(animal,area,[],0);
        [sampleContrasts testContrasts]=area_metadata(area);
        if strcmp(area,'v4_1')
            conds=[8 9 10];%test contrast of 31, 32, 33%
        elseif strcmp(area,'v1_1')
            conds=8;%test contrast of 32%
        elseif strcmp(area,'v1_2_1')||strcmp(area,'v1_2_2')||strcmp(area,'v1_2_3')
            conds=[4;9;8];%test contrast of 22, 35, 42% for respective sample contrast conditions
        end
        colormap(cool(length(sessionNums)));
        for sampleContrastsInd=1:length(sampleContrasts)
            testContrast=testContrasts(sampleContrastsInd,:);
            matName=['all_ch_sample_test_3conds_ROCs_',area];
            matPathName=fullfile(rootFolder,'PL',analysisType,animal,'coef_mat',matName);
            matFolderName=fullfile(rootFolder,'PL',analysisType,animal,'coef_mat');
            loadText=['load ',matPathName];
            eval(loadText);
            for subplotMultiple=1:length(conds(sampleContrastsInd,:))%1:3
                if animalInd==1&&areaInd==1
                    figChs=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
                    set(figChs,'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                    figInd(subplotMultiple)=figChs;
                else
                    figure(figInd(subplotMultiple));
                end
                subplot(2,2,animalInd+(areaInd-1)*2);
                for chInd=1:length(channels)
                    chCol=[1-(chInd)/length(channels) 0 (chInd)/length(channels)];
                    plot(1:length(sessionNums),allAUROC{subplotMultiple}(chInd,:),'Marker','o','MarkerFaceColor',chCol,'LineStyle','none','Color',chCol);
                    hold on
                end
                xlim([0 length(sessionNums)+1]);
                plot([0 length(sessionNums)+1],[0.5 0.5],'k--');
                title(sprintf('%s %s %d',animal,area,conds(subplotMultiple)));
                %plot mean AUROC (Averaged across channels) against time:
                meanAUROC=mean(allAUROC{subplotMultiple}(:,:));
                if animalInd==1&&areaInd==1&&subplotMultiple==1
%                     figMeanChs=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
                    figMeanChs=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.5, 0.4]); %
                    set(figMeanChs,'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                else
                    figure(figMeanChs);
                end
%                 subplot(2,2,animalInd+(areaInd-1)*2);
                subplot(1,2,areaInd);
                if strcmp(area,'v4_1')
                    condCol=[1 0 0;0.05 0.73 0.24;0 0 1];
                else
                    condCol=[0.05 0.73 0.24];
                end
                animalMarkerFill=[{'none'} {condCol(subplotMultiple,:)}];
                plot(1:length(sessionNums),meanAUROC,'Marker','o','MarkerFaceColor',animalMarkerFill{animalInd},'LineStyle','none','Color',condCol(subplotMultiple,:));
                hold on
                xlim([0 length(sessionNums)+1]);
                plot([0 length(sessionNums)+1],[0.5 0.5],'k--');
%                 if animalInd+(areaInd-1)*2==2
                if animalInd+(areaInd-1)*2==1
                    subplot(1,2,2);
                    hold on
                    ptext=sprintf('- %d%% test',testContrast(conds(subplotMultiple)));
                    text('Position',[14 0.455+0.005*subplotMultiple],'FontSize',14,'String',ptext,'Color',condCol(subplotMultiple,:));
%                     text('Position',[18 0.44+0.015*subplotMultiple],'FontSize',14,'String',ptext,'Color',condCol(subplotMultiple,:));
                end
                [r p]=corr([1:length(sessionNums)]',meanAUROC','type','Spearman');
                allStats=[allStats;{animal} {area} testContrast(conds(subplotMultiple)) r p meanAUROC];
                notPoint5=find(ps(subplotMultiple,:)<0.05);
                numNotPoint5=length(notPoint5);
                lessPoint5=sum(cis{subplotMultiple}(notPoint5,1)<0.5);%number of channels with AUROC values with mean of significantly of less than 0.5
                morePoint5=sum(cis{subplotMultiple}(notPoint5,1)>0.5);%number of channels with AUROC values with mean of significantly of more than 0.5
                alltStats=[alltStats;{animal} {area} testContrast(conds(subplotMultiple)) numNotPoint5 lessPoint5 morePoint5];
                [h p ci stat]=ttest(allAUROC{subplotMultiple}',0.5);
                allchAUROC=[];
                for chInd=1:length(channels)
                    allchAUROC=[allchAUROC allAUROC{subplotMultiple}(chInd,:)];
                end
                [h p ci stat]=ttest(allchAUROC,0.5);
                allChtStats=[allChtStats;{animal} {area} testContrast(conds(subplotMultiple)) h p ci(1) ci(2) {stat}];
            end
        end
    end
end
% figure(figMeanChs);
% subplot(2,2,1);title('Monkey 1');ylabel('V4');
% subplot(2,2,2);ylim([0.32 0.51]);title('Monkey 2');
% subplot(2,2,3);ylabel('V1');
% subplot(2,2,4);ylim([0.41 0.51])
allStats
alltStats
allChtStats
allChtStatsTable=[];
allChtStatsCITable=[];
for row=1:8
    allChtStatsCITable=[allChtStatsCITable;allChtStats{row,8}.df allChtStats{row,8}.tstat {sprintf('%f - %f',allChtStats{row,6},allChtStats{row,7})} allChtStats{row,5}];
end
for row=1:8
    allChtStatsTable=[allChtStatsTable;allChtStats{row,8}.df allChtStats{row,8}.tstat allChtStats{row,5}];
end