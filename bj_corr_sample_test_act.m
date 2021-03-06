function bj_corr_sample_test_act(animals)
%written by Xing 12/04/13
%Calculates correlations in activity to sample and test

saveSampleTestAct=1;
comparisonType=1;
switch(comparisonType)
    case(1)%all trials
        analysisType='ROC';
    case(2)%only correct trials
        analysisType='ROC_correct_only';
end
onExternalHD=0;
if onExternalHD==1
    rootFolder='G:\PL_backup_060413';
else
    rootFolder='F:';
end
plotDual=0;
plotFigs=0;
plotFigs2=0;
if nargin<1||isempty(animals)
    animals=[{'blanco'} {'jack'}];
end
areas=[{'v4_1'} {'v4_2'} {'v1_1'} {'v1_2'}];
% areas=[{'v4_2'} {'v1_2'}];
areas=[{'v1_2'} {'v1_2_1'} {'v1_2_2'}];
areas=[{'v4_1'}];
test_epochs={0 512 512*2 512*3};durSpon=150;
durSpon=150;%length of period prior to sample onset from which spontaneous rates are calculated. Can take on a value of up to 512 ms.
minTrials=10;%set value of minumum number of trials for inclusion of session
subPeriod=1;
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        channels=main_channels(animal,area);
        sessionNums = main_raw_sessions_final(animal,area,[],0);
        cellEpochTimes={0 512 512*2 512*3};%{[0 40 300] 529 [529*2 529*2+40 529*2+300] 529*3}
        [sampleContrasts testContrasts]=area_metadata(area);
        if plotFigs==1
            figSess=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
            set(figSess, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
        end
        for h=1:length(channels)
            channel=channels(h);
            for sampleContrastsInd=1:length(sampleContrasts)
                sampleContrast=sampleContrasts(sampleContrastsInd);
                testContrast=testContrasts(sampleContrastsInd,:);
                colmapText=colormap(jet(size(testContrast,2)));
                if plotFigs==1
                    colmapText=[colmapText(1,:);132/255 22/255 216/255;202/255 65/255 223/255;colmapText(3:7,:);157/255 212/255 61/255;colmapText(10:12,:);178/255 111/255 12/255;colmapText(end,:)];
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
                        for cond=1:size(matarray,1)
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
                        end
                        if saveSampleTestAct==1
                            if comparisonType==1
                                saveFolderName=['F:\PL\sample_test_activity\',animal,'_',area];
                            elseif comparisonType==2
                                saveFolderName=['F:\PL\sample_test_activity\',animal,'_',area,'\plotRedArtifacts\correct_trials_only'];
                            end
                            if ~exist(saveFolderName,'dir')
                                mkdir(saveFolderName);
                            end
                            if strcmp(area,'v4_1')||strcmp(area,'v1_1')
                                saveText=['save F:\PL\sample_test_activity\',animal,'_',area,'\ch',num2str(channels(h)),'_',num2str(sessionNums(i)),'_example_sample_test_act.mat epoch2 epoch4'];
                            else
                                saveText=['save F:\PL\sample_test_activity\',animal,'_',area,'\ch',num2str(channels(h)),'_',num2str(sessionNums(i)),'_',num2str(sampleContrast),'_example_sample_test_act.mat epoch2 epoch4'];
                            end
                            eval(saveText);
                            matActName=['F:\PL\ch',num2str(channels(h)),'_',num2str(sessionNums(i)),'_example_sample_test_act.mat'];
%                             txtActName=['F:\PL\ch',num2str(channels(h)),'_',num2str(sessionNums(i)),'_example_sample_test_act.txt'];
%                             movefile(matActName,txtActName);
                        end
                        if plotFigs==1
                            figure(figCh);
                            if plotDual==1
                                subplot(ceil(length(sessionNums)/5),10,i*2-1);
                            else
                                subplot(ceil(length(sessionNums)/5),5,i);
                            end
                            for cond=1:size(matarray,1)
                                plot(epoch4{cond,1},epoch2{cond,1},'x','Color',colmapText(cond,:));
                                hold on
                            end
                        end
                        diff=[];
                        for cond=1:size(matarray,1)
                            [r(i,cond),p(i,cond)]=corr(epoch2{cond,1}',epoch4{cond,1}');
                            meanDiff(i,cond)=mean(epoch4{cond,1}-epoch2{cond,1});%test minus sample
                            diffMeanAll(i,cond)=mean(epoch4{cond,1})-mean(epoch2{cond,1});%test minus sample
                        end
                        if plotDual==1&&plotFigs==1
                            subplot(ceil(length(sessionNums)/5),10,i*2);
                            for cond=1:size(matarray,1)
                                plot(cond,r(i,cond),'x','Color',colmapText(cond,:));
                                hold on
                            end
                            xlim([0 15]);
                        end
                    end
                end
                if plotFigs2==1
                    if h==1
                        figDiffsDiagonal=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
                        set(figDiffsDiagonal, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                    end
                    figure(figDiffsDiagonal)
                    for i=1:length(sessionNums)
                        for cond=1:length(testContrast)
                            plot(meanDiff(i,cond),diffMeanAll(i,cond),'Marker','o','LineStyle','none','Color',colmapText(cond,:));hold on
                        end
                    end
                    if h==length(testContrast)
                        axis square
                        xlims=get(gca,'XLim');
                        ylims=get(gca,'YLim');
                        squarelims=[min([xlims(1) ylims(1)]) max([xlims(2) ylims(2)])];
                        line([squarelims(1) squarelims(2)],[squarelims(1) squarelims(2)],'Color','k','LineStyle',':');
                    end
                    if h==1
                        figDiffs=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
                        set(figDiffs, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                    end
                    figure(figDiffs)
                    diffDiffs=meanDiff-diffMeanAll;
                    for i=1:length(sessionNums)
                        for cond=1:length(testContrast)
                            plot(cond,diffDiffs(i,cond),'Marker','o','LineStyle','none','Color',colmapText(cond,:));hold on
                        end
                    end
                    if h==length(testContrast)
                        axis square
                        xlims=get(gca,'XLim');
                        ylims=get(gca,'YLim');
                        squarelims=[min([xlims(1) ylims(1)]) max([xlims(2) ylims(2)])];
                        line([0 14],[0 0],'Color','k','LineStyle',':');
                    end
                end
                if plotFigs==1
                    if plotDual==1
                        imagename=['ch_',num2str(channel),'test_vs_sample_act_',area,'_',num2str(sampleContrast),'_splitPlots'];
                    else
                        imagename=['ch_',num2str(channel),'test_vs_sample_act_',area,'_',num2str(sampleContrast)];
                    end
                    pathname=fullfile(rootFolder,'PL',analysisType,animal,'ROC_diff',imagename);
                    foldername=fullfile(rootFolder,'PL',analysisType,animal,'ROC_diff');
                    if ~exist(foldername,'dir')
                        mkdir(foldername)
                    end
                    printtext=sprintf('print -dpng %s.png',pathname);
                    set(gcf,'PaperPositionMode','auto')
                    eval(printtext);
                    close(figCh);
                end
                smallPs=sum(sum(p(:,:)<0.05));
                percPs(h)=smallPs/(length(sessionNums)*size(matarray,1));
                posRs=sum(sum(r(:,:)>0));
                percRs(h)=posRs/(length(sessionNums)*size(matarray,1));
                matName=[num2str(channel),'_',num2str(sampleContrast),'_Rp_ROCdiff_',area];
                matPathName=fullfile(rootFolder,'PL',analysisType,animal,'ROC_diff',matName);
                matFolderName=fullfile(rootFolder,'PL',analysisType,animal,'ROC_diff');
                if ~exist(matFolderName,'dir')
                    mkdir(matFolderName);
                end
                saveText=['save ',matPathName,' p r'];
                eval(saveText);
                if plotFigs==1
                    figure(figSess);
                    subplot(ceil(length(channels)/5),5,h);
                    for rowInd=1:size(r,1)
                        for columnInd=1:size(r,2)
                            if p(rowInd,columnInd)<0.05
                                markerType='o';
                            else
                                markerType='x';
                            end
                            plot(columnInd,r(rowInd,columnInd),markerType,'Color',colmapText(columnInd,:));
                            hold on
                        end
                    end
                end
            end
        end
        if plotFigs==1
            imagename=['Rs_test_vs_sample_',area,'_',num2str(sampleContrast)];
            pathname=fullfile(rootFolder,'PL',analysisType,animal,'ROC_diff',imagename);
            printtext=sprintf('print -dpng %s.png',pathname);
            set(gcf,'PaperPositionMode','auto')
            eval(printtext);
        end
    end
end

