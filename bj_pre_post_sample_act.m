function bj_pre_post_sample_act(animals,readData)
%written by Xing 12/04/13
%Calculates correlations in activity to pre-sample and pre-test

analysisTypeText='preStim';
saveSampleTestAct=1;
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
% areas=[{'v4_2'} {'v1_2'}];
% areas=[{'v4_1'} {'v1_1'}];
test_epochs={0 512 512*2 512*3};durSpon=150;
durSpon=150;%length of period prior to sample onset from which spontaneous rates are calculated. Can take on a value of up to 512 ms.
minTrials=10;%set value of minumum number of trials for inclusion of session
subPeriod=1;
if readData
    for animalInd=1:length(animals)
        animal=animals{animalInd};
        for areaInd=1:length(areas)
            area=areas{areaInd};
            channels=main_channels(animal,area);
            sessionNums = main_raw_sessions_final(animal,area,[],0);
            cellEpochTimes={0 512 512*2 512*3};%{[0 40 300] 529 [529*2 529*2+40 529*2+300] 529*3}
            [sampleContrasts testContrasts]=area_metadata(area);
            for sampleContrastsInd=1:length(sampleContrasts)
                sampleContrast=sampleContrasts(sampleContrastsInd);
                testContrast=testContrasts(sampleContrastsInd,:);
                if plotFigs==1
                    figSess=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
                    set(figSess, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                end
                AUROC=[];
                hs=[];
                ps=[];
                cis=[];
                stats={[]};
                coefficients=[];
                for h=1:length(channels)
                    channel=channels(h);
                    colmapText=colormap(jet(size(testContrast,2)));
                    if plotFigs==1
                        colmapText=[colmapText(1,:);132/255 22/255 216/255;202/255 65/255 223/255;colmapText(3:7,:);157/255 212/255 61/255;colmapText(10:12,:);178/255 111/255 12/255;colmapText(end,:)];
                        figCh=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
                        set(figCh, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                    end
                    for i=1:length(sessionNums)
                        matFolder=['F:\PL\spikeData\',animal];
                        chStr=[num2str(channel),'_',num2str(sessionNums(i)),'_',num2str(sampleContrast),'.mat'];
                        matPath=fullfile(matFolder,chStr);
                        matExists=0;
                        if exist(matPath,'file')
                            matExists=1;
                        end
                        if matExists==1
                            valsText=['load ',matPath,' matarray'];
                            eval(valsText);
                            higherPost=0;%combine trials across conditions
                            higherPre=0;
                            preCount=1;
                            postCount=1;
                            for cond=1:size(matarray,1)
                                trials1=length(matarray{cond,1});
                                trials3=length(matarray{cond,3});
                                actList=[];
                                for n=1:min([trials1 trials3])
                                    temp1=matarray{cond,1}{n}>test_epochs{1}-256;%activity during sample presentation
                                    spikes=matarray{cond,1}{n}(temp1);
                                    actList1(preCount)=length(spikes)*1000/256;
                                    temp3=matarray{cond,3}{n}>test_epochs{3}-256;%activity during ISI
                                    spikes=matarray{cond,3}{n}(temp3);
                                    actList3(postCount)=length(spikes)*1000/256;
                                    if actList1(preCount)<actList3(postCount)
                                        higherPost=higherPost+1;
                                    elseif actList1(preCount)>actList3(postCount)
                                        higherPre=higherPre+1;
                                    end
                                    preCount=preCount+1;
                                    postCount=postCount+1;
                                end
                                if preCount~=postCount
                                    pauseHere=1;
                                end
                            end
                            if saveSampleTestAct==1
                                saveFolderName=['F:\PL\',analysisTypeText,'\',animal,'\',area];
                                if ~exist(saveFolderName,'dir')
                                    mkdir(saveFolderName);
                                end
                                saveText=['save F:\PL\',analysisTypeText,'\',animal,'\',area,'\ch',num2str(channels(h)),'_',num2str(sessionNums(i)),'_example_sample_test_act.mat actList1 actList3 higherPost higherPre'];
                                eval(saveText);
                            end
                            if plotFigs==1
                                figure(figCh);
                                if plotDual==1
                                    subplot(ceil(length(sessionNums)/5),10,i*2-1);
                                    if i*2-1==1
                                        ylabel('pre-sample');
                                        xlabel('pre-test');hold on
                                    end
                                else
                                    subplot(ceil(length(sessionNums)/5),5,i);
                                    if i==1
                                        xlabel('pre-sample');
                                        ylabel('pre-test');hold on
                                    end
                                end
                                for cond=1:size(matarray,1)%combine across conditions, as activity occurs before test ever appears
                                    plot(actList1,actList3,'kx','MarkerSize',4);%pre-test activity against pre-sample activity ,'Color',colmapText(cond,:)
                                    hold on
                                end
                                xlims=get(gca,'XLim');
                                ylims=get(gca,'YLim');
                                newLims=[min([xlims ylims]) max([xlims ylims])];
                                xlim([newLims(1) newLims(2)]);
                                ylim([newLims(1) newLims(2)]);
                                plot([newLims(1) newLims(2)],[newLims(1) newLims(2)],'k--');
                                axis square
                                title(sessionNums(i));
                            end
                            [hs(h,i),ps(h,i),cis(h,i,:),stat]=ttest(actList3,actList1);
                            stats{h,i,:}=stat;
                            
                            %                         diff=[];
                            %                         for cond=1:size(matarray,1)
                            %                             meanDiff(i,cond)=mean(epoch4{cond,1}-epoch2{cond,1});%test minus sample
                            %                             diffMeanAll(i,cond)=mean(epoch4{cond,1})-mean(epoch2{cond,1});%test minus sample
                            %                         end
                            %                         if plotDual==1&&plotFigs==1
                            %                             subplot(ceil(length(sessionNums)/5),10,i*2);
                            %                             for cond=1:size(matarray,1)
                            %                                 plot(cond,r(i,cond),'x','Color',colmapText(cond,:));
                            %                                 hold on
                            %                             end
                            %                             xlim([0 15]);
                            %                         end
                        end
                        AUROC(h,i)=higherPost/(higherPost+higherPre);
                    end
                    if plotFigs==1
                        figure(figCh);
                        if plotDual==1
                            imagename=['ch_',num2str(channel),'_presample_vs_pretest_act_',area,'_splitPlots'];
                        else
                            imagename=['ch_',num2str(channel),'_presample_vs_pretest_act_',area];
                        end
                        pathname=fullfile(rootFolder,'PL',analysisTypeText,animal,area,imagename);
                        foldername=fullfile(rootFolder,'PL',analysisTypeText,animal,area);
                        if ~exist(foldername,'dir')
                            mkdir(foldername)
                        end
                        printtext=sprintf('print -dpng %s.png',pathname);
                        set(gcf,'PaperPositionMode','auto')
                        eval(printtext);
                        close(figCh);
                        figure(figSess);
                        subplot(ceil(length(channels)/5),5,h);
                        if length(AUROC)~=length(sessionNums)
                            unequal=1;
                        end
                    end
                    [coefficients(1,h) coefficients(2,h)]=corr(AUROC(h,:)',[1:length(sessionNums)]','type','Spearman');%r and h values
                    if plotFigs==1
                        if coefficients(2,h)<0.05
                            pCol='g';
                        else
                            pCol='k';
                        end
                        plot(1:length(sessionNums),AUROC(h,:),'ko');
                        hold on
                        plot([0 length(sessionNums)+1],[0.5 0.5],'k--');
                        xlim([0 length(sessionNums)+1]);
                        ptext=sprintf(' r= %f  p= %f',coefficients(1,h),coefficients(2,h));
                        yLimVals=get(gca,'YLim');
                        xLimVals=get(gca,'XLim');
                        text('Position',[xLimVals(1) yLimVals(1)+0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext,'Color',pCol);
                        title(channels(h));
                    end
                end
                matName=['preStim_',area,'_',num2str(sampleContrast)];
                matPathName=fullfile(rootFolder,'PL',analysisTypeText,animal,matName);
                matFolderName=fullfile(rootFolder,'PL',analysisTypeText,animal);
                if ~exist(matFolderName,'dir')
                    mkdir(matFolderName);
                end
                saveText=['save ',matPathName,' hs ps cis stats coefficients AUROC'];
                eval(saveText);
                if plotFigs==1
                    figure(figSess);
                    set(gcf,'PaperPositionMode','auto')
                    subFolderName=[analysisTypeText,'_coef_images'];
                    %             if excludeSessHighSSE==0
                    coefImagename=[analysisTypeText,'_coefs_',area,'_',num2str(sampleContrast)];
                    %             elseif excludeSessHighSSE==1
                    %                 if excludeOutliers==0
                    %                     coefImagename=[chText,appendText,startEndTime,'_',analysisTypeText,'_coefs_',area,'goodSSE'];
                    %                 elseif excludeOutliers==1
                    %                     coefImagename=[chText,appendText,startEndTime,'_',analysisTypeText,'_coefs_',area,'goodSSE_no_outliers_sl',num2str(slSigmaMultiple),'_C50',num2str(c50SigmaMultiple)];
                    %                 end
                    %             end
                    coefFolder=fullfile('F:','PL',analysisTypeText,animal,subFolderName);
                    if ~exist(coefFolder,'dir')
                        mkdir(coefFolder)
                    end
                    coefPathname=fullfile(coefFolder,coefImagename);
                    printtext=sprintf('print -dpng %s.png',coefPathname);
                    eval(printtext);
                end
            end
            %         if plotFigs==1
            %             imagename=['Rs_test_vs_sample_',area];
            %             pathname=fullfile(rootFolder,'PL',analysisType,animal,'ROC_diff',imagename);
            %             printtext=sprintf('print -dpng %s.png',pathname);
            %             set(gcf,'PaperPositionMode','auto')
            %             eval(printtext);
            %         end
        end
    end
end

justNonRoving=1;
if justNonRoving==1
    areas=[{'v4_1'} {'v1_1'}];
end
allStats=[];
alltStats=[];
allChtStats=[];
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        channels=main_channels(animal,area);
        sessionNums = main_raw_sessions_final(animal,area,[],0);
        cellEpochTimes={0 512 512*2 512*3};%{[0 40 300] 529 [529*2 529*2+40 529*2+300] 529*3}
        [sampleContrasts testContrasts]=area_metadata(area);
        for sampleContrastsInd=1:length(sampleContrasts)
            sampleContrast=sampleContrasts(sampleContrastsInd);
            testContrast=testContrasts(sampleContrastsInd,:);
            if plotFigs==1
                figSess=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
                set(figSess, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
            end
            matName=['preStim_',area,'_',num2str(sampleContrast)];
            matPathName=fullfile(rootFolder,'PL',analysisTypeText,animal,matName);
            matFolderName=fullfile(rootFolder,'PL',analysisTypeText,animal);
            loadText=['load ',matPathName,' hs ps cis stats coefficients AUROC'];
            eval(loadText);
            [h p ci stat]=ttest(AUROC',0.5)
            notPoint5=find(p<0.05);
            numNotPoint5=length(notPoint5);
            lessPoint5=sum(ci(1,notPoint5)<0.5);%number of channels with AUROC values with mean of significantly of less than 0.5
            morePoint5=sum(ci(1,notPoint5)>0.5);%number of channels with AUROC values with mean of significantly of more than 0.5
            alltStats=[alltStats;{animal} {area} numNotPoint5 lessPoint5 morePoint5];
            allchAUROC=[];
            allCI=[];
            stats=[];
            for chInd=1:length(channels)
                allchAUROC=[allchAUROC AUROC(chInd,:)];
                [hs(chInd) ps(chInd) ci stat]=ttest(AUROC(chInd,:)',0.5);
                allCI(chInd,:)=ci;
                stats{chInd}=stat;
            end
            [h p ci stat]=ttest(allchAUROC,0.5);
            allChtStats=[allChtStats;{animal} {area} h p ci(1) ci(2) {stat}];
            meanAUROC=mean(AUROC);
            if animalInd==1&&areaInd==1
                figMeanChs=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.35, 0.5]); %
                set(figMeanChs,'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
            else
                figure(figMeanChs);
            end
            subplot(1,2,areaInd);
            animalMarkerFill=[{'none'} {'k'}];
            plot(1:length(sessionNums),meanAUROC,'Marker','o','MarkerFaceColor',animalMarkerFill{animalInd},'LineStyle','none','Color','k');
            hold on
            xlim([0 length(sessionNums)+1]);
            plot([0 length(sessionNums)+1],[0.5 0.5],'k--');
            [r p]=corr([1:length(sessionNums)]',meanAUROC','type','Spearman');
            allStats=[allStats;{animal} {area} r p meanAUROC];
        end
    end
end
figure(figMeanChs);
subplot(1,2,1);title('V4');
subplot(1,2,2);title('V1');
allStats%corr analysis across channels, check population AUROC changes with time
alltStats%t-test analysis to see if mean of AUROC values differ from 0.5 for individual channels, tally numbers
allChtStats%t-test combining AUROC values across channels, check if population AUROC values differ from 0.5
allChtStatsTable=[];
allChtStatsCITable=[];
for row=1:size(allChtStats,1)
    allChtStatsCITable=[allChtStatsCITable;allChtStats{row,7}.df allChtStats{row,7}.tstat {sprintf('%f - %f',allChtStats{row,5},allChtStats{row,6})} allChtStats{row,4}];
end
for row=1:size(allChtStats,1)
    allChtStatsTable=[allChtStatsTable;allChtStats{row,7}.df allChtStats{row,7}.tstat allChtStats{row,4}];
end