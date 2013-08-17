function bj_pre_post_sample_act2(animals,readData)
%written by Xing 25/07/13
%Modified from bj_pre_post_sample_act, calculates PROBMAT values with
%activity across channels, instead of finding mean PROBMAT across
%individual channels.
%Calculates correlations in activity to pre-sample and pre-test

analysisTypeText='preStim';
saveSampleTestAct=1;
allChstStatsAct=[];
r=[];
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
% areas=[{'v1_2_1'} {'v1_2_2'} {'v1_2_3'}];
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
                PROBMAT=[];
                colmapText=colormap(jet(size(testContrast,2)));
                if plotFigs==1
                    colmapText=[colmapText(1,:);132/255 22/255 216/255;202/255 65/255 223/255;colmapText(3:7,:);157/255 212/255 61/255;colmapText(10:12,:);178/255 111/255 12/255;colmapText(end,:)];
                end
                %read activity across channels into large new matrix
                for h=1:length(channels)
                    channel=channels(h);
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
                            allChMat{h,i}=matarray;
                        end
                    end
                end
                preCount=1;
                postCount=1;
                actList1=zeros(1000000,1);
                actList3=zeros(1000000,1);
                for i=1:length(sessionNums)%combine activity across channels, the calculate grand PROBMAT across trials and conditions
                    higherPost=0;
                    higherPre=0;
                    allSumAct1=[];
                    allSumAct3=[];
                    for cond=1:size(matarray,1)
                        trials1=length(allChMat{1,i}{cond,1});%same number of trials for each channel
                        trials3=length(allChMat{1,i}{cond,3});
                        if trials1~=trials3
                            notEqual=1;
                        end
                        sumAct1=zeros(trials1,1);
                        sumAct3=zeros(trials3,1);
                        for n=1:min([trials1 trials3])
                            for h=1:length(channels)
                                channel=channels(h);
                                temp1=allChMat{h,i}{cond,1}{n}>test_epochs{1}-256;%sum up activity across chs during sample presentation for each trial
                                spikes=allChMat{h,i}{cond,1}{n}(temp1);
                                sumAct1(n)=sumAct1(n)+length(spikes)*1000/256;
                                actList1(preCount)=length(spikes)*1000/256;
                                temp3=allChMat{h,i}{cond,3}{n}>test_epochs{3}-256;%activity during ISI
                                spikes=allChMat{h,i}{cond,3}{n}(temp3);
                                sumAct3(n)=sumAct3(n)+length(spikes)*1000/256;
                                actList3(postCount)=length(spikes)*1000/256;
                                preCount=preCount+1;
                                postCount=postCount+1;
                            end
                        end
                        allSumAct1=[allSumAct1;sumAct1];%concatenate all trials, across conditions and sessions
                        allSumAct3=[allSumAct3;sumAct3];
                    end
                    %calculate PROBMAT values for each session:
                    for trialCount=1:length(allSumAct1)
                        if allSumAct1(trialCount)<allSumAct3(trialCount)%after summing up activity across channels, compare between pre-stim periods
                            higherPost=higherPost+1;
                        elseif allSumAct1(trialCount)>allSumAct3(trialCount)
                            higherPre=higherPre+1;
                        end
                    end
                    PROBMAT(i)=higherPost/(higherPost+higherPre);                    
                end
                %difference between activity across all session, trials,
                %conditions, based on population-wide activity:
                actList1=actList1(1:preCount);
                actList3=actList3(1:postCount);
                [hs,ps,cis,stat]=ttest(actList3,actList1);
                stats{:}=stat;
                allChstStatsAct=[allChstStatsAct;{animal} {area} stat.df stat.tstat ps];
                %check for changes in PROBMAT over time:
                [coef1 coef2]=corr(PROBMAT',[1:length(sessionNums)]','type','Spearman');%r and h values
                r=[r;coef1 coef2];
                matName=['preStim_',area,'_',num2str(sampleContrast),'_act'];
                matPathName=fullfile(rootFolder,'PL',analysisTypeText,animal,matName);
                matFolderName=fullfile(rootFolder,'PL',analysisTypeText,animal);
                if ~exist(matFolderName,'dir')
                    mkdir(matFolderName);
                end
                saveText=['save ',matPathName,' hs ps cis stats allChstStatsAct r PROBMAT'];
                eval(saveText);
            end
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
            matName=['preStim_',area,'_',num2str(sampleContrast),'_act'];
            matPathName=fullfile(rootFolder,'PL',analysisTypeText,animal,matName);
            matFolderName=fullfile(rootFolder,'PL',analysisTypeText,animal);
            loadText=['load ',matPathName,' hs ps cis stats allChstStatsAct r PROBMAT'];
            eval(loadText);
            [h p ci stat]=ttest(PROBMAT',0.5)
            notPoint5=find(p<0.05);
            numNotPoint5=length(notPoint5);
            lessPoint5=sum(ci(1,notPoint5)<0.5);%number of channels with AUROC values with mean of significantly of less than 0.5
            morePoint5=sum(ci(1,notPoint5)>0.5);%number of channels with AUROC values with mean of significantly of more than 0.5
            alltStats=[alltStats;{animal} {area} numNotPoint5 lessPoint5 morePoint5];
            allchAUROC=[];
            allCI=[];
            stats=[];
            if animalInd==1&&areaInd==1
                figMeanChs=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.35, 0.5]); %
                set(figMeanChs,'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
            else
                figure(figMeanChs);
            end
            subplot(1,2,areaInd);
            animalMarkerFill=[{'none'} {'k'}];
            plot(1:length(sessionNums),PROBMAT,'Marker','o','MarkerFaceColor',animalMarkerFill{animalInd},'LineStyle','none','Color','k');
            hold on
            xlim([0 length(sessionNums)+1]);
            plot([0 length(sessionNums)+1],[0.5 0.5],'k--');
            [r p]=corr([1:length(sessionNums)]',PROBMAT','type','Spearman');
            allStats=[allStats;{animal} {area} r p PROBMAT];
        end
    end
end
figure(figMeanChs);
subplot(1,2,1);title('V4');
xlabel('session number');
ylabel('PROBMAT');
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