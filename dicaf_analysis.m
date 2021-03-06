function dicaf_analysis
%Written by Xing 06/02/14.
%Calculates noise correlations across channels, for each trial, with a
%different (and I think incorrect) method- calculates z-scores based on
%each channel, rather than across channels.

sampleTestSamePlot=1;
onExternalHD=0;
if onExternalHD==1
    rootFolder='G:\PL_backup_060413';
else
    rootFolder='F:';
end
originalSessList=1;
animals=[{'blanco'} {'jack'}];
areas=[{'v4_1'} {'v1_1'}];
writeMeanVar=0;
if writeMeanVar==1%calculate mean and var statistics for real data
    for animalInd=1:length(animals)
        animal=animals{animalInd};
        for areaInd=1:length(areas)
            area=areas{areaInd};
            allSessR=[];
            folderName=['F:\PL\sample_test_activity\',animal,'_',area];
            channels = main_channels(animal,area);
            sessionNums = main_raw_sessions_final(animal,area,[],0);
            [sampleContrasts testContrasts]=area_metadata(area);
            if originalSessList==1
                if strcmp(animal,'blanco')&&strcmp(area,'v4_1')
                    sessionNums=[307 308 311 313 314 318 320 321 329 330 331:1:341];% blanco V4
                end
            end
            for cond=1:length(testContrasts)
                figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.05, 0.01, 0.9, 0.9]);
                for sessionInd=1:length(sessionNums)
                    if sessionNums(sessionInd)~=322
                        for chInd=1:length(channels)
                            fileName=['ch',num2str(channels(chInd)),'_',num2str(sessionNums(sessionInd)),'_example_sample_test_act'];
                            matName=fullfile(folderName,fileName);
                            load(matName);
                            chMeanActS(chInd,sessionInd,cond)=mean(epoch2{cond});%calculate mean sample-evoked response for each channel and session
                            chVarActS(chInd,sessionInd,cond)=var(epoch2{cond});
                            chMeanActT(chInd,sessionInd,cond)=mean(epoch4{cond});%calculate mean test-evoked response for each channel and session
                            chVarActT(chInd,sessionInd,cond)=var(epoch4{cond});
                        end
                    end
                end
                chSessCount=1;
                for chInd=1:length(channels)
                    for sessInd=1:length(sessionNums)
                        chSessCount=chSessCount+0.1;
                        plot(chSessCount,chMeanActS(chInd,sessInd,cond),'ro','MarkerSize',5,'MarkerFaceColor','r');hold on
                        plot([chSessCount chSessCount],[chMeanActS(chInd,sessInd,cond)-sqrt(chVarActS(chInd,sessionInd,cond)) chMeanActS(chInd,sessInd,cond)+sqrt(chVarActS(chInd,sessionInd,cond))],'r');
                        if sampleTestSamePlot==1
                            plot(chSessCount+0.05,chMeanActT(chInd,sessInd,cond),'bo','MarkerSize',5,'MarkerFaceColor','b');hold on
                            plot([chSessCount+0.05 chSessCount+0.05],[chMeanActT(chInd,sessInd,cond)-sqrt(chVarActT(chInd,sessionInd,cond)) chMeanActT(chInd,sessInd,cond)+sqrt(chVarActT(chInd,sessionInd,cond))],'b');
                        end
                    end
                    plot(chSessCount,chMeanActS(chInd,sessInd,cond),'ko','MarkerSize',5,'MarkerFaceColor','k');hold on
                end
                title([animal,' ',area]);
                %         subplot(2,2,animalInd+(areaInd-1)*2);
                xlabel('channel and session number');
                ylabel('mean sample act (spikes/s)');
                printFigs=1;
                if printFigs==1
                    folderPrint=fullfile(rootFolder,'PL','AUROC_vs_PROBMAT','real_data');
                    formats=[{'epsc'} {'png'}];
                    for i=1:length(formats)
                        format=formats{i};
                        if sampleTestSamePlot==0
                            figName=[folderPrint,'\',animal,'_',area,'_mean_sample_act_cond_',num2str(cond)];
                        elseif sampleTestSamePlot==1
                            figName=[folderPrint,'\',animal,'_',area,'_mean_sample_test_act_cond_',num2str(cond)];
                        end
                        printtext=sprintf('print -d%s %s -r600',format,figName);
                        set(gcf, 'PaperPositionMode', 'auto');
                        eval(printtext);
                    end
                end
                if sampleTestSamePlot==0
                    chSessCount=1;
                    figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.05, 0.01, 0.9, 0.9]);
                    for chInd=1:length(channels)
                        for sessInd=1:length(sessionNums)
                            chSessCount=chSessCount+0.1;
                            plot(chSessCount,chMeanActT(chInd,sessInd,cond),'ko','MarkerSize',5,'MarkerFaceColor','b');hold on
                            plot([chSessCount chSessCount],[chMeanActT(chInd,sessInd,cond)-sqrt(chVarActT(chInd,sessionInd,cond)) chMeanActT(chInd,sessInd,cond)+sqrt(chVarActT(chInd,sessionInd,cond))],'r');
                        end
                        plot(chSessCount,chMeanActT(chInd,sessInd,cond),'ko','MarkerSize',5,'MarkerFaceColor','k');hold on
                    end
                    title([animal,' ',area]);
                    xlabel('channel and session number');
                    ylabel('mean test act (spikes/s)');
                    printFigs=1;
                    if printFigs==1
                        folderPrint=fullfile(rootFolder,'PL','AUROC_vs_PROBMAT','real_data');
                        formats=[{'epsc'} {'png'}];
                        for i=1:length(formats)
                            format=formats{i};
                            figName=[folderPrint,'\',animal,'_',area,'_mean_test_act_cond_',num2str(cond)];
                            printtext=sprintf('print -d%s %s -r600',format,figName);
                            set(gcf, 'PaperPositionMode', 'auto');
                            eval(printtext);
                        end
                    end
                end
            end
            saveText=['save F:\PL\AUROC_vs_PROBMAT\real_data\',animal,'_',area,'_mean_var.mat chMeanActS chMeanActT chVarActS chVarActT'];
            eval(saveText);
        end
    end
end

%create dummy data for 2 simulated channels
animal='blanco';
area='v4_1';
loadText=['load F:\PL\AUROC_vs_PROBMAT\real_data\',animal,'_',area,'_mean_var.mat chMeanActS chMeanActT chVarActS chVarActT'];
eval(loadText);
cond=12;
ch1=1;%set example ch to be first ch
ch2=2;%set example ch to be second ch
sess=1;%example session

%select means and variances for simulated chs
simChMeanActS1=chMeanActS(ch1,sess,cond);
simChMeanActS2=chMeanActS(ch2,sess,cond);
diffST1=chMeanActT(ch1,:,cond)-chMeanActS(ch1,:,cond);%difference between sample- & test-evoked activities
diffST2=chMeanActT(ch2,:,cond)-chMeanActS(ch2,:,cond);%
varDiffST1=var(diffST1);%calculate variance of differences in means for first ch
varDiffST2=var(diffST2);%calculate variance of differences in means for second ch
numSimTrials=100;
numSimRuns=1000;
simDiff1=normrnd(mean(diffST1),sqrt(varDiffST1),numSimTrials,numSimRuns);%generate random sample-test difference in activity levels
simDiff2=normrnd(mean(diffST2),sqrt(varDiffST2),numSimTrials,numSimRuns);%

% simChMeanActT1=chMeanActT(ch1,sess,cond);
% simChMeanActT2=chMeanActT(ch2,sess,cond);

simChVarActS1=chVarActT(ch1,sess,cond);
simChVarActS2=chVarActT(ch2,sess,cond);
% simChVarActT1=chVarActT(ch1,sess,cond);
% simChVarActT2=chVarActT(ch2,sess,cond);

simActS1=normrnd(simChMeanActS1,sqrt(simChVarActS1),numSimTrials,numSimRuns);%generate random sample-evoked activity levels for sim ch 1
simActS2=normrnd(simChMeanActS2,sqrt(simChVarActS2),numSimTrials,numSimRuns);%generate random sample-evoked activity levels for sim ch 2
simActT1=simActS1+simDiff1;%generate simulated test-evoked act by adding a simulated difference value
simActT2=simActS2+simDiff2;
% simActT1=normrnd(simChMeanActT1,sqrt(simChVarActT1),numSimTrials);%generate random sample-evoked activity levels for sim ch 1
% simActT2=normrnd(simChMeanActT2,sqrt(simChVarActT2),numSimTrials);%generate random sample-evoked activity levels for sim ch 2
simActST1=[simActS1;simActT1];%combine sample- and test-evoked act together to calculate R between sim chs 1 & 2
simActST2=[simActS2;simActT2];
for runInd=1:numSimRuns
    tempR=corrcoef([simActST1(:,runInd) simActST2(:,runInd)]);%calculate R (correlation) between the 2 sim chs, for each run
    r12(runInd)=tempR(2);
    AUROC1(runInd)=sglroc3(simActT1(:,runInd)',simActS1(:,runInd)');
    AUROC2(runInd)=sglroc3(simActT2(:,runInd)',simActS2(:,runInd)');
    AUROC12(runInd)=sglroc3(simActT1(:,runInd)'+simActT2(:,runInd)',simActS1(:,runInd)'+simActS2(:,runInd)');%add act across both sim chs, to calculate AUROC on combined ch act
    higher1=sum(simActT1(:,runInd)>simActS1(:,runInd));%number of trials where test act > sample act for ch 1
    higher2=sum(simActT2(:,runInd)>simActS2(:,runInd));%number of trials where test act > sample act for ch 2
    lower1=sum(simActT1(:,runInd)<simActS1(:,runInd));%number of trials where test act < sample act for ch 1
    lower2=sum(simActT2(:,runInd)<simActS2(:,runInd));%number of trials where test act < sample act for ch 2
    DICAF1(runInd)=higher1/(higher1+lower1);%calculate DICAF/PROBMAT/COBAM value
    DICAF2(runInd)=higher2/(higher2+lower2);
    higher12=sum((simActT1(:,runInd)+simActT2(:,runInd))>(simActS1(:,runInd)+simActS2(:,runInd)));%number of trials where test act > sample act for combined ch activities
    lower12=sum((simActT1(:,runInd)+simActT2(:,runInd))<(simActS1(:,runInd)+simActS2(:,runInd)));%number of trials where test act > sample act for combined ch activities
    DICAF12(runInd)=(higher12)/(higher12+lower12);%add act across both sim chs, to calculate AUROC on combined ch act
end
figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.05, 0.01, 0.9, 0.9]);%population data (combined across chs)
plot(r12,AUROC12,'ko');hold on
plot(r12,DICAF12,'ko','MarkerFaceColor','k');
title('simulated population data');
xlabel('R between two chs');
ylabel('AUROC or DICAF value');
printFigs2=1;
if printFigs2==1
    folderPrint=fullfile(rootFolder,'PL','AUROC_vs_PROBMAT','simulated_data');
    formats=[{'epsc'} {'png'}];
    for i=1:length(formats)
        format=formats{i};
        figName=[folderPrint,'\',animal,'_',area,'_AUROCDICAF_vs_R_chs',num2str(ch1),'_and_',num2str(ch2),'_sess',num2str(sess),'_cond_',num2str(cond)];
        printtext=sprintf('print -d%s %s -r600',format,figName);
        set(gcf, 'PaperPositionMode', 'auto');
        eval(printtext);
    end
end

samePlot=1;
figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.05, 0.01, 0.9, 0.9]);%individual ch data
if samePlot==1
    title('simulated individual ch 1 data');
    xlabel('AUROC');
    ylabel('DICAF');
else
    subplot(1,2,1);
end
plot(AUROC1,DICAF1,'ko');hold on
if samePlot==1
    plot(AUROC2,DICAF2,'ro');
    title('simulated individual ch data');
    xlabel('AUROC');
    ylabel('DICAF');
else
    subplot(1,2,2);
    plot(AUROC2,DICAF2,'ko');
    title('simulated individual ch 2 data');
    xlabel('AUROC');
    ylabel('DICAF');
end
xlimvals=get(gca,'xlim');
ylimvals=get(gca,'ylim');
% lineOfEquality=[max([xlimvals(1) ylimvals(1)]) min([xlimvals(2) ylimvals(2)])];
% plot([lineOfEquality(1) lineOfEquality(2)],[lineOfEquality(1) lineOfEquality(2)],'LineStyle','--','Color','k');
plot([0.5 0.5],[ylimvals(1) ylimvals(2)],'LineStyle','--','Color','k');
xlim(xlimvals);

lineOfEquality=[min([xlimvals(1) ylimvals(1)]) max([xlimvals(2) ylimvals(2)])];
plot([lineOfEquality(1) lineOfEquality(2)],[lineOfEquality(1) lineOfEquality(2)],'LineStyle','--','Color','k');
xlim(lineOfEquality);
ylim(lineOfEquality);
axis square
printFigs3=1;
if printFigs3==1&&samePlot==1
    folderPrint=fullfile(rootFolder,'PL','AUROC_vs_PROBMAT','simulated_data');
    formats=[{'epsc'} {'png'}];
    for i=1:length(formats)
        format=formats{i};
        figName=[folderPrint,'\',animal,'_',area,'_AUROC_vs_DICAF_chs',num2str(ch1),'_and_',num2str(ch2),'_sess',num2str(sess),'_cond_',num2str(cond)];
        printtext=sprintf('print -d%s %s -r600',format,figName);
        set(gcf, 'PaperPositionMode', 'auto');
        eval(printtext);
    end
end

saveForLater=1;
if saveForLater==0
    
simRangeS=[min(min(chMeanActS(:,:,cond))) max(max(chMeanActS(:,:,cond)))];%find range for mean sample-evoked act, across chs
simRangeT=[min(min(chMeanActT(:,:,cond))) max(max(chMeanActS(:,:,cond)))];%find range for mean test-evoked act, across chs



for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        allSessR=[];
        folderName=['F:\PL\sample_test_activity\',animal,'_',area];
        channels = main_channels(animal,area);
        sessionNums = main_raw_sessions_final(animal,area,[],0);
        [sampleContrasts testContrasts]=area_metadata(area);
        if originalSessList==1
            if strcmp(animal,'blanco')&&strcmp(area,'v4_1')
                sessionNums=[307 308 311 313 314 318 320 321 329 330 331:1:341];% blanco V4
            end
        end
        
        for sessionInd=1:length(sessionNums)
            if sessionNums(sessionInd)~=322
                zChAct=cell(1,length(channels));
                rvals=[];
                for chPair1=1:length(channels)-1
                    fileName=['ch',num2str(channels(chPair1)),'_',num2str(sessionNums(sessionInd)),'_example_sample_test_act'];
                    matName=fullfile(folderName,fileName);
                    load(matName);
                    ch1Act=epoch2;
                    for chPair2=chPair1+1:length(channels)
                        fileName=['ch',num2str(channels(chPair2)),'_',num2str(sessionNums(sessionInd)),'_example_sample_test_act'];
                        matName=fullfile(folderName,fileName);
                        load(matName);
                        ch2Act=epoch2;
                        zCh1Act=[];
                        zCh2Act=[];
                        for condInd=1:length(testContrasts)
                            zCh1Act=[zCh1Act zscore(ch1Act{condInd,1})];
                            zCh2Act=[zCh2Act zscore(ch2Act{condInd,1})];
                        end
                        tempR=corrcoef(zCh1Act,zCh2Act);
                        rvals=[rvals tempR(2)];%for each pairwise comparison between channels
                    end
                end
                allSessR(sessionInd,:)=rvals;%compile across sessions
            end
        end
        numCombinedSess=5;
%         numCombinedSess=floor(length(sessionNums)/2);
        early1=[];
        late1=[];
        for sessInd=1:numCombinedSess
            early1=[early1 allSessR(sessInd,:)];
            late1=[late1 allSessR(length(sessionNums)-sessInd+1,:)];
        end
        mean(early1)
        mean(late1)
        bins=[-0.95:0.01:0.95];
        subplot(2,2,animalInd+(areaInd-1)*2);
        hist(early1(:),bins);
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','r','facealpha',0.5)
        hold on
        hist(late1(:),bins);
        h = findobj(gca,'Type','patch');
        set(h,'facealpha',0.5)
        xlim([-0.3 0.5]);
        allSessR=[sessionNums' allSessR];
        saveText=['save F:\PL\noise_trial_corr\',animal,'_',area,'_noise_R.mat allSessR'];
        eval(saveText);
    end
end

end