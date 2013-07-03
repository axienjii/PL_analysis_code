function bj_plot_trials_across_sessions(animals,areas)
%Written by Xing 01/07/2013.
%Reads activity and response correctness (correct or incorrect) for each
%trial, across sessions, from sample_test_activity & TAS folders
%respectively. Combines activity levels across channels to generate AUROC
%values, plots against trial number, combining sessions into single huge plot.
%Writes TAS values (AUROC values and correctness of responses) to mat file:
%session number, time epochs, AUROC & correctness, for each condition.
%
%Prior to calling this function, run bj_SE_batch_write_trial_resp (batch
%file which calls bj_SE_batch_write_resp).
plotActFigs=0;
if nargin<1||isempty(animals)
    animals=[{'blanco'} {'jack'}];
end
if nargin<2||isempty(areas)
    areas=[{'v4_1'} {'v4_2'} {'v1_1'} {'v1_2'}];
    areas=[{'v4_1'} {'v1_1'} {'v1_2_1'} {'v1_2_2'} {'v1_2_3'}];
    areas=[{'v4_1'} {'v1_1'}];
end
test_epochs={0 512 512*2 512*3};%{[0 40 300] 529 [529*2 529*2+40 529*2+300] 529*3}
notEqual=[];
adjustSessions=[311 318 333 352];
adjustSessions=[333 352 398 451];
if plotActFigs==1
    for animalInd=1:length(animals)
        animal=animals{animalInd};
        for areaInd=1:length(areas)
            area=areas{areaInd};
            channels = main_channels(animal,area);
            sessionNums = main_raw_sessions_final(animal,area,[],0);
            [sampleContrasts testContrasts]=area_metadata(area);
            for epoch=4:4%only examine test presentation period
                if epoch==1
                    periods=[-durSpon 0];
                else
                    periods=[test_epochs{epoch-1} test_epochs{epoch}(1)];
                end
                for subPeriod=1:length(periods)-1
                    startTime=periods(subPeriod);
                    endTime=periods(subPeriod+1);
                    for sampleContrastsInd=1:length(sampleContrasts)
                        sampleContrast=sampleContrasts(sampleContrastsInd);
                        testContrast=testContrasts(sampleContrastsInd,:);
                        allSessActS=cell(length(testContrast),1);
                        allSessActT=cell(length(testContrast),1);
                        TASmat=cell(length(testContrast),1);
                        nSess=zeros(1,size(testContrasts,2));
                        allTrialsCount=zeros(1,size(testContrasts,2))+1;%keep track of number of correct plus incorrect trials during plotting of diff
                        allTrialsSTCount=zeros(1,size(testContrasts,2))+1;%keep track of number of correct plus incorrect trials during plotting of sample & test act
                        for i=1:length(sessionNums)
                            if sessionNums(i)~=355&&sessionNums(i)~=405&&sessionNums(i)~=435
                                allChActS=cell(length(testContrast),1);
                                allChActT=cell(length(testContrast),1);
                                respName=[area,'_',num2str(sessionNums(i)),'_',num2str(sampleContrast)];
                                respDataFolder=fullfile('F:','PL','TAS',animal);
                                respDataFolder=fullfile(respDataFolder,respName);
                                loadBehavText=['load ',respDataFolder,'.mat respArray'];
                                eval(loadBehavText);%read list of incorrect and correct trials
                                notEqual=0;
                                for h=1:length(channels)
                                    if strncmp(area,'v1_2',4)
                                        loadText=['load F:\PL\sample_test_activity\',animal,'_v1_2\ch',num2str(channels(h)),'_',num2str(sessionNums(i)),'_example_sample_test_act.mat epoch2 epoch4'];
                                    else
                                        loadText=['load F:\PL\sample_test_activity\',animal,'_',area,'\ch',num2str(channels(h)),'_',num2str(sessionNums(i)),'_example_sample_test_act.mat epoch2 epoch4'];
                                    end
                                    eval(loadText);%load channel activity
                                    for condInd=1:length(testContrast)
                                        if size(epoch2{condInd},2)~=size(epoch4{condInd},2)||size(epoch2{condInd},2)~=size(respArray{condInd},1)
                                            notEqual=1;
                                        else
                                            if isempty(allChActS{condInd})
                                                allChActS{condInd}=epoch2{condInd,1};
                                            else
                                                allChActS{condInd}=allChActS{condInd}+epoch2{condInd,1};%one column per trial, calculate sum over all channels
                                            end
                                            if isempty(allChActT{condInd})
                                                allChActT{condInd}=epoch4{condInd,1};
                                            else
                                                allChActT{condInd}=allChActT{condInd}+epoch4{condInd,1};
                                            end
                                        end
                                    end
                                end
                                if notEqual==0%if number of trials matches
                                    for condInd=1:length(testContrast)%find mean across channels
                                        allChActS{condInd}=allChActS{condInd}/length(channels);
                                        allChActT{condInd}=allChActT{condInd}/length(channels);
                                        diffSampTest{condInd}=allChActT{condInd}-allChActS{condInd};%find difference in activity- if positive, indicates higher response to test than to sample
                                        if i==1%generate one figure for each condition when processing the first session
                                            figTASconds(condInd)=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.5]); %
                                            set(figTASconds(condInd), 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                                        else
                                            figure(figTASconds(condInd));
                                        end
                                        for n=1:size(respArray{condInd},1)%for each trial, plot difference between sample & test act
                                            if respArray{condInd}(n)==1
                                                plot(allTrialsCount(condInd),diffSampTest{condInd}(n),'k.');hold on%correct response
                                            elseif respArray{condInd}(n)==-1
                                                plot(allTrialsCount(condInd),diffSampTest{condInd}(n),'r.');hold on%incorrect response
                                            end
                                            allTrialsCount(condInd)=allTrialsCount(condInd)+1;%keep track of trial number on x-axis (across all sessions)
                                        end
                                        nSess(i,condInd)=allTrialsCount(condInd)-0.5;%list of endpoints for each session
                                        allSessActS{condInd}=[allSessActS{condInd} allChActS{condInd}];
                                        allSessActT{condInd}=[allSessActT{condInd} allChActT{condInd}];
                                        TASmat{condInd}=[TASmat{condInd} diffSampTest{condInd}];
                                        if i==1%generate one figure for each condition when processing the first session
                                            figSTconds(condInd)=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.5]); %
                                            set(figSTconds(condInd), 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                                        else
                                            figure(figSTconds(condInd));
                                        end
                                        for n=1:size(respArray{condInd},1)%for each trial, plot absolute activity levels for sample & test
                                            if respArray{condInd}(n)==1
                                                plot(allTrialsSTCount(condInd),allChActS{condInd}(n),'k.');hold on%sample, correct response
                                                plot(allTrialsSTCount(condInd),allChActT{condInd}(n),'kx');hold on%test, correct response
                                            elseif respArray{condInd}(n)==-1
                                                plot(allTrialsSTCount(condInd),allChActS{condInd}(n),'r.');hold on%sample, incorrect response
                                                plot(allTrialsSTCount(condInd),allChActT{condInd}(n),'rx');hold on%test, incorrect response
                                            end
                                            allTrialsSTCount(condInd)=allTrialsSTCount(condInd)+1;%keep track of trial number on x-axis (across all sessions)
                                        end
                                    end
                                end
                            end
                        end
                        %print figure for each condition
                        startEndTime=['_',num2str(periods(subPeriod)),'_to_',num2str(periods(subPeriod+1))];
                        for condInd=1:length(testContrast)
                            figure(figTASconds(condInd));
                            ylimVals=get(gca,'Ylim');
                            xlimVals=get(gca,'Xlim');
                            line([nSess(:,condInd) nSess(:,condInd)],[ylimVals(1) ylimVals(2)],'Color','k','LineStyle',':');
                            line([xlimVals(1) xlimVals(2)],[0 0],'Color','k','LineStyle',':');
                            xlim([0 nSess(end,condInd)]);
                            TASName=['TAS_diff_cond_',num2str(condInd),'_',num2str(sampleContrast),startEndTime,'.mat'];
                            TASFolder=fullfile('F:','PL','TAS',animal,area);
                            if ~exist(TASFolder,'dir')
                                mkdir(TASFolder);
                            end
                            TASPath=fullfile('F:','PL','TAS',animal,area,TASName);
                            printtext=sprintf('print -dpng %s.png',TASPath);
                            set(gcf,'PaperPositionMode','auto')
                            eval(printtext);
                            figure(figSTconds(condInd));
                            ylimVals=get(gca,'Ylim');
                            xlimVals=get(gca,'Xlim');
                            line([nSess(:,condInd) nSess(:,condInd)],[ylimVals(1) ylimVals(2)],'Color','k','LineStyle',':');
                            line([xlimVals(1) xlimVals(2)],[0 0],'Color','k','LineStyle',':');
                            ylim(ylimVals);
                            xlim([0 nSess(end,condInd)]);
                            TASName2=['TAS_act_cond_',num2str(condInd),'_',num2str(sampleContrast),startEndTime,'.mat'];
                            TASPath2=fullfile('F:','PL','TAS',animal,area,TASName2);
                            printtext2=sprintf('print -dpng %s.png',TASPath2);
                            set(gcf,'PaperPositionMode','auto')
                            eval(printtext2);
                        end
                        TASmatName=['TAS_',num2str(sampleContrast),startEndTime,'.mat'];
                        TASmatFolder=fullfile('F:','PL','TAS',animal,area);
                        if ~exist(TASmatFolder,'dir')
                            mkdir(TASmatFolder);
                        end
                        TASmatPath=fullfile('F:','PL','TAS',animal,area,TASmatName);
                        if exist(TASmatPath,'file')
                            loadText=['load ',TASmatPath,' TASmat'];
                            eval(loadText);
                            for condInd=1:length(testContrast)
                                %                             TASmat{condInd}=[TASmat{condInd};diffSampTest{condInd}];%combine across all sessions
                                %                             allSessActS{condInd}=[allSessActS{condInd};allChActS{condInd}];
                                %                             allSessActT{condInd}=[allSessActT{condInd};allChActT{condInd}];
                            end
                        end
                        saveText=['save ',TASmatPath,' TASmat allSessActS allSessActT'];
                        eval(saveText);
                    end
                end
            end
        end
    end
end

%calculate AUROC values in bins, across all trials for all sessions:
separateImages=0;
binWidths=[50 150 200 250];
for binWidthInd=1:length(binWidths)
    binWidth=binWidths(binWidthInd);%number of trials per bin, to calculate AUROC values
    for animalInd=1:length(animals)
        animal=animals{animalInd};
        for areaInd=1:length(areas)
            area=areas{areaInd};
            [sampleContrasts testContrasts]=area_metadata(area);
            for epoch=4:4%only examine test presentation period
                if epoch==1
                    periods=[-durSpon 0];
                else
                    periods=[test_epochs{epoch-1} test_epochs{epoch}(1)];
                end
                for subPeriod=1:length(periods)-1
                    startTime=periods(subPeriod);
                    endTime=periods(subPeriod+1);
                    startEndTime=['_',num2str(periods(subPeriod)),'_to_',num2str(periods(subPeriod+1))];
                    for sampleContrastsInd=1:length(sampleContrasts)
                        sampleContrast=sampleContrasts(sampleContrastsInd);
                        testContrast=testContrasts(sampleContrastsInd,:);
                        colmapText=colormap(jet(size(testContrast,2)));
                        colmapText=[colmapText(1,:);132/255 22/255 216/255;202/255 65/255 223/255;colmapText(3:7,:);157/255 212/255 61/255;colmapText(10:12,:);178/255 111/255 12/255;colmapText(end,:)];
                        markerTexts='+x';
                        markerText=markerTexts(2);markerS=8;
                        if separateImages==1
                            figAUROCconds=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.5]); %
                            set(figAUROCconds, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                        else
                            if animalInd==1&&areaInd==1
                                figAUROCconds=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.5, 0.5]); %
                                set(figAUROCconds, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                            end
                            subplot(2,2,animalInd+2*(areaInd-1));
                        end
                        TASmatName=['TAS_',num2str(sampleContrast),startEndTime,'.mat'];
                        TASmatFolder=fullfile('F:','PL','TAS',animal,area);
                        if ~exist(TASmatFolder,'dir')
                            mkdir(TASmatFolder);
                        end
                        TASmatPath=fullfile('F:','PL','TAS',animal,area,TASmatName);
                        if exist(TASmatPath,'file')
                            loadText=['load ',TASmatPath,' TASmat allSessActS allSessActT'];
                            eval(loadText);
                            AUROC=cell(length(testContrast),1);
                            for condInd=1:length(testContrast)
                                for binStart=1:length(allSessActS{condInd})-binWidth+1;%lower edge of each bin
                                    higher=0;%count proportion of trials with higher test act than sample act and vice versa
                                    lower=0;
                                    for trialCount=binStart:binStart+binWidth-1
                                        if allSessActS{condInd}(trialCount)<allSessActT{condInd}(trialCount)
                                            higher=higher+1;
                                        elseif allSessActS{condInd}(trialCount)>allSessActT{condInd}(trialCount)
                                            lower=lower+1;
                                        end
                                    end
                                    AUROC{condInd}=[AUROC{condInd} higher/(higher+lower)];%calculate ROC value within each bin of width binWidth and slide window along by 1 trial at a time
                                end
                                realTrialNums=1:length(allSessActS{condInd})-binWidth+1;
                                translatedTrialNums=realTrialNums./realTrialNums(end);
                                plot(translatedTrialNums,AUROC{condInd},'Color',colmapText(condInd,:));
                                hold on
                            end
                            if separateImages==1||separateImages==0&&animalInd==2
                                xlim([0 1]);
                                ylim([0 1]);
                                for condInd=1:length(testContrast)
                                    yLimVals=get(gca,'ylim');
                                    xLimVals=get(gca,'xlim');
                                    unitSpace=(yLimVals(2)-yLimVals(1))/30;
                                    text('Position',[xLimVals(2)+(xLimVals(2)-xLimVals(1))/25 yLimVals(1)+unitSpace*condInd*2],'FontSize',9,'String',[markerText,'  ',num2str(testContrast(condInd)),'%'],'Color',colmapText(condInd,:));
                                end                                
                            end
                            if separateImages==1||separateImages==0&&animalInd+2*(areaInd-1)==1
                                xlabel('proportion of total number of trials per condition')
                                ylabel('AUROC value')
                            end
                            AUROCName=['TAS_AUROC_',num2str(sampleContrast),startEndTime,'_binwidth',num2str(binWidth)];
                            if separateImages==1
                                AUROCPath=fullfile('F:','PL','TAS',animal,area,AUROCName);
                            elseif separateImages==0&&animalInd+2*(areaInd-1)==4
                                AUROCPath=fullfile('F:','PL','TAS',AUROCName);
                            end
                            if separateImages==1||separateImages==0&&animalInd+2*(areaInd-1)==4
                                printtext3=sprintf('print -dpng %s.png',AUROCPath);
                                set(gcf,'PaperPositionMode','auto')
                                eval(printtext3);
                            end
                        end
                    end
                end
            end
        end
    end
end