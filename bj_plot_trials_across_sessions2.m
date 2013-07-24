function bj_plot_trials_across_sessions2(animals,areas,roving)
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
readData=0;
separateImages=0;
plotAUROC30Fig=1;
continuousLine=0;
drawAUROCsliding=0;
calcStatsEL=1;
if nargin<1||isempty(animals)
    animals=[{'blanco'} {'jack'}];
end
if nargin<2||isempty(areas)
    if roving==0
        areas=[{'v4_1'} {'v1_1'}];
    elseif roving==1
        areas=[{'v1_2'}];
        areas=[{'v1_2_1'} {'v1_2_2'}];
    end
end
test_epochs={0 512 512*2 512*3};%{[0 40 300] 529 [529*2 529*2+40 529*2+300] 529*3}
notEqual=[];
adjustSessions=[311 318 333 352];
adjustSessions=[333 352 398 451];
if readData==1
    for animalInd=1:length(animals)
        animal=animals{animalInd};
        for areaInd=1:length(areas)
            area=areas{areaInd};
            channels = main_channels(animal,area);
            sessionNums = main_raw_sessions_final(animal,area,[],0);
            if find(sessionNums==451)
                sessionNums=sessionNums(find(sessionNums~=451));
            end
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
                        AUROC30=cell(length(testContrasts),1);
                        allAct30S=cell(2,length(testContrast),1);
                        allAct30T=cell(2,length(testContrast),1);
                        sessTrials30=cell(1,length(testContrast),1);
                        for i=1:length(sessionNums)
%                             if sessionNums(i)~=405&&sessionNums(i)~=435
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
                                        loadText=['load F:\PL\sample_test_activity\',animal,'_',area,'\ch',num2str(channels(h)),'_',num2str(sessionNums(i)),'_',num2str(sampleContrast),'_example_sample_test_act.mat epoch2 epoch4'];
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
                                        if plotActFigs==1
                                            if i==1%generate one figure for each condition when processing the first session
                                                figTASconds(condInd)=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.5]); %
                                                set(figTASconds(condInd), 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                                            else
                                                figure(figTASconds(condInd));
                                            end
                                        end
                                        for n=1:size(respArray{condInd},1)%for each trial, plot difference between sample & test act
                                            if plotActFigs==1
                                                if respArray{condInd}(n)==1
                                                    plot(allTrialsCount(condInd),diffSampTest{condInd}(n),'k.');hold on%correct response
                                                elseif respArray{condInd}(n)==-1
                                                    plot(allTrialsCount(condInd),diffSampTest{condInd}(n),'r.');hold on%incorrect response
                                                end
                                            end
                                            allTrialsCount(condInd)=allTrialsCount(condInd)+1;%keep track of trial number on x-axis (across all sessions)
                                        end
                                        nSess(i,condInd)=allTrialsCount(condInd)-0.5;%list of endpoints for each session
                                        allSessActS{condInd}=[allSessActS{condInd} allChActS{condInd}];
                                        allSessActT{condInd}=[allSessActT{condInd} allChActT{condInd}];
                                        
                                        %calculate AUROC values for each
                                        %session using first & last 30% of
                                        %trials:
                                        first30=1:floor(size(respArray{condInd},1)*0.3);%30% of trials per session
                                        last30=size(respArray{condInd},1)-floor(size(respArray{condInd},1)*0.3)+1:size(respArray{condInd},1);%last 30% of trials
                                        totalTrials(1)=sum(allChActS{condInd}(first30)<allChActT{condInd}(first30))+sum(allChActS{condInd}(first30)>allChActT{condInd}(first30));
                                        AUROC30{condInd}(1,i)=sum(allChActS{condInd}(first30)<allChActT{condInd}(first30))/totalTrials(1);%proportion of first 30% of trials with higher test act
                                        totalTrials(2)=sum(allChActS{condInd}(last30)<allChActT{condInd}(last30))+sum(allChActS{condInd}(last30)>allChActT{condInd}(last30));
                                        AUROC30{condInd}(2,i)=sum(allChActS{condInd}(last30)<allChActT{condInd}(last30))/totalTrials(2);%proportion of first 30% of trials with higher test act
                                        allAct30S{1,condInd}=[allAct30S{1,condInd} allChActS{condInd}(first30)];%array containing sampleactivity for each trial within the first 30% of trials, across conditions and sessions
                                        allAct30S{2,condInd}=[allAct30S{2,condInd} allChActS{condInd}(last30)];
                                        sessTrials30{1,condInd}=[sessTrials30{1,condInd} zeros(1,floor(size(respArray{condInd},1)*0.3))+condInd];%array containing session number, one value for each trial, across conditions and sessions
                                        allAct30T{1,condInd}=[allAct30T{1,condInd} allChActT{condInd}(first30)];%array containing activity for each trial within the first 30% of trials, across conditions and sessions
                                        allAct30T{2,condInd}=[allAct30T{2,condInd} allChActT{condInd}(last30)];
                                        
                                        TASmat{condInd}=[TASmat{condInd} diffSampTest{condInd}];
                                        if plotActFigs==1
                                            if i==1%generate one figure for each condition when processing the first session
                                                figSTconds(condInd)=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.5]); %
                                                set(figSTconds(condInd), 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                                            else
                                                figure(figSTconds(condInd));
                                            end
                                        end
                                        for n=1:size(respArray{condInd},1)%for each trial, plot absolute activity levels for sample & test
                                            if plotActFigs==1
                                                if respArray{condInd}(n)==1
                                                    plot(allTrialsSTCount(condInd),allChActS{condInd}(n),'k.');hold on%sample, correct response
                                                    plot(allTrialsSTCount(condInd),allChActT{condInd}(n),'kx');hold on%test, correct response
                                                elseif respArray{condInd}(n)==-1
                                                    plot(allTrialsSTCount(condInd),allChActS{condInd}(n),'r.');hold on%sample, incorrect response
                                                    plot(allTrialsSTCount(condInd),allChActT{condInd}(n),'rx');hold on%test, incorrect response
                                                end
                                            end
                                            allTrialsSTCount(condInd)=allTrialsSTCount(condInd)+1;%keep track of trial number on x-axis (across all sessions)
                                        end
                                    end
                                end
%                             end
                        end
                        %print figure for each condition
                        if plotActFigs==1
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
                        end
                        %print figure for each condition
                        if plotAUROC30Fig==1
                            colmapText=colormap(jet(size(testContrast,2)));
                            colmapText=[colmapText(1,:);132/255 22/255 216/255;202/255 65/255 223/255;colmapText(3:7,:);157/255 212/255 61/255;colmapText(10:12,:);178/255 111/255 12/255;colmapText(end,:)];
                            markerTexts='+x';
                            markerText=markerTexts(2);markerS=8;
                            startEndTime=['_',num2str(periods(subPeriod)),'_to_',num2str(periods(subPeriod+1))];
                            if separateImages==1
                                figAUROC30conds=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.5]); %
                                set(figAUROC30conds, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                            else
                                if animalInd==1&&areaInd==1
                                    figAUROC30conds=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.5, 0.5]); %
                                    set(figAUROC30conds, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                                end
                                subplot(2,2,animalInd+2*(areaInd-1));
                            end
                            for condInd=1:length(testContrast)
                                if continuousLine==0
                                for sessionInd=1:length(sessionNums)
                                    plot([sessionInd sessionInd+0.4],[AUROC30{condInd}(1,sessionInd) AUROC30{condInd}(2,sessionInd)],'Marker','o','Color',colmapText(condInd,:));hold on%plot 2 datas points, for early and late trials
                                    plot([1:length(sessionNums)]+0.4,AUROC30{condInd}(2,:),'Marker','o','MarkerFaceColor',colmapText(condInd,:),'Color',colmapText(condInd,:),'LineStyle','none');%fill in markers for late trials
                                end
                                elseif continuousLine==1
%                                 plot(1:length(sessionNums),AUROC30{condInd}(1,:),'Marker','o','Color',colmapText(condInd,:));hold on
%                                 plot(1:length(sessionNums),AUROC30{condInd}(2,:),'Marker','o','MarkerFaceColor',colmapText(condInd,:),'Color',colmapText(condInd,:));
                                end
                            end
                            if separateImages==1||separateImages==0&&animalInd==2
                                for condInd=1:length(testContrast)
                                    yLimVals=get(gca,'ylim');
                                    xLimVals=get(gca,'xlim');
                                    unitSpace=(yLimVals(2)-yLimVals(1))/30;
                                    text('Position',[xLimVals(2)+(xLimVals(2)-xLimVals(1))/25 yLimVals(1)+unitSpace*condInd*2],'FontSize',9,'String',[markerText,'  ',num2str(testContrast(condInd)),'%'],'Color',colmapText(condInd,:));
                                end                                
                            end
                            if separateImages==1||separateImages==0&&animalInd+2*(areaInd-1)==1
                                xlabel('session number')
                                ylabel('AUROC value for 30% of trials')
                            end
                            if continuousLine==1
                                AUROC30Name=['TAS_AUROC30_',num2str(sampleContrast),startEndTime];
                            elseif continuousLine==0
                                AUROC30Name=['TAS_AUROC30_earlylate_',num2str(sampleContrast),startEndTime];
                            end
                            if separateImages==1
                                AUROC30Path=fullfile('F:','PL','TAS',animal,area,AUROC30Name);
                            elseif separateImages==0&&animalInd+2*(areaInd-1)==4
                                AUROC30Path=fullfile('F:','PL','TAS',AUROC30Name);
                            end
                            if separateImages==1||separateImages==0&&animalInd+2*(areaInd-1)==4
                                printtext4=sprintf('print -dpng %s.png',AUROC30Path);
                                set(gcf,'PaperPositionMode','auto')
                                eval(printtext4);
                            end
                        end
                        TASmatName=['TAS_',num2str(sampleContrast),startEndTime,'.mat'];
                        TASmatFolder=fullfile('F:','PL','TAS',animal,area);
                        if ~exist(TASmatFolder,'dir')
                            mkdir(TASmatFolder);
                        end
                        TASmatPath=fullfile('F:','PL','TAS',animal,area,TASmatName);
                        if exist(TASmatPath,'file')
                            loadText=['load ',TASmatPath];
                            eval(loadText);
                            for condInd=1:length(testContrast)
                                %                             TASmat{condInd}=[TASmat{condInd};diffSampTest{condInd}];%combine across all sessions
                                %                             allSessActS{condInd}=[allSessActS{condInd};allChActS{condInd}];
                                %                             allSessActT{condInd}=[allSessActT{condInd};allChActT{condInd}];
                            end
                        end
                        saveText=['save ',TASmatPath,' TASmat allSessActS allSessActT nSess AUROC30 allAct30S allAct30T sessTrials30'];
                        eval(saveText);
                    end
                end
            end
        end
    end
end

if drawAUROCsliding==1
    %calculate AUROC values in bins, across all trials for all sessions:
    separateImages=0;
    if roving==1
        binWidths=[50 100 150 200;25 50 75 100];
    else
        binWidths=[50 100 150 200 250;50 100 150 200 250];
    end
    for binWidthInd=1:size(binWidths,2)
        for animalInd=1:length(animals)
            animal=animals{animalInd};
            binWidth=binWidths(animalInd,binWidthInd);%number of trials per bin, to calculate AUROC values
            for areaInd=1:length(areas)
                area=areas{areaInd};
                if strcmp(area,'v1_2_1')
                    binWidths=[50 100 150 200;25 50 75 100];
                elseif strcmp(area,'v1_2_2')
                    binWidths=[50 100 150 200;50 100 150 200];
                end
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
                            if strncmp(area,'v1_2',4)
                                if separateImages==1
                                    figAUROCconds=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.5]); %
                                    set(figAUROCconds, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                                else
                                    if animalInd==1&&sampleContrastsInd==1
                                        figAUROCconds(areaInd+size(binWidths,2)*(binWidthInd-1))=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.5, 0.5]); %
                                        set(figAUROCconds(areaInd+size(binWidths,2)*(binWidthInd-1)), 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                                    else
                                        figure(figAUROCconds(areaInd+size(binWidths,2)*(binWidthInd-1)));
                                    end
                                    subplot(3,2,animalInd+2*(sampleContrastsInd-1));
                                end
                            else
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
                                    if binWidth+1<length(allSessActS{condInd})
                                        realTrialNums=1:length(allSessActS{condInd})-binWidth+1;
                                        translatedTrialNums=realTrialNums./realTrialNums(end);
                                        plot(translatedTrialNums,AUROC{condInd},'Color',colmapText(condInd,:));
                                        hold on
                                    else
                                        checkThis=1;
                                    end
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
                                if strncmp(area,'v1_2',4)
                                    if separateImages==1||separateImages==0&&animalInd+2*(sampleContrastsInd-1)==1
                                        xlabel('proportion of total number of trials per condition')
                                        ylabel('AUROC value')
                                    end
                                else
                                    if separateImages==1||separateImages==0&&animalInd+2*(areaInd-1)==1
                                        xlabel('proportion of total number of trials per condition')
                                        ylabel('AUROC value')
                                    end
                                end
                                if strncmp(area,'v1_2',4)
                                    if animalInd==1&&sampleContrastsInd==1
                                        title('Monkey 1');
                                    elseif animalInd==2&&sampleContrastsInd==1
                                        title('Monkey 2');
                                    end
                                else
                                    if animalInd+2*(areaInd-1)==1
                                        title('Monkey 1');
                                    elseif animalInd+2*(areaInd-1)==2
                                        title('Monkey 2');
                                    end
                                end
                                if roving==0
                                    AUROCName=['TAS_AUROC_nonroving_',num2str(sampleContrast),startEndTime,'_binwidth',num2str(binWidth)];
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
                        if roving==1
                            figure(figAUROCconds(areaInd+size(binWidths,2)*(binWidthInd-1)));
                            AUROCName=['TAS_AUROC_roving_',area,'_',startEndTime,'_binwidth',num2str(binWidths(1,binWidthInd)),'_',num2str(binWidths(2,binWidthInd))];
                            if separateImages==1
                                AUROCPath=fullfile('F:','PL','TAS',animal,area,AUROCName);
                            elseif separateImages==0&&animalInd==2&&sampleContrastsInd==3
                                AUROCPath=fullfile('F:','PL','TAS',AUROCName);
                            end
                            if separateImages==1||separateImages==0&&animalInd==2&&sampleContrastsInd==3
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

statsTable=[];
%read early and late session AUROCs, calculate stats
if calcStatsEL==1
    for animalInd=1:length(animals)
        animal=animals{animalInd};
        for areaInd=1:length(areas)
            area=areas{areaInd};
            sessionNums = main_raw_sessions_final(animal,area,[],0);
            [sampleContrasts testContrasts]=area_metadata(area);
            for epoch=4:4%only examine test presentation period
                if epoch==1
                    periods=[-durSpon 0];
                else
                    periods=[test_epochs{epoch-1} test_epochs{epoch}(1)];
                end
                for subPeriod=1:length(periods)-1
                    startEndTime=['_',num2str(periods(subPeriod)),'_to_',num2str(periods(subPeriod+1))];
                    for sampleContrastsInd=1:length(sampleContrasts)
                        sampleContrast=sampleContrasts(sampleContrastsInd);
                        testContrast=testContrasts(sampleContrastsInd);
                        TASmatName=['TAS_',num2str(sampleContrast),startEndTime,'.mat'];
                        TASmatFolder=fullfile('F:','PL','TAS',animal,area);
                        if ~exist(TASmatFolder,'dir')
                            mkdir(TASmatFolder);
                        end
                        TASmatPath=fullfile('F:','PL','TAS',animal,area,TASmatName);
                        if exist(TASmatPath,'file')
                            loadText=['load ',TASmatPath];
                            eval(loadText);
                            for condInd=1:length(testContrast)
                                %                             TASmat{condInd}=[TASmat{condInd};diffSampTest{condInd}];%combine across all sessions
                                %                             allSessActS{condInd}=[allSessActS{condInd};allChActS{condInd}];
                                %                             allSessActT{condInd}=[allSessActT{condInd};allChActT{condInd}];
                            end
                        end
                        loadText=['load ',TASmatPath,' TASmat allSessActS allSessActT nSess AUROC30'];
                        eval(loadText);
                        allAct30ArrS=[];
                        allAct30ArrT=[];
                        sessTrials30Arr=[];
                        elArr=[];
                        condArr=[];
                        for condInd=1:length(testContrasts)
                            [h,p,ci,stats]=ttest(AUROC30{condInd}(1,:),AUROC30{condInd}(2,:));
                            if roving==0
                                h0{animalInd,areaInd}(condInd)=h;
                                p0{animalInd,areaInd}(condInd)=p;
                                cis0{animalInd,areaInd}{condInd}=ci;
                                stats0{animalInd,areaInd}{condInd}=stats;
                            elseif roving==1
                                h0{areaInd,sampleContrastsInd}(animalInd,condInd)=h;
                                p0{areaInd,sampleContrastsInd}(animalInd,condInd)=p;
                                cis0{areaInd,sampleContrastsInd}{animalInd,condInd}=ci;
%                                 stats0{areaInd,sampleContrastsInd}(animal
%                                 Ind,condInd)={stats};%figure out how to
%                                 store this, if analysis is necessary
                            end
                            
                            %sample activity within-session changes, combine acros conditions:
                            allAct30ArrS=[allAct30ArrS allAct30S{1,condInd} allAct30S{2,condInd}];
                            sessTrials30Arr=[sessTrials30Arr sessTrials30{1,condInd} sessTrials30{1,condInd}];
                            elArr=[elArr zeros(1,length(sessTrials30{1,condInd}))+1 zeros(1,length(sessTrials30{1,condInd}))+2];%early and late trials
                            condArr=[condArr zeros(1,length([allAct30S{1,condInd} allAct30S{2,condInd}]))+condInd];
                            %test activity within-session changes:
                            allAct30ArrT=[allAct30ArrT allAct30T{1,condInd} allAct30T{2,condInd}];
                        end
                        [p,f,stats]=anovan(allAct30ArrS,{elArr,sessTrials30Arr,condArr});
                        figure;
                        [comparison,means1,h,gnames]=multcompare(stats)%means contains mean in first column, SD in second. early period in first row, late period in second
                        if roving==0
                            p1{animalInd,areaInd}=p(1);
                            f1{animalInd,areaInd}=f{2,6};
                            df1{animalInd,areaInd}=f{5,3};
                            stats1{animalInd,areaInd}=stats;
                        elseif roving==1
                            p1{animalInd,sampleContrastsInd}=p(1);
                            f1{animalInd,sampleContrastsInd}=f{2,6};
                            df1{animalInd,sampleContrastsInd}=f{5,3};
                            stats1{animalInd,sampleContrastsInd}=stats;
                        end
                        [p,f,stats]=anovan(allAct30ArrT,{elArr,sessTrials30Arr,condArr});
                        figure;
                        [comparison,means2,h,gnames]=multcompare(stats)
                        if roving==0
                            p2{animalInd,areaInd}=p(1);
                            f2{animalInd,areaInd}=f{2,6};
                            df2{animalInd,areaInd}=f{5,3};
                            stats2{animalInd,areaInd}=stats;
                        elseif roving==1
                            p2{animalInd,sampleContrastsInd}=p(1);
                            f2{animalInd,sampleContrastsInd}=f{2,6};
                            df2{animalInd,sampleContrastsInd}=f{5,3};
                            stats2{animalInd,sampleContrastsInd}=stats;
                        end
                        if roving==0
                            statsTable=[statsTable;means1(1,:) means1(2,:) df1{animalInd,areaInd} f1{animalInd,areaInd} p1{animalInd,areaInd} means2(1,:) means2(2,:) df2{animalInd,areaInd} f2{animalInd,areaInd} p2{animalInd,areaInd}];
                        elseif roving==1
                            statsTable=[statsTable;means1(1,:) means1(2,:) df1{animalInd,sampleContrastsInd} f1{animalInd,sampleContrastsInd} p1{animalInd,sampleContrastsInd} means2(1,:) means2(2,:) df2{animalInd,sampleContrastsInd} f2{animalInd,sampleContrastsInd} p2{animalInd,sampleContrastsInd}];
                        end
                    end
                end
            end
        end
    end                        
end
statsTable%for non-roving: remember to swap second and third rows. original order: first row- blanco V4, second row- blanco V1, third row- jack V4, fourth row- jack V1