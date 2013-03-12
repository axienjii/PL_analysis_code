function batch_identify_trialcorr_artifact_example_figures

animal='jack';area='v1_1';
sessions=[73 75];
initialThreshDrawing=1;
if initialThreshDrawing==1
    plotFigs=1;
    allNumRemoveTrials=sessions';
    allRemoveTrialsStats=[];
    histoFig=figure;
    if strcmp(animal,'blanco')
        if strncmp(area,'v4',2)
            cutoffs=0.54:0.01:0.6;
            cutoffs=0.50:0.01:0.51;
            xlimTail=[0.4 1];
        elseif strncmp(area,'v1',2)
            cutoffs=0.54:0.01:0.59;
            cutoffs=0.6:0.01:0.65;
            xlimTail=[0.4 1];
        end
    elseif strcmp(animal,'jack')
        if strncmp(area,'v4',2)
            cutoffs=0.42:0.01:0.50;
            xlimTail=[0.3 0.9];
        elseif strncmp(area,'v1',2)
            cutoffs=0.65:0.01:0.73;
            cutoffs=0.67;
            xlimTail=[0.6 1];
        end
    end
    for j=1:length(cutoffs)
        cutoff=cutoffs(j);
        %     cutoffFig=figure;
        for i=1:length(sessions)
            output_fname = trial_corr_artifact_finder(animal,sessions(i));%csv files containing list of trials (identified by timestamp for 'NLX_TRIAL_START' encode) with movement artifacts, for each session
            loadText=['load ',output_fname,' removeTrialsTimestamps rlist'];
            eval(loadText)
            if plotFigs==1
                figure(histoFig);
                subplot(1,2,i);
                if j==1
                    hist(rlist(:,1),300);
                end
                title(num2str(sessions(i)));
                %     ylim([-100 500]);
                %     xlim([0.4 0.7]);
                %     grandMean=mean(rlist(:,1));
                %     grandstd=std(rlist(:,1));
                %     cutoff=grandMean+grandstd;
                ylimVals=get(gca,'YLim');
%                 line([cutoff cutoff],[0 ylimVals(2)],'Color','r');
            end
            tooHighInd=(rlist(:,1)>cutoff);
            removeTrialsTimestamps=[unique(rlist(tooHighInd,4)) unique(rlist(tooHighInd,5))];
            output_fname = trial_corr_artifact_finder(animal,sessions(i));%csv files containing list of trials (identified by timestamp for 'NLX_TRIAL_START' encode) with movement artifacts, for each session
            saveText=['save ',output_fname,' removeTrialsTimestamps rlist'];
            eval(saveText)
            numRemoveTrials(i)=size(removeTrialsTimestamps,1);
            
            %     testCutoffs=grandMean+2*grandstd:0.001:1;
            %     numTrialspercutoff=zeros(length(testCutoffs),1);
            %     for j=1:length(testCutoffs)
            %         testCutoff=testCutoffs(j);
            %         tooHighInd=(rlist(:,1)>testCutoff);
            %         removeTrialsTimestamps=[unique(rlist(tooHighInd,4)) unique(rlist(tooHighInd,5))];
            %         numTrialspercutoff(j,1)=sum(tooHighInd);
            %     end
            %     figure(cutoffFig);
            %     subplot(ceil(length(sessions)/5),5,i);
            %     plot(testCutoffs,numTrialspercutoff);
            %     title(num2str(sessions(i)));
        end
        
        for i=1:length(sessions)
            subplot(1,2,i);
            xlim([0 1]);
            %         ylim([-100 2500]);
        end
        for i=1:length(sessions)
            subplot(1,2,i);
            xlim([0.3 0.8]);
            ylim([-50 500]);
        end
        allNumRemoveTrials=[allNumRemoveTrials numRemoveTrials'];
    end
    totalNumTrials=zeros(length(sessions),1);
    for i=1:length(sessions)
        loadText=['load F:\PL\vals_perf\',animal,'\vals_',num2str(sessions(i)),'.mat'];
        eval(loadText)
        if strcmp(area,'v1_2')
        totalNumTrials(i,1)=size(vals{1},1)+size(vals{2},1)+size(vals{3},1);
        else
        totalNumTrials(i,1)=size(vals,1);
        end
    end
    allNumRemoveTrials
    allTotalNumTrials=[sessions' totalNumTrials]
    for cutoff=2:length(cutoffs)+1
        allRemoveTrialsStats=[allRemoveTrialsStats allNumRemoveTrials(:,cutoff)./totalNumTrials*100];
    end
    allRemoveTrialsStats=[allNumRemoveTrials allTotalNumTrials(:,2) allRemoveTrialsStats];
end

%plot lines corresponding to desired thresholds:
subplot(1,2,1);
line([0.43 0.43],[0 ylimVals(2)],'Color','r');
subplot(1,2,2);
line([0.47 0.47],[0 ylimVals(2)],'Color','r');

%set axis limits for zoomed-in version:
subplot(1,2,1);
xlim([0.3 0.8])
ylim([-50 500])
subplot(1,2,2);
xlim([0.3 0.8])
ylim([-50 500])

%set axis limits for zoomed-out version:
subplot(1,2,1);
xlim([0 0.8])
ylim([-50 1600])
subplot(1,2,2);
xlim([0 0.8])
ylim([-50 1900])