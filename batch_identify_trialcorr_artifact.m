% **** Batch file ****
function batch_identify_trialcorr_artifact(animal, area, sessions, icopy, ncopies, istart)

channels = main_channels(animal,area);
if nargin<3 || isempty(sessions)
    sessions = main_raw_sessions(animal,area);
end
if nargin<5
    icopy = 1;
end
if nargin<6
    ncopies = 1;
end
if nargin<7
    istart = 1;
end

verbose = 10;

generateInitialList=0;
if generateInitialList==1
    parfor iter = icopy:ncopies:length(sessions)
        
        if iter<istart
            continue;
        end
        
        session = sessions(iter);
        
        if verbose; fprintf('%d(%d): %d/%d Processing session %.1f\n',...
                icopy, ncopies, iter, length(sessions), session);
        end
        
        close all;
        
%         try
            output_artifact_trial_list = trial_corr_artifact_finder(animal, session);%csv files containing list of trials (identified by timestamp for 'NLX_TRIAL_START' encode) with movement artifacts, for each session
            if exist(output_artifact_trial_list, 'file'); continue; end
            identify_trial_corr_artifact(animal, area, session, channels, max(0,verbose-1))
%         catch ME
%             animal
%             session
%             
%             disp(ME);
%             for i=1:length(ME.stack);
%                 last_error_made = ME.stack(i);
%                 fprintf('In ==> %s at %d.\n', last_error_made.file, last_error_made.line);
%             end
%         end
        
        fclose all;
    end
end

initialThreshDrawing=0;
if initialThreshDrawing==1
    plotFigs=1;
    allNumRemoveTrials=sessions';
    allRemoveTrialsStats=[];
    histoFig=figure;
    if strcmp(animal,'blanco')
        if strcmp(area,'v4')||strcmp(area,'v4_1')||strncmp(area,'v4',2)
            cutoffs=0.44:0.01:0.49;
            cutoffs=0.64:0.001:0.65;
        elseif strcmp(area,'v4_2')
            cutoffs=0.45:0.01:0.8;
            cutoffs=0.45:0.01:0.51;
            xlimTail=[0.4 1];
        elseif strncmp(area,'v1',2)
            cutoffs=0.54:0.01:0.59;
            cutoffs=0.79:0.01:0.85;
            cutoffs=0.76:0.01:0.78;
            xlimTail=[0.4 1];
        end
    elseif strcmp(animal,'jack')
        if strncmp(area,'v4',2)||strncmp(area,'v4',2)
            cutoffs=0.42:0.01:0.50;
            cutoffs=0.51:0.01:0.64;
            cutoffs=0.43:0.005:0.45;
            xlimTail=[0.3 0.9];
        elseif strncmp(area,'v1',2)
            cutoffs=0.75:0.01:0.85;
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
                subplot(ceil(length(sessions)/5),5,i);
                if j==1
                    hist(rlist(:,1),300);
                end
                title(num2str(sessions(i)));
                %     ylim([-100 500]);
                %     xlim([0.4 0.7]);
                ylimVals=get(gca,'YLim');
                line([cutoff cutoff],[0 ylimVals(2)],'Color','r');
            end
%             grandMean=mean(rlist(:,1));
%             grandstd=std(rlist(:,1));
%             cutoff=grandMean+grandstd;
            tooHighInd=(rlist(:,1)>cutoff);
            removeTrialsTimestamps=[unique(rlist(tooHighInd,4)) unique(rlist(tooHighInd,5))];
            output_fname = trial_corr_artifact_finder(animal,sessions(i));%csv files containing list of trials (identified by timestamp for 'NLX_TRIAL_START' encode) with movement artifacts, for each session
            saveText=['save ',output_fname,' removeTrialsTimestamps rlist'];
            eval(saveText)
            numRemoveTrials(i)=size(removeTrialsTimestamps,1);
            
%                 testCutoffs=grandMean+2*grandstd:0.001:1;
%                 numTrialspercutoff=zeros(length(testCutoffs),1);
%                 for k=1:length(testCutoffs)
%                     testCutoff=testCutoffs(k);
%                     tooHighInd=(rlist(:,1)>testCutoff);
%                     removeTrialsTimestamps=[unique(rlist(tooHighInd,4)) unique(rlist(tooHighInd,5))];
%                     numTrialspercutoff(k,1)=sum(tooHighInd);
%                 end
            %     figure(cutoffFig);
            %     subplot(ceil(length(sessions)/5),5,i);
            %     plot(testCutoffs,numTrialspercutoff);
            %     title(num2str(sessions(i)));
        end
        
%         for i=1:length(sessions)
%             subplot(ceil(length(sessions)/5),5,i);
%             xlim([0 1]);
%             %         ylim([-100 2500]);
%         end
        for i=1:length(sessions)
            subplot(ceil(length(sessions)/5),5,i);
%             xlim(xlimTail);
            ylim([-100 500]);
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

%Once histograms have been generated and scrutinised, and threshold values
%manually selected and written to file 'sessThresholds.mat':
plotFigs=1;
allNumRemoveTrials=sessions';
histoFig=figure;
sessThresholdsPath=fullfile('F:','PL','pl_corr_art_trials',animal,'sessThresholds.mat');
loadText=['load ',sessThresholdsPath,' sessThresholds'];
eval(loadText);
if strcmp(animal,'blanco')
    if strncmp(area,'v4',2)
        xlimTail=[0.4 1];
    elseif strncmp(area,'v1',2)
        xlimTail=[0.4 1];
    end
elseif strcmp(animal,'jack')
    if strncmp(area,'v4',2)
        xlimTail=[0.3 0.9];
    elseif strncmp(area,'v1',2)
        xlimTail=[0.6 1];
    end
end
%     cutoffFig=figure;
for i=1:length(sessions)
    cutoff=sessThresholds(find(sessThresholds(:,1)==sessions(i)),2);
    output_fname = trial_corr_artifact_finder(animal,sessions(i));%csv files containing list of trials (identified by timestamp for 'NLX_TRIAL_START' encode) with movement artifacts, for each session
    loadText=['load ',output_fname,' removeTrialsTimestamps rlist'];
    eval(loadText)
    if plotFigs==1
        figure(histoFig);
        subplot(ceil(length(sessions)/5),5,i);
        hist(rlist(:,1),300);
        title(num2str(sessions(i)));
        %     ylim([-100 500]);
        %     xlim([0.4 0.7]);
        %     grandMean=mean(rlist(:,1));
        %     grandstd=std(rlist(:,1));
        %     cutoff=grandMean+grandstd;
        ylimVals=get(gca,'YLim');
        line([cutoff cutoff],[0 ylimVals(2)],'Color','r');
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
    subplot(ceil(length(sessions)/5),5,i);
    xlim([0 1]);
%             ylim([-100 2500]);
end
for i=1:length(sessions)
    subplot(ceil(length(sessions)/5),5,i);
    xlim(xlimTail);
    ylim([-100 500]);
end
allNumRemoveTrials=[allNumRemoveTrials numRemoveTrials'];
totalNumTrials=zeros(length(sessions),1);
for i=1:length(sessions)
    loadText=['load F:\PL\vals_perf\',animal,'\vals_',num2str(sessions(i)),'.mat'];
    eval(loadText)
    if strcmp(area,'v1_2')
        vals=[vals{1};vals{2};vals{3}];
    end
    totalNumTrials(i,1)=size(vals,1);
end
allRemoveTrialsStats=[allNumRemoveTrials totalNumTrials allNumRemoveTrials(:,2)./totalNumTrials*100]