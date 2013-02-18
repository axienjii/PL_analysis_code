function identify_trial_corr_artifact(animal, area, session, channels, verbose)
%Written by Xing 11/2/13.
%Reads MUA files with non-rectified, downsampled MUA data. Checks for
%degree of correlation of MUA data between every possible pair of channels,
%for each trial (e.g. if there are 10 channels, calculates R for 45 pairs.)
%End up with number of possible pairs x number of trials of R values (e.g.
%for 10 channels and 1000 trials per session, have 45*1000 R values). Plot
%distribution of all R values and see if there exists a subpopulation of R
%values that are much higher than the rest (e.g. most R values are
%distributed around a mean of 0.4, but some are much higher at 0.6 and
%above). If so, identify trials that correspond to outlyingly high R
%values, and obtain time stamp of encode of 'NLX_TRIAL_START,' writes these
%timestamps to list and save in .csv file.

drawWaitbar=0;
%load vals from .mat file, to find pertinent time stamps for each trial
%from session of interest
loadText=['load V:\thielelab\Groups\ThieleGroup\monkey_data\blanco\_grid\vals_perf\',animal,'\vals_',num2str(session),'.mat'];
eval(loadText);
[file_of_int, test_contrasts, sample_contrasts, expt_type]=session_metadata(session,animal);
if expt_type==2
    vals=[vals{1};vals{2};vals{3}];
end
nev_fname=nev_finder(animal,area,session);
[on_time,off_time] = get_on_off_times(file_of_int, nev_fname);

%Create list of channel pairs
rowInd=1;
channelPairList=zeros((length(channels)*length(channels)-length(channels))/2,2);
for channel1=1:length(channels)
    for channel2=channel1+1:length(channels)
        channelPairList(rowInd,:)=[channels(channel1) channels(channel2)];
        rowInd=rowInd+1;
    end
end

rlist=zeros(size(channelPairList,1)*size(vals,1),5);
rlistCount=1;
removeRow=zeros(size(channelPairList,1)*size(vals,1),1);
removeRowCount=1;
%For each pair of channels, read MUA files with non-rectified, downsampled
%MUA data.
if drawWaitbar==1
    waitbarText=['calculating R values for session ',num2str(session)];
    waitbarHandle=waitbar(0,waitbarText);
end
for channel1=1:length(channels)
    % Extract data from the NCS file for the first channel in the pair
    if ispc
        input_fname_ch1=nonRectified_MUA_ncs_finder(animal,session,channels(channel1));
        [TimeStamps1 samplesCh1]=Nlx2MatCSC(input_fname_ch1,[1 0 0 0 1],0,4,[on_time off_time]);
    end
    for channel2=channel1+1:length(channels)
        input_fname_ch2=nonRectified_MUA_ncs_finder(animal,session,channels(channel2));
        
        % Extract data from the NCS file for the second channel in the pair
        if ispc
            if channels(channel2)==55
                x=1;
            end
            [TimeStamps2 samplesCh2]=Nlx2MatCSC(input_fname_ch2,[1 0 0 0 1],0,4,[on_time off_time]);
        end
        for trialInd=1:size(vals,1)
            if drawWaitbar==1
                waitbar(rlistCount/(size(channelPairList,1)*size(vals,1)));
            end
            trialTimeStamps1=TimeStamps1((TimeStamps1>=vals(trialInd,1))+(TimeStamps1<=vals(trialInd,21))==2); %1: NLX_TRIAL_START 21: NLX_TRIAL_END
            trialSamples1=samplesCh1(:,(TimeStamps1>=vals(trialInd,1))+(TimeStamps1<=vals(trialInd,21))==2);
            trialTimeStamps2=TimeStamps2((TimeStamps2>=vals(trialInd,1))+(TimeStamps2<=vals(trialInd,21))==2);
            trialSamples2=samplesCh2(:,(TimeStamps2>=vals(trialInd,1))+(TimeStamps2<=vals(trialInd,21))==2);
            %Check for degree of correlation of MUA data, for each trial. Note down
            %timestamp along with R value.
            if size(trialSamples1,1)~=size(trialSamples2,1)||size(trialSamples1,2)~=size(trialSamples2,2)||isempty(trialSamples1)||isempty(trialSamples2)
                removeRow(removeRowCount)=rlistCount;
                removeRowCount=removeRowCount+1;
            else
                r=corr(reshape(trialSamples1,size(trialSamples1,1)*size(trialSamples1,2),1),reshape(trialSamples2,size(trialSamples2,1)*size(trialSamples2,2),1));
                rlist(rlistCount,:)=[r channel1 channel2 vals(trialInd,1) vals(trialInd,21)];%R, first channel of the pair, second channel of the pair, timestamps corresponding to start and end of trial
            end
            rlistCount=rlistCount+1;
        end
    end
end
removeRow=removeRow(1:removeRowCount-1);
rlist(removeRow,:)=nan;
rlist=reshape(rlist(~isnan(rlist)),[],5);
if drawWaitbar==1
    close(waitbarHandle);
end

%Pool all R values together and plot their distribution to check whether
%there exists a subpopulation of R values that are much higher than the rest.
hist(rlist(:,1),300);
grandMean=mean(rlist(:,1));
grandstd=std(rlist(:,1));

%If so, set threshold for R value such that trials that correspond to
%outlyingly high R values above this value are identified as artifact-contaminated trials.
cutoff=grandMean+grandstd;
ylimVals=get(gca,'YLim');
line([cutoff cutoff],[0 ylimVals(2)],'Color','r');

%Collect unique time stamps of encode of 'NLX_TRIAL_START,' write these
%timestamps to list and save in .csv file.
tooHighInd=(rlist(:,1)>cutoff);
removeTrialsTimestamps=[unique(rlist(tooHighInd,4)) unique(rlist(tooHighInd,5))];
output_fname = trial_corr_artifact_finder(animal,session);%csv files containing list of trials (identified by timestamp for 'NLX_TRIAL_START' encode) with movement artifacts, for each session
saveText=['save ',output_fname,' removeTrialsTimestamps rlist'];
eval(saveText)
