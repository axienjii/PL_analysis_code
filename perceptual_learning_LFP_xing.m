function[]=perceptual_learning_LFP_xing(animal,area,channel, session, CELL,cycles,split_spike_timepoints, varargin)
%  drug_NLX_jones({'grcjdr1.'},[46],'grcjdr1.46 on_grcjdr1.46 offSc1.nse', 1,'V:\thielelab\Groups\ThieleGroup\monkey_data\Jones\Recordings\pen040\NLX_control\2012-10-23_15-22-24',[1],1, 'drugtext')
% NAME: 'grcjdr1.'
% exten: [46]
% file: sorted spike file
% CELL: cell number in sorted spikes file
% NLX_Path: full path (boring)
% cycles: determines which cycles are to be used for the analysis,
%           important for the analysis of e.g. different drug level
% split_spike_timepoints: if 1 it only uses trials from frist spike in
%   sorted spike file until last spike of interest in sorted spike file
%   use_spikes_for_start_end_time: if 0 use on_time and off_time as determined
%   by search for filename, otherwise use only events that occur after the
%   first sorted spike and before the last sorted spike
% drugtext: indicate clearly which drug was used (e.g. scopolamine...).




% Neuralynx encodes in columns 1 & 2:
% NLX_TRIAL_START and NLX_TRIAL_END, respectively.
%
% A copy of the full set of .mat files with timestamps for excluded trials is here:
% M:\Xing\pl_corr_art_trials_copy
%
% So for session 307 for Blanco, for example, the .mat file is M:\Xing\pl_corr_art_trials_copy\blanco\307_corrtrialartifact.mat
%
% ------------------------
%
% This is how I read trials from .nev files- the 'correct sacacde' or 'incorrect saccade' encode must be present; furthermore, the number of encodes recorded for a given trial has to be correct.
%
% Length of 39 for correct trials (encode '43'), and length of 37 for incorrect trials (encode '44').
% This is very important, as sometimes encodes get messed up/duplicated, and if our trial inclusion differs (e.g. you don't check that the length of encodes on a given trial is correct), then the trials won't match up:
%

%
% The final .mat spike files do not have timestamp IDs (though I could generate them if necessary)- I read the .nse files and make sure that only artefact-free trials are included, and the spike data is written to .mat files.
%
% A copy of the full set of .mat files with spike timestamps from good trials is here:
% M:\Xing\spikeData
%
% Each .mat file corresponds to a channel and a session.
%
% The spike timestamps are in the variable 'matarray' and each row corresponds to a condition (in ascending order), while each column corresponds to a time epoch (spontaneous; sample, ISI, test, post-test).
%
% Each cell contains rows of trials; for example, if condition 1 had 20 trials, then there would be 20 rows in each of the five cells in the first row of cells of 'matarray.'
%
% Within each trial is a list of timestamps- the rate is calculated by finding the number of timestamps that occur within a given epoch.





spontanSampOnly=0;%set to 1 to only include spontan and sample epochs. set to 0 to additionally include ISI and test epochs

%%%%%%%%%%%%%%%%% cell classification

close all hidden;

read_save_data=0;

if read_save_data==0
    %%%%%% flags defining what is to be saved for population
    %%%%%% analysis after initial analyis is over
    write_LFP_data=1;
    filter_50Hz=1;
    %filter_75Hz=0;
    
    params.tapers=[3 5];
    params.pad=0;
    params.Fs=1008;
    params.trialave=1;
    params.fpass=[0 500];
    
    %%%%%%%%%%%%% if block changes are to be tracked
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LFP_file=['M:\Xing\pl_LFP\',animal,'\',num2str(session),'\lfp_CSC',num2str(channel),'.ncs'];
    artifact_file=['M:\Xing\pl_corr_art_trials_copy\',animal,'\',num2str(session),'_corrtrialartifact.mat'];
    load (artifact_file)
    % subtract_spont=0; % analyse with or without spontaneous subtracted, not sure it is currently implements anywhere
    [first_file, test_contrasts, sample_contrasts, expt_type, rotated, areaDummy, rovingType,second_file] = session_metadata_alex(session, animal);
    numcond = length(test_contrasts);
    conditions=1:1:length(test_contrasts);
    numblock=1;
    
%     evaldir=sprintf('%s',NLX_Path); % this simply ensures that we are in the subdirectory where data are
%     eval(evaldir);
    
    if strcmp(animal,'blanco')
        eventsFilePath=['V:\thielelab\Groups\ThieleGroup\monkey_data\blanco\_grid\',area,'_events_files\',num2str(session),'\Events.Nev'];
    elseif strcmp(animal,'jack')
        eventsFilePath=['V:\thielelab\Groups\ThieleGroup\monkey_data\Jack\_jackgrid\j_',area,'_events_files\',num2str(session),'\Events.Nev'];
    end
    
    %%% get the timestamps that correspond to the period of interest
    disp(['First cortex file: ', first_file]);
    [on_time1, off_time1]=find_NLXFILE_start_xing(first_file,eventsFilePath);
    [on_time2,off_time2]=find_NLXFILE_start_xing(second_file,eventsFilePath);
    off_time2=off_time2; 
    
        
    ExtractMode = 4;
    ModeArray(1) = on_time1(1);
    ModeArray(2) = off_time2(end);     
    [TimeStamps_EV,  event_arr_NLX, EventStrings ] = Nlx2MatEV(eventsFilePath, [1 0 1 0 1],0, ExtractMode, ModeArray ); %#ok<NASGU> %
    
    %%% The encodes are buffered by zeros, ALWIN THINKS THAT IS MORE SECURE and he is probably correct,
    %%% and here we eliminate the zeros
%     NonZeroEvents = find(event_arr_NLX~=0);
%     event_arr_NLX = event_arr_NLX(NonZeroEvents);
%     TimeStamps_EV = TimeStamps_EV(NonZeroEvents);
    
    %%%%%%%% Read the header of the continuous file (can be commented out since it is not used currently)
    ExtractHeader = 1;
    
    
    [NlxHeader] = Nlx2MatCSC(LFP_file, [0 0 0 0 0 ], ExtractHeader,1 )
    %[Timestamp, ChanNum, SampleFrequency, NumValSamples, Samples, NlxHeader] = Nlx2MatCSC( LFP_file, [1 1 1 1 1 ], ExtractHeader,1,1 );
    
    %%%%%%%%%%%%% read the LFP in a format where each LFP data point is
    %%%%%%%%%%%%% associated with a time stamp, that is not the case in the
    %%%%%%%%%%%%% original format, only needed if LFP analysis or SFC  is done.
    [time_Arr_LFP, value_Arr, SampFreq]=read_nlx_CSC (LFP_file);
    
    prl=find(time_Arr_LFP>on_time1(1) & time_Arr_LFP<off_time2(end));
    time_Arr_LFP=time_Arr_LFP(prl);
    value_Arr=value_Arr(prl);
    for jk=1:size(NlxHeader,1)
        prl=strmatch('-ADBitVolts', NlxHeader{jk});
        if prl==1
            prl=find(NlxHeader{jk}==' ');
            prl2=NlxHeader{jk};
            prl2=str2double(prl2(prl+1:end))
            Header.ADBitVolts=prl2
            value_Arr=value_Arr*(-1)*Header.ADBitVolts;
        end
    end
    
    prl=ceil(time_Arr_LFP/1000);
    prl0=prl-prl(1);
    length(time_Arr_LFP)-length(prl)
    prl2=prl(2:end)-prl(1:end-1);
    prl3=find(prl2~=0);
    prl4=find(prl2==0);
    time_Arr_LFP=time_Arr_LFP(prl3);
    value_Arr=value_Arr(prl3);
    SampFreq=1000;
    
    
    
    params.Fs=SampFreq;
    
    
    if filter_50Hz==1;
        [b,a] = butter(3,[49.7  50.3]/(SampFreq/2),'bandpass');
        y2 = filtfilt(b,a,value_Arr);
        value_Arr=value_Arr+(y2*-1);
    end
    
    disp(['Number of blocks ',int2str(numblock)]);
    disp(['Number of conditions ',int2str(numcond)]);
    
    
    %     %%%% open a figure into which to plot
    %     figinfo=[NLX_Path ' ' name ' ' num2str(exten) ' Cell ' num2str(CELL) ' Cycles ' num2str(cycles)];
    %     pos=[0.01, 0.01, 0.85, 0.85];
    %     fighandle1=  figure('Color',[1,1,1],'Units', 'Normalized', 'Position',pos); %
    %     set(fighandle1, 'NumberTitle', 'off', 'Name', figinfo);
    %     set(fighandle1, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.0 0.0 29.305 20.65]);
    %
    %
    
    %%%%%%% trial counters, which keep track of how often a condition/block
    %%%%%%% combination actually occurred. Important info, for e.g. average
    %%%%%%% firing rates for a condition etc.
    ndir=zeros(1,numcond);
    %%%%%%%%%%%% variables to contain the spiking information
    
    
    
    %     curr_file=sprintf('%s%d',name,exten(1)); %current file of interest in case there are multiple
    %     [on_time, off_time]=find_NLXFILE_start(curr_file);  % this finds the time stamps in neuralynx, which correspond to those data files, necessary to match spikes and LFP time stamps
    disp(['on_time ',int2str(on_time1)]);
    disp(['off_time ',int2str(off_time2)]);
    curr_EV=event_arr_NLX(TimeStamps_EV>=on_time1 & TimeStamps_EV<=off_time2); % now selecting the events
    curr_TimeEV=TimeStamps_EV(TimeStamps_EV>=on_time1 & TimeStamps_EV<=off_time2); %time staps of events
    
    
    
    NLX_Trial_PARAM_START=find(curr_EV==253); % get a representation (starting point) of all the different trials that occurred
    
    %%%%%%% get a counter for where to start in the neuralynxfile and
    %%%%%%% cortexfile.
    trialcounts=[1 length(NLX_Trial_PARAM_START)];
    
    PARAMBASE =10000;
    TRIALPARAM_START =300;
    TRIALPARAM_END =301;
    STIMPARAM_START =302;
    STIMPARAM_END =303;
    STIM_SWITCH =304;
    REWARDPARAM_START =305;
    REWARDPARAM_END =306;
    STIM1_ON =307; %sample
    STIM1_OFF =308;
    STIM2_ON =309; %test
    STIM2_OFF =310;
    STIM3_ON =311; %targets
    FIX_COL_CHANGE =312;
    SACCADE_CORR =313;
    SACCADE_DIST =314;
    SACCADE_TEST =315;
    FIX_BREAK =316;
    REW_NUM_BASE =320;
    SACCADE_ONSET =5006;
    NLX_TRIAL_START =255;
    NLX_RECORD_START =2;
    NLX_SUBJECT_START =4;
    NLX_STIM_ON =8;
    NLX_STIM_OFF =16;
    NLX_TEST_ON =102;
    NLX_TEST_OFF =103;
    NLX_SUBJECT_END =32;
    NLX_RECORD_END =64;
    NLX_TRIAL_END =254;
    NLX_FIX_BREAK =39;
    NLX_TARGET_ON =40;
    NLX_TARGET_OFF =41;
    NLX_FIX_COL_CHANGE =42;
    NLX_SACCADE_CORR =43;
    NLX_SACCADE_ERROR =44;
    NLX_SACCADE_TEST =45;
    NLX_REW_NUM =46;
    NLX_FIX_SPOT_ON =59;
    NLX_FIX_START =60;
    NLX_TESTDIMMED =17;
    NLX_DISTDIMMED =18;
    NLX_BARRELEASED =19;
    NLX_CUE_ON =20;
    NLX_CUE_OFF =21;
    NLX_EVENT_1 =9;
    NLX_EVENT_2 =10;
    NLX_EVENT_3 =11;
    NLX_EVENT_4 =12;
    NLX_EVENT_5 =13;
    NLX_EVENT_6 =14;
    NLX_EVENT_7 =15;
    NLX_READ_DATA =128;
    %condition parameter encodes (send as 1 bytes)
    NLX_TRIALPARAM_START =253;
    NLX_TRIALPARAM_END =252;
    %Stimparameter encodes (send as 2 bytes)
    NLX_STIMPARAM_START =251;
    NLX_STIMPARAM_END =250;
    
    sgl_trial_LFP=zeros(100,8,1024);
    
    corr_resp=zeros(1,max(conditions));
    error_resp=zeros(1,max(conditions));
    prl=find(curr_EV==255);
    count_resp_trials=0;
    max_event_length=0;
    event_arr=zeros(1,max(conditions));
    time_arr=zeros(max(conditions),100);
    for k=1:length(prl)-1%look through all trials
        for j=1:numcond
            if curr_EV(prl(k)+6)==conditions(j)%condition number (Alex previously set to 4, as count started from NLX_TRIALPARAM_START instead of NLX_TRIAL_START
                encode_arr=curr_EV(prl(k):prl(k+1));%first value is 255 for start of current trial, ends with 255 of next trial
                corr_resp(1,j)=corr_resp(1,j)+length(find(encode_arr==43));%#define NLX_SACCADE_CORR
                error_resp(1,j)=error_resp(1,j)+length(find(encode_arr==44));%#define NLX_SACCADE_ERROR
                if ~isempty(find(encode_arr==43))&&length(encode_arr)==39||~isempty(find(encode_arr==44))&&length(encode_arr)==37
                    count_resp_trials=count_resp_trials+1;
                    if length(encode_arr)>max_event_length
                        max_event_length=length(encode_arr);
                    end
%                     if ~(length(encode_arr)==37||length(encode_arr)==39)
%                         length(encode_arr)
%                     end
                    event_arr(1,j)=event_arr(1,j)+1;%each trial occupies 1 column, trial number is value in first row
                    time_arr(j,event_arr(1,j))=curr_TimeEV(prl(k));
                end
            end
        end
    end
    LFP_data=zeros(max(conditions),max(event_arr),1025); %%% downsampled to 1 kHz
    corr_resp=zeros(1,max(conditions));
    error_resp=zeros(1,max(conditions));
    prl=find(curr_EV==255);
    count_resp_trials=0;
    max_event_length=0;
    event_arr=zeros(1,max(conditions));
    time_arr=zeros(max(conditions),100);
    skipTrialTime=[];
    for k=1:length(prl)-1%look through all trials
        for j=1:numcond
            if curr_EV(prl(k)+6)==conditions(j)%condition number (Alex previously set to 4, as count started from NLX_TRIALPARAM_START instead of NLX_TRIAL_START
                encode_arr=curr_EV(prl(k):prl(k+1));%first value is 255 for start of current trial, ends with 255 of next trial
                t_arr=curr_TimeEV(prl(k):prl(k+1));
                find_252=find(encode_arr==252);
                new_encode=encode_arr(find_252+1:end);
                new_time=t_arr(find_252+1:end);
                corr_resp(1,j)=corr_resp(1,j)+length(find(encode_arr==43));%#define NLX_SACCADE_CORR
                error_resp(1,j)=error_resp(1,j)+length(find(encode_arr==44));%#define NLX_SACCADE_ERROR
                if ~isempty(find(encode_arr==43))&&length(encode_arr)==39||~isempty(find(encode_arr==44))&&length(encode_arr)==37
                    count_resp_trials=count_resp_trials+1;
                    if length(encode_arr)>max_event_length
                        max_event_length=length(encode_arr);
                    end
                    %                     if ~(length(encode_arr)==37||length(encode_arr)==39)
                    %                         length(encode_arr)
                    %                     end
                    sample_on=new_time(new_encode==8);%spontaneous and sample act
                    if spontanSampOnly==1
                        LFP=value_Arr(time_Arr_LFP>sample_on-512*10^3 & time_Arr_LFP<sample_on+(1033-512)*10^3);
                    elseif spontanSampOnly==0
                        test_on=new_time(new_encode==102);%spontan, sample, ISI, and test act
                        LFP1=value_Arr(time_Arr_LFP>sample_on-512*10^3 & time_Arr_LFP<sample_on+(1024)*10^3);
                        LFP2=value_Arr(time_Arr_LFP>test_on & time_Arr_LFP<test_on+(1033-512)*10^3);
                        LFP=[LFP1;LFP2];
                    end
                    event_arr(1,j)=event_arr(1,j)+1;%each trial occupies 1 column, trial number is value in first row
                    time_arr(j,event_arr(1,j))=curr_TimeEV(prl(k));
                    %t_p=time_Arr_LFP(time_Arr_LFP>stim_on-512*10^3 & time_Arr_LFP<stim_on-502*10^3);
                    LFP_data(j,event_arr(1,j),1)=t_arr(1);%first 512 entries correspond to spontan act, second half correspond to post-stim onset
                    if spontanSampOnly==1
                        if ~isempty(LFP)&&length(LFP)>=1024
                            LFP_data(j,event_arr(1,j),2:1025)=LFP(1:1024);
                        else
                            skipTrialTime=[skipTrialTime;t_arr(1)];
                        end
                    elseif spontanSampOnly==0
                        if ~isempty(LFP)&&length(LFP)>=2048%512*4=2048
                            LFP_data(j,event_arr(1,j),2:2049)=LFP(1:2048);
                        else
                            skipTrialTime=[skipTrialTime;t_arr(1)];
                        end
                    end
                end
            end
        end
    end
    if ~isempty(skipTrialTime)
        skipTrialTime
        channel
        session
    end
    if spontanSampOnly==1%spontaneous and sample act
        saveText=['save M:\Xing\pl_LFP\',animal,'\',num2str(session),'\SGL_trial_LFP_ch',num2str(channel),'.mat LFP_data event_arr'];
    elseif spontanSampOnly==0%spontan, sample, ISI, and test act
        saveText=['save M:\Xing\pl_LFP\',animal,'\',num2str(session),'\SGL_trial_LFP_-512_1536_ch',num2str(channel),'.mat LFP_data event_arr'];
    end
    eval(saveText)
end
if spontanSampOnly==1%spontaneous and sample act
    loadText=['load M:\Xing\pl_LFP\',animal,'\',num2str(session),'\SGL_trial_LFP_ch',num2str(channel),'.mat LFP_data event_arr'];
elseif spontanSampOnly==0%spontan, sample, ISI, and test act
    loadText=['load M:\Xing\pl_LFP\',animal,'\',num2str(session),'\SGL_trial_LFP_-512_1536_ch',num2str(channel),'.mat LFP_data event_arr'];
end
eval(loadText)
for j=1:14
    mean_LFP=(squeeze(LFP_data(j,1:event_arr(1,j),2:1025)));
    %mean_LFP=mean(mean_LFP')
    test=subplot(1,1,1);
    for k=1:size(mean_LFP,1)
        plot(mean_LFP(k,:));
        hold on    
    end
    mean_LFP=mean(mean_LFP);
    plot(mean_LFP,'r');
    clf
end

