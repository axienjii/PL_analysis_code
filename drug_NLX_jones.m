function[ATTrocs,ARMI, DrugROC, DrugMI,mean_act, Fanos, P, cell_text, drugtext, unit_identifier, corr_prctile, ndir]=drug_NLX_jones(NAME,exten,SE_file, CELL,NLX_Path,cycles,split_spike_timepoints, drugtext, unit_identifier)
%  drug_NLX_jones({'grcjdr1.'},[46],'grcjdr1.46 on_grcjdr1.46 offSc1.nse', 1,'V:\thielelab\Groups\ThieleGroup\monkey_data\Jones\Recordings\pen040\NLX_control\2012-10-23_15-22-24',[1],1, 'drugtext')
% NAME: 'grcjdr1.'
% exten: [46]
% file: sorted spike file
% CELL: cell number in sorted spikes file
% NLX_Path: full path (boring)
% cycles: determines which cycles are to be used for the analysis,
%           important for the analysis of e.g. different drug level  
% split_spike_timepoints: if 1 it only uses trials from frist spike in
%   sorted spike file until last spike of ispsnterest in sorted spike file
%   use_spikes_for_start_end_time: if 0 use on_time and off_time as determined
%   by search for filename, otherwise use only events that occur after the
%   first sorted spike and before the last sorted spike
% drugtext: indicate clearly which drug was used (e.g. scopolamine...).


close all hidden;


oldFolder = cd(NLX_Path); % change to the NLX directory since we are using the full path anyway

eye_win=3; %could be used to restrict analysis to good eye position trials, but not currently done
t_win=[100 400; 500 0];% first time window to be used for analysis of spikes after stimulus and cue onset, second for analysis of spikes before dimming

%%%%%% flags defining what is to be saved for population
%%%%%% analysis after initial analyis is over
write_histos=0;
write_data=0;
write_Tres_FF=0;
write_RT=0;
write_eps=0;
Read_ReactionTimes=0; %%% determines whether the high resolution bar file is loaded (if present), loading and analysing it is slow, so not needed when basic visual assessment of responses 

%%%%%%%%%%%%% if block changes are to be tracked
last_blk=-1; %%% a variable that tracks whether a block change has just occurred
elim_first_blkCHG_trials=6; %%% number of trials to be eliminated after a block change
blk_CHG_trial_counter=0; %%% tracks how many trials since last block change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SE_file=SE_file{:};

cutoff_act=3; %%%%  mean+ cutoff_act*standard deviation to eliminate trials where some activity went wrong

% subtract_spont=0; % analyse with or without spontaneous subtracted, not sure it is currently implements anywhere
normalize=0; % 1:=normalize activity, not sure implemented
plot_rasters=1; % if 1 it plots rasters, otherwise just histograms

numfiles=size(exten,2); %checks how many cortex files were used to do no-drug vs. drug condition, probably only 1 for pascoe and jone, but for pascoe we are currently not sure.

evaldir=sprintf('cd %s',NLX_Path); % this simple ensures that we are in the subdirectory where data are
eval(evaldir);
% prl=dir; % find out what is in this subdirectory
% prl=sprintf('cd %s',prl(3).name);
% eval(prl);

% %%% This block only shows the coordinates of the receptive field
% %%% it can be commented out without affecting any other code
% idx = max(strfind(NLX_Path,'\'));  % find the last directory seperator 
% Base_path = NLX_Path(1:idx-1);  % go one directory back in the path
% idx = max(strfind(Base_path,'\'));  % find the last directory seperator 
% Base_path = Base_path(1:idx);  % go one directory back in the path
% CTXFileName = [Base_path, 'cortex\'];
% Dir_Content = dir(CTXFileName);
% CTXFileName = [Base_path, 'cortex\',Dir_Content(3).name,'\',NAME{1},int2str(exten(1))];
% Coordinates = Read_RF_Coordinates(CTXFileName);
% 
% disp('Coordinates for the Receptive field and the out points');
% disp(['RF  = (',num2str(Coordinates.RF_X,2),',',num2str(Coordinates.RF_Y,2),')' ]);
% disp(['Out1= (',num2str(Coordinates.Out1_X,2),',',num2str(Coordinates.Out1_Y,2),')' ]);
% disp(['Out2= (',num2str(Coordinates.Out2_X,2),',',num2str(Coordinates.Out2_Y,2),')' ]);
% %%% end of coordinate block




%%% get the timestamps that correspond to the period of interest
name=NAME{1};
first_file=sprintf('%s%d',name,exten(1));
disp(['First cortex file: ', first_file]);
[on_time1, off_time2]=find_NLXFILE_start(first_file);
if size(NAME)>1 % in case the off encode was forgotten, enter a second name in function line which gives the encode of the next cortex file that was run. it will take the on time of that for the off time of our file of interest
    name2=NAME{2};
    second_file=sprintf('%s%d',name2,exten(2)) ;
    [on_time2,~]=find_NLXFILE_start(second_file);
    off_time2=on_time2;
end


%%%%%%% %%% neuralynx specific, ensures all the relevant events are read
% FieldSelection(1) = 1;
% FieldSelection(2) = 1;
% FieldSelection(3) = 1;
% FieldSelection(4) = 1;
% FieldSelection(5) = 1;
% FieldSelection(6) = 1;

ExtractMode = 4;
ModeArray(1) = on_time1(1);
ModeArray(2) = off_time2(end);
[TimeStamps_EV,  event_arr_NLX, EventStrings ] = Nlx2MatEV( 'Events.Nev', [1 0 1 0 1],0, ExtractMode, ModeArray ); %#ok<NASGU> %

%%% The encodes are buffered by zeros, ALWIN THINKS THAT IS MORE SECURE and he is probably correct,
%%% and here we eliminate the zeros
NonZeroEvents = find(event_arr_NLX~=0);
event_arr_NLX = event_arr_NLX(NonZeroEvents);
TimeStamps_EV = TimeStamps_EV(NonZeroEvents);

%%%%%%%%% Read the header of the continuous file (can be commented out since it is not used currently)
% ExtractHeader = 1;
% cont_files=dir;
% for k=3:size(cont_files,1)
%     prl=strmatch(cont_files(k).name,'LFP1.ncs');
%     if prl==1
%         contin_name=cont_files(k).name;
%     end
%     prl=strmatch(cont_files(k).name,'CSC1.Ncs');
%     if prl==1
%         contin_name=cont_files(k).name;
%     end
% end
%
% [NlxHeader] = Nlx2MatCSC(contin_name, [0 0 0 0 0 ], ExtractHeader,1 )


% %%%%%%%%%%%%% read the LFP in a format where each LFP data point is
% %%%%%%%%%%%%% associated with a time stamp, that is not the case in the
% %%%%%%%%%%%%% original format, only needed if LFP analysis or SFC  is done.
% [time_Arr_LFP, value_Arr, SampFreq]=read_nlx_CSC ('LFP1.ncs');
% prl=find(time_Arr_LFP>on_time1(1) & time_Arr_LFP>off_time2(end));
% time_Arr_LFP=time_Arr_LFP(prl);
% value_Arr=value_Arr(prl);
% prl=find(NlxHeader{15}==' ');
% prl2=NlxHeader{15};
% prl2=str2double(prl2(prl+1:end))
% Header.ADBitVolts=prl2
% value_Arr=value_Arr*(-1)*Header.ADBitVolts;

% load the bar file if it exists
if(exist([NLX_Path,'\CSCBar.ncs'],'file')==2 && Read_ReactionTimes==1)
    [Bar_time_Arr, Bar_value_Arr, Bar_SampFreq]=read_nlx_CSC ([NLX_Path,'\CSCBar.ncs']);
    BarFileExists = true;
else
    BarFileExists = false;
    disp('No bar file found. The time resolution for the reaction time will be limited to the refresh rate of the screen');
end

%%% Get the spikes from the sorted and split .NSE file e.g. SE_file='rvco4.50_bfm157.20.Nse'
[TimeStamps_SC,  CellNumbers] = Nlx2MatSpike(SE_file, [1 0 1 0 0],  0, ExtractMode, ModeArray );

%%%%% get the number of conditions and blocks %%%%%%%%
numblock=0;
numcond=0;

NLX_Trial_PARAM_START=find(event_arr_NLX==253);
for TRIAL=1:length(NLX_Trial_PARAM_START)-1
    EventsInTrial = event_arr_NLX(NLX_Trial_PARAM_START(TRIAL):NLX_Trial_PARAM_START(TRIAL+1));
    if EventsInTrial(5)==252 %& ~isempty(find(cycles==(header(2,TRIAL)+1))))
        cond_num=EventsInTrial(3);
        block = EventsInTrial(2);
        if cond_num>numcond
            numcond = cond_num;
        end
        if block>numblock
            numblock = block;
        end
    end
end



disp(['Number of blocks ',int2str(numblock)]);
disp(['Number of conditions ',int2str(numcond)]);

if numcond ==36 && numblock==4
    %%%%%%
    %     Dir 1:
    %     RF: 1:6
    %     Out1: 7:12
    %     Out2: 13:18
    %
    %     Dir2:
    %     RF: 19:24
    %     Out1: 25:30
    %     Out2: 31:36
    %%%%%%
    cnds=[1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4 5 5 5 5 5 5 6 6 6 6 6 6]; %this helps the assignment to variable locations later on, as cond 1-6 are identical, 7-12 are identical etc.
    blk_sorter=[1 2 1 2];
    blk_count=zeros(1,4);
    ass_blk_num=2;
    numcond=6;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% GET ON_TIME AND OFF_TIME IN CASE SPLIT IS NEEDED DEPENDENING ON CYCLE NUMBER
%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
currdir=cd;
if size(cycles,2)>1
    for file_no=1:numfiles %%% this cycles through the different data files, which is needed if drug low and high was done in different files, rather than blocks
        curr_file=sprintf('%s%d',name,exten(file_no)); %current file of interest in case there are multiple
        [on_time, off_time]=find_NLXFILE_start(curr_file);  % this finds the time stamps in neuralynx, which correspond to those data files, necessary to match spikes and LFP time stamps
        disp(['on_time ',int2str(on_time)]);
        disp(['off_time ',int2str(off_time)]);
        prl=find(currdir=='\'); %% find all subdirectories in the path (delineated by '\')
        prl=currdir(1:prl(end-1)-1);
        prl=sprintf('cd %s\\cortex',prl); %% move into the cortex subdirectory
        eval(prl); %% move into the cortex subdirectory
        %%%%%%%%%%%%%%% the next bit is specific to Jones, will be different
        %%%%%%%%%%%%%%% for Pascoe and has to be accomodated.......
        found_jones=strfind(lower(NLX_Path), 'jones');
        if ~isempty(found_jones)
            add_dir=dir;
            add_dir=sprintf('cd %s', add_dir(3).name); %%% moves one directory deeper as that is where Jones data hold the cortex recording file.
            eval(add_dir);%%% moves one directory deeper as that is where Jones data hold the cortex recording file.
            [cortex_time_arr,cortex_event_arr,eog_arr,epp_arr,eppnum,channel,header,trial]=readcortepp(curr_file); % %%%% read the cortex file, to cross check with neuralynx regarding number of trials, but will also be necessary to anlayse eye movements as these are not stored in neuralynx file.
        else
            [cortex_time_arr,cortex_event_arr,eog_arr,epp_arr,eppnum,channel,header,trial]=readcortepp(curr_file); % %%%% read the cortex file, to cross check with neuralynx regarding number of trials, but will also be necessary to anlayse eye movements as these are not stored in neuralynx file.
        end
        
        prl=sprintf('cd %s',currdir); %%% move back to Neuralynx data files
        eval(prl);
        curr_EV=event_arr_NLX(TimeStamps_EV>=on_time & TimeStamps_EV<=off_time); % now selecting the events
        curr_TimeEV=TimeStamps_EV(TimeStamps_EV>=on_time & TimeStamps_EV<=off_time); %time staps of events
        CurrT_SP1=TimeStamps_SC(TimeStamps_SC>=on_time & TimeStamps_SC<=off_time & CellNumbers==CELL(1)); % spike times of spikes of interest for cell 1
        
        if split_spike_timepoints==1
            on_time=CurrT_SP1(1);
            off_time=CurrT_SP1(end);
        end
        
        %%%% if there are two rows in the variable CELL find the second one,
        %%%% this is needed for e.g. noise correlation analysis
        if size(CELL,1)>1
            CurrT_SP2=TimeStamps_SC(TimeStamps_SC>=on_time & TimeStamps_SC<=off_time & CellNumbers==CELL(2)); % spike times of spikes of interest for cell 2
        end
        
        NLX_Trial_PARAM_START=find(curr_EV==253); % get a representation (starting point) of all the different trials that occurred
        
        %%%%%%% get a counter for where to start in the neuralynxfile and
        %%%%%%% cortexfile.
        if split_spike_timepoints==1
            prl1=find(curr_TimeEV(NLX_Trial_PARAM_START)>on_time);
            prl2=find(curr_TimeEV(NLX_Trial_PARAM_START)<off_time);
            trialcounts=[prl1(1) prl2(end)];
        else
            trialcounts=[1 length(NLX_Trial_PARAM_START)];
        end
        new_ontime=[];
        new_offtime=[];
        for TRIAL=trialcounts(1):trialcounts(2)-1  %cycle through all the trials (other than the last, for bracketing reasons)
            prl_ev= curr_EV(NLX_Trial_PARAM_START(TRIAL):NLX_Trial_PARAM_START(TRIAL+1)); %%% find the events that relate to current trials
            prl_EV_t=curr_TimeEV(NLX_Trial_PARAM_START(TRIAL):NLX_Trial_PARAM_START(TRIAL+1)); %%% find the events times that relate to current events
            if prl_ev(5)==252 %%% this makes sure that the condition/block encoding has the right number of events, otherwise skip the trial, as indeterminate
                %%%%% assign block and condition number
                %%% prl_ev(2) encodes the block
                %%% prl_ev(3) encodes the condition
                blk=blk_sorter(prl_ev(2));
                if last_blk==-1
                    last_blk=blk;
                    blk_CHG_trial_counter=elim_first_blkCHG_trials+1;%% ensures no trials are eliminated on experiment start
                    cycle=min(blk_count)+1;
                else
                    if blk~=last_blk
                        last_blk=blk; %%% theblock has changed
                        blk_CHG_trial_counter=1; %%% reset variable to one
                        blk_count(prl_ev(2))=blk_count(prl_ev(2))+1;
                        cycle=min(blk_count)+1;
                    else
                        blk_CHG_trial_counter=blk_CHG_trial_counter+1; %%%% increase by 1
                        cycle=min(blk_count)+1;
                    end
                end
                if isempty(new_ontime) && cycle==cycles(1)
                    new_ontime=prl_EV_t(1);
                end
                if isempty(new_offtime) && cycle>cycles(2)
                    new_offtime=prl_EV_t(1);
                end
            end
        end
    end
end
%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% end checking on and off time in case split for cycles
%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    

%%%% open a figure into which to plot
figinfo=[NLX_Path ' ' name ' ' num2str(exten) ' Cell ' num2str(CELL) ' Cycles ' num2str(cycles)];
pos=[0.01, 0.01, 0.85, 0.85];
fighandle1=  figure('Color',[1,1,1],'Units', 'Normalized', 'Position',pos); %
set(fighandle1, 'NumberTitle', 'off', 'Name', figinfo);
set(fighandle1, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.0 0.0 29.305 20.65]);



%%%%%%% trial counters, which keep track of how often a condition/block
%%%%%%% combination actually occurred. Important info, for e.g. average
%%%%%%% firing rates for a condition etc.
ndir=zeros(ass_blk_num,numcond+1);
ndir2=zeros(numfiles,numcond+1,numcond+1); % not sure where and if these are needed
ndir3=zeros(numfiles,numcond+1,numcond+1);
%%%%%%%%%%%% variables to contain the spiking information
sps=zeros(4,ass_blk_num,numcond+1,100); %% three spike analyses [3: at the beginning], one for after stim onset, one for after cue onset, one for before dimming; numcond +1 could be used to store spont activity
sps_varT=zeros(3,ass_blk_num,numcond+1,20,100); %% three spike analyses [3: at the beginning], one for after stim onset, one for after cue onset, one for before dimming; then 20 time resolved windows for time resolved FF (or Conditional variance

% akt_dir=zeros(ass_blk_num,numcond+1,100);
spontan=zeros(ass_blk_num,200);
spontan_n=zeros(ass_blk_num,1);

%%%%%%%%%%%%%% histograms not sure the most sensbel way of doing it, works
%%%%%%%%%%%%%% for now
histo_stim=zeros(ass_blk_num,numcond+1,100,760);
histo_cue=zeros(ass_blk_num,numcond+1,100,760);
histo_dim=zeros(ass_blk_num,numcond+1,100,760);
%%%%%%%%%%%%% Reaction time variable, stroed in CELLS, as likewise it does
RT=cell(1,36);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% go throught data once to get number of trials and be able to
%%%%%%%%%%% eliminate trials of extreme activity....
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for file_no=1:numfiles %%% this cycles through the different data files, which is needed if drug low and high was done in different files, rather than blocks
    curr_file=sprintf('%s%d',name,exten(file_no)); %current file of interest in case there are multiple
    [on_time, off_time]=find_NLXFILE_start(curr_file);  % this finds the time stamps in neuralynx, which correspond to those data files, necessary to match spikes and LFP time stamps
    disp(['on_time ',int2str(on_time)]);
    disp(['off_time ',int2str(off_time)]);
    prl=find(currdir=='\'); %% find all subdirectories in the path (delineated by '\')
    prl=currdir(1:prl(end-1)-1);
    prl=sprintf('cd %s\\cortex',prl); %% move into the cortex subdirectory
    eval(prl); %% move into the cortex subdirectory
    %%%%%%%%%%%%%%% the next bit is specific to Jones, will be different
    %%%%%%%%%%%%%%% for Pascoe and has to be accomodated.......
    found_jones=strfind(lower(NLX_Path), 'jones');
    prl=(find(curr_file==','));
    if ~isempty(prl)
        curr_file(prl)='.';
    end
    if ~isempty(found_jones)
        add_dir=dir;
        add_dir=sprintf('cd %s', add_dir(3).name); %%% moves one directory deeper as that is where Jones data hold the cortex recording file.
        eval(add_dir);%%% moves one directory deeper as that is where Jones data hold the cortex recording file.
        [cortex_time_arr,cortex_event_arr,eog_arr,epp_arr,eppnum,channel,header,trial]=readcortepp(curr_file); % %%%% read the cortex file, to cross check with neuralynx regarding number of trials, but will also be necessary to anlayse eye movements as these are not stored in neuralynx file.
    else
        [cortex_time_arr,cortex_event_arr,eog_arr,epp_arr,eppnum,channel,header,trial]=readcortepp(curr_file); % %%%% read the cortex file, to cross check with neuralynx regarding number of trials, but will also be necessary to anlayse eye movements as these are not stored in neuralynx file.
    end
    
    prl=sprintf('cd %s',currdir); %%% move back to Neuralynx data files
    eval(prl);
    curr_EV=event_arr_NLX(TimeStamps_EV>=on_time & TimeStamps_EV<=off_time); % now selecting the events
    curr_TimeEV=TimeStamps_EV(TimeStamps_EV>=on_time & TimeStamps_EV<=off_time); %time staps of events
    CurrT_SP1=TimeStamps_SC(TimeStamps_SC>=on_time & TimeStamps_SC<=off_time & CellNumbers==CELL(1)); % spike times of spikes of interest for cell 1
    
    if split_spike_timepoints==1
        on_time=CurrT_SP1(1);
        off_time=CurrT_SP1(end);
        
    end
    if size(cycles,2)>1
        if  new_ontime > on_time
            on_time=new_ontime;
        end
        if new_offtime < off_time
             off_time=new_offtime;
        end
    end
        
    
    %%%% if there are two rows in the variable CELL find the second one,
    %%%% this is needed for e.g. noise correlation analysis
    if size(CELL,1)>1
        CurrT_SP2=TimeStamps_SC(TimeStamps_SC>=on_time & TimeStamps_SC<=off_time & CellNumbers==CELL(2)); % spike times of spikes of interest for cell 2
    end
    
    NLX_Trial_PARAM_START=find(curr_EV==253); % get a representation (starting point) of all the different trials that occurred
    
    %%%%%%% get a counter for where to start in the neuralynxfile and
    %%%%%%% cortexfile.
    if split_spike_timepoints==1
        prl1=find(curr_TimeEV(NLX_Trial_PARAM_START)>on_time);
        prl2=find(curr_TimeEV(NLX_Trial_PARAM_START)<off_time);
        trialcounts=[prl1(1) prl2(end)];
    else
        trialcounts=[1 length(NLX_Trial_PARAM_START)];
    end
    
    for TRIAL=trialcounts(1):trialcounts(2)-1  %cycle through all the trials (other than the last, for bracketing reasons)
        prl_ev= curr_EV(NLX_Trial_PARAM_START(TRIAL):NLX_Trial_PARAM_START(TRIAL+1)); %%% find the events that relate to current trials
        prl_EV_t=curr_TimeEV(NLX_Trial_PARAM_START(TRIAL):NLX_Trial_PARAM_START(TRIAL+1)); %%% find the events times that relate to current events
        if prl_ev(5)==252 %%% this makes sure that the condition/block encoding has the right number of events, otherwise skip the trial, as indeterminate
            %%%%% assign block and condition number
            %%% prl_ev(2) encodes the block
            %%% prl_ev(3) encodes the condition
            blk=blk_sorter(prl_ev(2));
            if last_blk==-1
                last_blk=blk;
                blk_CHG_trial_counter=elim_first_blkCHG_trials+1;%% ensures no trials are eliminated on experiment start
            else
                if blk~=last_blk
                    last_blk=blk; %%% theblock has changed
                    blk_CHG_trial_counter=1; %%% reset variable to one
                    blk_count(prl_ev(2))=blk_count(prl_ev(2))+1;
                    cycle=min(blk_count)+1;
                else
                    blk_CHG_trial_counter=blk_CHG_trial_counter+1; %%%% increase by 1
                    %%%blk_count(blk)=blk_count(blk)+1;
                    cycle=min(blk_count)+1;
                end
            end
            if blk_CHG_trial_counter>elim_first_blkCHG_trials %%%%% only use trials if they are a certain distance from block change
                cond_num=cnds(1,prl_ev(3));
                %%% get spike times for this trial, needs to be done, before
                %%% aligning times of cortex and NLX files
                SP_times1=CurrT_SP1(CurrT_SP1>=prl_EV_t(1) &  CurrT_SP1<=prl_EV_t(end)); %  spike times in this trial
                if size(CELL,1)>1
                    SP_times2=CurrT_SP1(CurrT_SP2>=prl_EV_t(1) &  CurrT_SP2<=prl_EV_t(end)); %  spike times of second cell in this trial
                end
                %%%%% rather tedious, but we need to align the times of the cortex
                %%%%% file and the Neuralynxfile, i.e. also the times of the spikes
                %%%%% and LFP obviously, it seems that the second encode in the
                %%%%% cortex, corresponds to the first time stamp in the current
                %%%%% trial of the NLX file,
                time0 = prl_EV_t(1);
                prl_EV_t2=round((prl_EV_t-prl_EV_t(1))/1000);
                %cortex_time_arr(:,TRIAL)=cortex_time_arr(:,TRIAL)-cortex_time_arr(2,TRIAL);
                
                SP_times1=(SP_times1-prl_EV_t(1))/1000;
                if size(CELL,1)>1
                    SP_times2=(SP_times2-prl_EV_t(1))/1000;
                end
                prl_EV_t=prl_EV_t2;
                prl_ev=prl_ev(6:end);
                prl_EV_t=prl_EV_t(6:end);
                %%%% correct trials now (encoded by 5002, needs to be checked for different file types experiments %%%%%
                %if ~isempty(find(cortex_event_arr(:,TRIAL)==5002)) %%% it might be usefull to change this to neuralynx encodes, as then a mismatche between cortex and NLX doesn matter (unless eye movements are analysed)
                if ~isempty(find(prl_ev==104) ) %%% this is now using the NLX encodes rather than Cortex, correct trial indicator
                    stim_on_NLX=prl_EV_t(prl_ev==8);
                    bar_release=find(prl_ev==104);
                    target_dim=prl_EV_t(bar_release-1);
                    bar_release=prl_EV_t(bar_release);
                    if(BarFileExists)
                        bar_release_NLXtime = (bar_release*1000) + time0;   % calculate the bar release time in NLX time
                        Bar = ( (Bar_time_Arr>(bar_release_NLXtime-25000)) & (Bar_time_Arr<(bar_release_NLXtime+6000)) );   % get positions in bar channel around the bar release
                        Selected_Bar_Values = Bar_value_Arr(Bar);   % Get the amplitudes of the bar channel
                        Selected_Bar_times = Bar_time_Arr(Bar);     % Get the timestamps of the bar channel
                        BarMaxPoint = find(Selected_Bar_Values(:)>=(max(Selected_Bar_Values)/2),1,'first');     % find the first point of the bar channel where the maplitude exeedes half the maxium
                        Selected_Bar_times(BarMaxPoint);            % Get the timestamp for the BarMaxPoint
                        bar_release = (Selected_Bar_times(BarMaxPoint)-time0)/1000;
                    end
                    stim_on=min(prl_EV_t(prl_ev==8));
                    cue_on=min(prl_EV_t(prl_ev==20));
                    first_dim=min(prl_EV_t(prl_ev==25));
                    if ~isempty(stim_on) && ~isempty(cue_on) && ~isempty(first_dim)
                        ndir(blk,cond_num)=ndir(blk,cond_num)+1;
                        
                        x=SP_times1(SP_times1>(stim_on+t_win(1,1)) & SP_times1<(stim_on+t_win(1,2)) );
                        sps(1,blk, cond_num, ndir(blk,cond_num))=length(x);%*(1000/(t_win(1,2)-t_win(1,1)));
                        
                        x=SP_times1(SP_times1>(cue_on+t_win(1,1)) & SP_times1<(cue_on+t_win(1,2)) );
                        sps(2,blk, cond_num, ndir(blk,cond_num))=length(x);%*(1000/(t_win(2)-t_win(1)));
                        
                        x=SP_times1(SP_times1>(first_dim-t_win(2,1)) & SP_times1<(first_dim-t_win(2,2)) );
                        sps(3,blk, cond_num, ndir(blk,cond_num))=length(x);%*(1000/(t_win(2)-t_win(1)));
                        
                        x=SP_times1(SP_times1<(first_dim) & SP_times1>(first_dim-300));
                        sps(4,blk, cond_num, ndir(blk,cond_num))=length(x);%*(1000/(t_win(2)-t_win(1)));
                    end
                end
            end
        end
    end
end;
%%%%%%%%%%%%%%%%%%%
%%% done this first run through
%%%%%%%%%%%%%%%%%%
t_period_act_params=zeros(4,ass_blk_num,numcond, 2);
ndir_elim=zeros(ass_blk_num,numcond+1,150);
ndir_elim2=zeros(ass_blk_num,numcond+1);%%%% this one is for second round, to determine whether current trial is a trial that needs elimination...

for kl=1:4
    for lm=1:ass_blk_num
        for mn=1:numcond
            t_period_act_params(kl,lm, mn,1)=mean(sps(kl,lm,mn,1:ndir(lm,mn)));
            t_period_act_params(kl,lm, mn,2)=std(sps(kl,lm,mn,1:ndir(lm,mn)));
            prl=find(sps(kl,lm,mn,1:ndir(lm,mn))>t_period_act_params(kl,lm, mn,1)+cutoff_act*t_period_act_params(kl,lm, mn,2));
            if ~isempty(prl)
                ndir_elim(lm,mn,prl)=1;
            end
        end
    end
end

ndir=zeros(ass_blk_num,numcond+1);
ndir2=zeros(numfiles,numcond+1,numcond+1); % not sure where and if these are needed
ndir3=zeros(numfiles,numcond+1,numcond+1);
last_blk=-1; %%% reset last_blk, as we start from file beginning again.
blk_CHG_trial_counter=1; %%% reset counter
%%%%%%%%%%%% variables to contain the spiking information
sps=zeros(4,ass_blk_num,numcond+1,100); %% three spike analyses [3: at the beginning], one for after stim onset, one for after cue onset, one for before dimming; numcond +1 could be used to store spont activity

currdir=cd;
for file_no=1:numfiles %%% this cycles through the different data files, which is needed if drug low and high was done in different files, rather than blocks
    curr_file=sprintf('%s%d',name,exten(file_no)); %current file of interest in case there are multiple
    [on_time, off_time]=find_NLXFILE_start(curr_file);  % this finds the time stamps in neuralynx, which correspond to those data files, necessary to match spikes and LFP time stamps
    disp(['on_time ',int2str(on_time)]);
    disp(['off_time ',int2str(off_time)]);
    prl=find(currdir=='\'); %% find all subdirectories in the path (delineated by '\')
    prl=currdir(1:prl(end-1)-1);
    prl=sprintf('cd %s\\cortex',prl); %% move into the cortex subdirectory
    eval(prl); %% move into the cortex subdirectory
    %%%%%%%%%%%%%%% the next bit is specific to Jones, will be different
    %%%%%%%%%%%%%%% for Pascoe and has to be accomodated.......
    found_jones=strfind(lower(NLX_Path), 'jones');
    prl=(find(curr_file==','));
    if ~isempty(prl)
        curr_file(prl)='.';
    end
    if ~isempty(found_jones)
        add_dir=dir;
        add_dir=sprintf('cd %s', add_dir(3).name); %%% moves one directory deeper as that is where Jones data hold the cortex recording file.
        eval(add_dir);%%% moves one directory deeper as that is where Jones data hold the cortex recording file.
        [cortex_time_arr,cortex_event_arr,eog_arr,epp_arr,eppnum,channel,header,trial]=readcortepp(curr_file); % %%%% read the cortex file, to cross check with neuralynx regarding number of trials, but will also be necessary to anlayse eye movements as these are not stored in neuralynx file.
    else
        [cortex_time_arr,cortex_event_arr,eog_arr,epp_arr,eppnum,channel,header,trial]=readcortepp(curr_file); % %%%% read the cortex file, to cross check with neuralynx regarding number of trials, but will also be necessary to anlayse eye movements as these are not stored in neuralynx file.
    end
    
    %      for jk=1:size(header,2)
    %         if (header(13,jk)==0 & ~isempty(find(cortex_event_arr(:,jk)==4001))  )%& ~isempty(find(cycles==(header(2,i)+1))))
    %             cond_num=cnds(1,header(1,jk)+1);
    %             blk_num=header(3,jk)+1;
    %             Tdim=cortex_time_arr((find(cortex_event_arr(:,jk)==4002)),jk);
    %             %%% reaction times %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %             RT=cortex_time_arr((find(cortex_event_arr(:,jk)==5002)),jk)-Tdim;
    %        end
    %     end
    %%%%%%
    
    prl=sprintf('cd %s',currdir); %%% move back to Neuralynx data files
    eval(prl);
    curr_EV=event_arr_NLX(TimeStamps_EV>=on_time & TimeStamps_EV<=off_time); % now selecting the events
    curr_TimeEV=TimeStamps_EV(TimeStamps_EV>=on_time & TimeStamps_EV<=off_time); %time staps of events
    CurrT_SP1=TimeStamps_SC(TimeStamps_SC>=on_time & TimeStamps_SC<=off_time & CellNumbers==CELL(1)); % spike times of spikes of interest for cell 1
    %%%%%%%%%%%%% 5 if you want to split the file according to
    %%%%%%%%%%%%% sorting, because the cell did not come into
    %%%%%%%%%%%%% existence until later in the file or ceased
    %%%%%%%%%%%%% before the end then this needs special
    %%%%%%%%%%%%% attention, i.e. reselection the ontime and
    %%%%%%%%%%%%% offtime, reselecting the trials of interest,
    %%%%%%%%%%%%% etc. this also needs to be done for the cortex
    %%%%%%%%%%%%% file, as it is currently used for certain encodes.
    
    if split_spike_timepoints==1
        on_time=CurrT_SP1(1);
        off_time=CurrT_SP1(end);
    end
    if size(cycles,2)>1
        if  new_ontime > on_time
            on_time=new_ontime;
        end
        if new_offtime < off_time
            off_time=new_offtime ;
        end
    end
    %%%% if there are two rows in the variable CELL find the second one,
    %%%% this is needed for e.g. noise correlation analysis
    if size(CELL,1)>1
        CurrT_SP2=TimeStamps_SC(TimeStamps_SC>=on_time & TimeStamps_SC<=off_time & CellNumbers==CELL(2)); % spike times of spikes of interest for cell 2
    end
    
    %%% LFP analysis currently not done, so commented out%%%%%%%%%%%%%%%%%%%%%%
    %     CurrTLFP=time_Arr_LFP(time_Arr_LFP>=on_time & time_Arr_LFP<=off_time);
    %     CurrVLFP=value_Arr(time_Arr_LFP>=on_time & time_Arr_LFP<=off_time);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    NLX_Trial_PARAM_START=find(curr_EV==253); % get a representation (starting point) of all the different trials that occurred
    if(size(cortex_time_arr,2)~=length(NLX_Trial_PARAM_START))
        disp('cortex-NLX trial mismatch --- careful!!!!')
        %pause
    end
    %%%%%%% get a counter for where to start in the neuralynxfile and
    %%%%%%% cortexfile.
    if split_spike_timepoints==1
        prl1=find(curr_TimeEV(NLX_Trial_PARAM_START)>on_time);
        prl2=find(curr_TimeEV(NLX_Trial_PARAM_START)<off_time);
        trialcounts=[prl1(1) prl2(end)];
    else
        trialcounts=[1 length(NLX_Trial_PARAM_START)];
    end
    
    sigma=10; %%% used for smoothing of histograms
    count=0;
    t_stim=(-250:1000); %%%
    t_dim=(-1000:250);
    y=ones(1,length(t_stim));
    
    PlotDataX = cell(24,2);
    PlotDataY = cell(24,2);
    
    for ij=1:24
        PlotDataX{ij,1} = [];
        PlotDataX{ij,2} = [];
        PlotDataY{ij,1} = [];
        PlotDataY{ij,2} = [];
    end
    
    
    for TRIAL=trialcounts(1):trialcounts(2)-1  %cycle through all the trials (other than the last, for bracketing reasons)
        prl_ev= curr_EV(NLX_Trial_PARAM_START(TRIAL):NLX_Trial_PARAM_START(TRIAL+1)); %%% find the events that relate to current trials
        prl_EV_t=curr_TimeEV(NLX_Trial_PARAM_START(TRIAL):NLX_Trial_PARAM_START(TRIAL+1)); %%% find the events times that relate to current events
        
        
        if prl_ev(5)==252 %%% this makes sure that the condition/block encoding has the right number of events, otherwise skip the trial, as indeterminate
            %%%%% assign block and condition number
            %%% prl_ev(2) encodes the block
            %%% prl_ev(3) encodes the condition
            blk=blk_sorter(prl_ev(2));
            
            if last_blk==-1
                blk_CHG_trial_counter=elim_first_blkCHG_trials+1;%% ensures no trials are eliminated on experiment start
                last_blk=blk;
            else
                if blk~=last_blk
                    last_blk=blk; %%% theblock has changed
                    blk_CHG_trial_counter=1; %%% reset variable to one
                else
                    blk_CHG_trial_counter=blk_CHG_trial_counter+1; %%%% increase by 1
                end
            end
            if blk_CHG_trial_counter>elim_first_blkCHG_trials %%%%% only use trials if they are a certain distance from block change
                
                
                
                
                RT_condnum=prl_ev(3);
                cond_num=cnds(1,prl_ev(3));
                %%% get spike times for this trial, needs to be done, before
                %%% aligning times of cortex and NLX files
                SP_times1=CurrT_SP1(CurrT_SP1>=prl_EV_t(1) &  CurrT_SP1<=prl_EV_t(end)); %  spike times in this trial
                if size(CELL,1)>1
                    SP_times2=CurrT_SP1(CurrT_SP2>=prl_EV_t(1) &  CurrT_SP2<=prl_EV_t(end)); %  spike times of second cell in this trial
                end
                %%%%% rather tedious, but we need to align the times of the cortex
                %%%%% file and the Neuralynxfile, i.e. also the times of the spikes
                %%%%% and LFP obviously, it seems that the second encode in the
                %%%%% cortex, corresponds to the first time stamp in the current
                %%%%% trial of the NLX file,
                time0 = prl_EV_t(1);
                prl_EV_t2=round((prl_EV_t-prl_EV_t(1))/1000);
                %cortex_time_arr(:,TRIAL)=cortex_time_arr(:,TRIAL)-cortex_time_arr(2,TRIAL);
                
                SP_times1=(SP_times1-prl_EV_t(1))/1000;
                if size(CELL,1)>1
                    SP_times2=(SP_times2-prl_EV_t(1))/1000;
                end
                prl_EV_t=prl_EV_t2;
                %             [prl_EV_t; prl_ev]
                %             [cortex_time_arr(:,TRIAL) cortex_event_arr(:,TRIAL)]
                %             pause
                %%%% now remove the first five encoders, as the info there can be
                %%%% read as events if not appropriately controlled for,
                prl_ev=prl_ev(6:end);
                prl_EV_t=prl_EV_t(6:end);
                %%%% in principle all the sorting can now be done based on the cortex encoded events.....
                %%%%% NLX reaction time events have a real problem as they have
                %%%%% a resolution of 10ms which is way too coars. need to sort
                %%%%% why that is, but try to use Cortex events for now.
                
                %%%% correct trials now (encoded by 5002, needs to be checked for different file types experiments %%%%%
                %if ~isempty(find(cortex_event_arr(:,TRIAL)==5002)) %%% it might be usefull to change this to neuralynx encodes, as then a mismatche between cortex and NLX doesn matter (unless eye movements are analysed)
                if ~isempty(find(prl_ev==104) ) %%% this is now using the NLX encodes rather than Cortex, correct trial indicator
                    stim_on_NLX=prl_EV_t(prl_ev==8);
                    bar_release=find(prl_ev==104);
                    target_dim=prl_EV_t(bar_release-1);
                    bar_release=prl_EV_t(bar_release);
                    
                    if(BarFileExists)
                        bar_release_NLXtime = (bar_release*1000) + time0;   % calculate the bar release time in NLX time
                        Bar = ( (Bar_time_Arr>(bar_release_NLXtime-25000)) & (Bar_time_Arr<(bar_release_NLXtime+6000)) );   % get positions in bar channel around the bar release
                        Selected_Bar_Values = Bar_value_Arr(Bar);   % Get the amplitudes of the bar channel
                        Selected_Bar_times = Bar_time_Arr(Bar);     % Get the timestamps of the bar channel
                        BarMaxPoint = find(Selected_Bar_Values(:)>=(max(Selected_Bar_Values)/2),1,'first');     % find the first point of the bar channel where the maplitude exeedes half the maxium
                        Selected_Bar_times(BarMaxPoint);            % Get the timestamp for the BarMaxPoint
                        bar_release = (Selected_Bar_times(BarMaxPoint)-time0)/1000;
                    end
                    
                    React_t=bar_release-target_dim;
                    RT{RT_condnum}=[RT{RT_condnum} React_t];
                    stim_on=min(prl_EV_t(prl_ev==8));
                    cue_on=min(prl_EV_t(prl_ev==20));
                    first_dim=min(prl_EV_t(prl_ev==25));
                    if ~isempty(stim_on) && ~isempty(cue_on) && ~isempty(first_dim)
                        ndir_elim2(blk,cond_num)= ndir_elim2(blk,cond_num)+1;
                        if ndir_elim(blk,cond_num,ndir_elim2(blk,cond_num))~=1; %%% this eliminates all trials that have activity exceeding the allowed cutoff
                            
                            %stim_on_CTX=cortex_time_arr(cortex_event_arr(:,TRIAL)==4001,TRIAL);
                            %stim_on_CTX=stim_on_CTX(1);
                            %%%% variable that encodes the number of correct trials,
                            %%%% error rate thus the # of errortrials divided by this
                            %%%% variable
                            %%% the different ndir are used for differential couning of how
                            %%% often a stimulus has been presented etc.
                            ndir(blk,cond_num)=ndir(blk,cond_num)+1;
                            ndir2(file_no,blk,cond_num)=ndir(blk,cond_num);
                            ndir3(file_no,blk,cond_num)=ndir3(file_no,blk,cond_num)+1;
                            
                            %%% Event 2 stands for record start in Cortex, which is when eye pdata are acquired,
                            %%%T useful time zero if ever an eye data analysis is being done
                            
                            %%%%% in cortex a stimulus on was encoded as 4001, thus
                            %%%%% first_4001 stands for stimulus onset. In Neuralynx that
                            %%%%% is usually encoded by 8, here I call it stim_on;
                            
                            
                            % get the spontaneous activity (500 ms before stimulus
                            % onset
                            %                 s=SP_times1(SP_times1<stim_on & SP_times1<stim_on-500);
                            %                 % define colors for att RF att away and no drug - drug....
                            %                 if  blk==1 && (cond_num==1 || cond_num==4)
                            %                     col=[1.0 0.0 0.0];
                            %                 elseif blk==1 && (cond_num==2 || cond_num==5 || cond_num==3 || cond_num==6)
                            %                     col=[0.0 0.0 1.0];
                            %                 elseif blk==2 && (cond_num==1 || cond_num==4 )
                            %                     col=[0.0 1.0 0.0];
                            %                 elseif blk==2 && (cond_num==2 || cond_num==5 || cond_num==3 || cond_num==6)
                            %                     col=[0.0 0.0 0.0];
                            %                 end
                            
                            %%% target dim time
                            % #define NLX_TESTDIMMED         17
                            % #define NLX_DISTDIMMED         18
                            % #define NLX_BARRELEASED        19
                            % #define NLX_CUE_ON             20
                            % #define NLX_CUE_OFF            21
                            % #define NLX_DIST1DIMMED        22
                            % #define NLX_DIST2DIMMED        23
                            % #define NLX_SACCADE_START      24
                            % #define NLX_DIMMING1	         25
                            % #define	NLX_DIMMING2	       26
                            % #define	NLX_DIMMING3           27
                            % #define NLX_MICRO_STIM	       28
                            
                            
                            
                            x=SP_times1(SP_times1>(stim_on+t_win(1,1)) & SP_times1<(stim_on+t_win(1,2)) );
                            sps(1,blk, cond_num, ndir(blk,cond_num))=length(x);%*(1000/(t_win(1,2)-t_win(1,1)));
                            for jkl=1:20
                                x=SP_times1(SP_times1>(stim_on-300+(jkl-1)*60) & SP_times1<stim_on-300+(jkl)*60 );
                                sps_varT(1,blk, cond_num,jkl, ndir(blk,cond_num))=length(x);
                            end
                            x=SP_times1(SP_times1>(cue_on+t_win(1,1)) & SP_times1<(cue_on+t_win(1,2)) );
                            sps(2,blk, cond_num, ndir(blk,cond_num))=length(x);%*(1000/(t_win(2)-t_win(1)));
                            for jkl=1:20
                                x=SP_times1(SP_times1>(cue_on-300+(jkl-1)*60) & SP_times1<cue_on-300+(jkl)*60 );
                                sps_varT(2,blk, cond_num,jkl, ndir(blk,cond_num))=length(x);
                            end
                            x=SP_times1(SP_times1>(first_dim-t_win(2,1)) & SP_times1<(first_dim-t_win(2,2)) );
                            sps(3,blk, cond_num, ndir(blk,cond_num))=length(x);%*(1000/(t_win(2)-t_win(1)));
                            for jkl=1:20
                                x=SP_times1(SP_times1>(first_dim-900+(jkl-1)*60) & SP_times1<first_dim-900+(jkl)*60 );
                                sps_varT(3,blk, cond_num,jkl, ndir(blk,cond_num))=length(x);
                            end
                            
                            
                            
                            
                            
                            x=SP_times1(SP_times1<(stim_on) & SP_times1>(stim_on-300));
                            spontan_n(blk,1)=spontan_n(blk,1)+1;
                            spontan(blk,spontan_n(blk,1))=length(x);
                            
                            %%%% these are just relevant to plot spikes and generate
                            %%%% histograms
                            spikes_stim=SP_times1(SP_times1>(stim_on-250) & SP_times1<(stim_on+1000))-stim_on;
                            spikes_cue=SP_times1(SP_times1>(cue_on-250) & SP_times1<(cue_on+1000))-cue_on;
                            spikes_dim=SP_times1(SP_times1>(first_dim-1000) & SP_times1<(first_dim+250))-first_dim;
                            
                            if (~isempty(spikes_stim))
                                count=count+1;
                                %%%%% 18 main plots (attention*direction) with 2 different
                                %%%%% histograms (drug/no drug), for stim on, for
                                %%%%% cue on, and for first dim aligned spikes
                                
                                %%%%%%%% First do stim related activity
                                trialdata=zeros(1,length(t_stim));
                                for j=1:length(spikes_stim)
                                    if (plot_rasters==1 && (blk==1 || blk==2))
                                        PlotDataX{cond_num,blk} = [PlotDataX{cond_num,blk},spikes_stim(j)];
                                        PlotDataY{cond_num,blk} = [PlotDataY{cond_num,blk},ndir(blk,cond_num)*2];
                                    end
                                    spike=(y*(1/(sigma*sqrt(2*pi))).*exp(-(((t_stim-spikes_stim(j)).^2)/(2*sigma^2))));
                                    trialdata(1,:)=trialdata(1,:)+spike;
                                end;
                                %%% histrograms are initially generated for each trial
                                %%% individually, so it is possible to calculate SEM
                                %%% etc if ever needed
                                histo_stim(blk,cond_num,ndir(blk,cond_num),1:length(trialdata(1,:)))=trialdata(1,:);
                                
                                %%%%%%%%%%%%%%%% now do cue related activity
                                trialdata=zeros(1,length(t_stim));
                                for j=1:length(spikes_cue)
                                    if (plot_rasters==1&& (blk==1 || blk==2))
                                        PlotDataX{cond_num+(numcond+2),blk} = [PlotDataX{cond_num+(numcond+2),blk},spikes_cue(j)];
                                        PlotDataY{cond_num+(numcond+2),blk} = [PlotDataY{cond_num+(numcond+2),blk},ndir(blk,cond_num)*2];
                                    end
                                    spike=(y*(1/(sigma*sqrt(2*pi))).*exp(-(((t_stim-spikes_cue(j)).^2)/(2*sigma^2))));
                                    trialdata(1,:)=trialdata(1,:)+spike;
                                end;
                                histo_cue(blk,cond_num,ndir(blk,cond_num),1:length(trialdata(1,:)))=trialdata(1,:);
                                
                                %%%%%%%%%%%%%% And now do dim related activity
                                trialdata=zeros(1,length(t_dim));
                                for j=1:length(spikes_dim)
                                    if (plot_rasters==1 && (blk==1 || blk==2))
                                        PlotDataX{cond_num+2*(numcond+2),blk} = [PlotDataX{cond_num+2*(numcond+2),blk},spikes_dim(j)];
                                        PlotDataY{cond_num+2*(numcond+2),blk} = [PlotDataY{cond_num+2*(numcond+2),blk},ndir(blk,cond_num)*2];
                                    end
                                    spike=(y*(1/(sigma*sqrt(2*pi))).*exp(-(((t_dim-spikes_dim(j)).^2)/(2*sigma^2))));
                                    trialdata(1,:)=trialdata(1,:)+spike;
                                end;
                                histo_dim(blk,cond_num,ndir(blk,cond_num),1:length(trialdata(1,:)))=trialdata(1,:);
                            end
                        end
                    end
                end
            end
        end
    end
end;

%%%%%%%%%%%%%%%% do drift analysis
%%%%% correlation first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_vals_correlation=ones(3,2,numcond);
for kl=1:3
    for lm=1:ass_blk_num
        for mn=1:numcond
            curr_akt=squeeze(sps(kl,lm,mn,1:ndir(lm,mn)));
            prl=[1:1:length(curr_akt)];
            [r,p]=corrcoef(curr_akt',prl);
            P_vals_correlation(kl,lm,mn)=p(2);
        end
    end
end
prl=P_vals_correlation(:);
corr_prctile=prctile(prl,[25]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% plot the subplots and the raster plots
y_off=100;
for PlotNum = 1:(3*(numcond+2))
    test=subplot(3,numcond+2,PlotNum);
    hold on
    if (plot_rasters==1)
        switch PlotNum
            case {1,4,9,12,17,20}
                col=[1.0 0.0 0.0];
                col2=[0.0 1.0 0.0];
            otherwise
                col=[0.0 0.0 1.0];
                col2=[0.0 0.0 0.0];
        end
        
        for k=1:length(PlotDataX{PlotNum,1})
            line( [PlotDataX{PlotNum,1}(k),PlotDataX{PlotNum,1}(k)], [PlotDataY{PlotNum,1}(k)-0.5+y_off,PlotDataY{PlotNum,1}(k)+0.5+y_off] ,'color',col);
        end
        
        if (PlotNum>16)
            cond_num = PlotNum - 16;
        elseif (PlotNum>8)
            cond_num = PlotNum - 8;
        else
            cond_num = PlotNum;
        end
        
        if (cond_num < 7)
            yShift = max(ndir)*2 + 10;
            yShift=yShift(1);
        else
            yShift = 0;
        end
        
        for k=1:length(PlotDataX{PlotNum,2})
            line( [PlotDataX{PlotNum,2}(k),PlotDataX{PlotNum,2}(k)], [PlotDataY{PlotNum,2}(k)-0.5+y_off+yShift,PlotDataY{PlotNum,2}(k)+0.5+y_off+yShift] ,'color',col2);
        end
    end
end




disp('done with spikes')


%%%%%%% do mean analyses, histograms ANOVAs etc
meandata=zeros(3,2,numcond,length(t_stim));
col=['b' 'r' 'k' 'g'];
norm=0;

for DRUG=1:2 %%% separately for drug off (h=1)and drug on (h=2)...
    for COND=1:numcond
        nHight = ndir(DRUG,COND);
        nLength = length(t_stim);
        plotdata(1:nHight,:)=histo_stim(DRUG,COND,1:nHight,1:nLength);
        
        %%% calculate mean for this block/condition histrogram STIM related activity
        plotdata=mean(plotdata);
        if sum(plotdata)>0
            meandata(1,DRUG,COND,1:nLength)=plotdata;
            gnorm=max(plotdata);
            %%%% needed to normalize all histos to appropriate height (i.e.
            %%%% approp. scaling)
            if gnorm>norm
                norm=gnorm;
            end
        else
            meandata(1,DRUG,COND,1:nLength)=0;
        end
        
        %%% calculate mean for this block/condition histrogram CUE related activity
        plotdata(1:nHight,:)=histo_cue(DRUG,COND,1:nHight,1:nLength);
        plotdata=mean(plotdata);
        if sum(plotdata)>0
            meandata(2,DRUG,COND,1:nLength)=plotdata;
            gnorm=max(plotdata);
            if gnorm>norm
                norm=gnorm;
            end
        else
            meandata(2,DRUG,COND,1:nLength)=0;
        end
        
        %%% calculate mean for this block/condition histrogram DIM related activity
        plotdata(1:nHight,:)=histo_dim(DRUG,COND,1:nHight,1:nLength);
        plotdata=mean(plotdata);
        if sum(plotdata)>0
            meandata(3,DRUG,COND,1:nLength)=plotdata;
            gnorm=max(plotdata);
            if gnorm>norm
                norm=gnorm;
            end
        else
            meandata(3,DRUG,COND,1:nLength)=0;
        end
    end
end

% norm
% prl=max(squeeze(meandata(:)))

%%%% do the nornmalization
meandata=(meandata/norm)*100;

% prl=max(squeeze(meandata(:)))
% pause
col3=['b*';'r*'];

%%% calculate the mean spontaneous activities
mean_sp_1=mean(spontan(1,1:spontan_n(1,1)))/(300/1000);
mean_sp_2=mean(spontan(2,1:spontan_n(2,1)))/(300/1000);


blk_text=[{'att out lo ACh'};{'att in lo ACh'};{'att out hi ACh'};{'att in hi ACh'}];

FF_Tres=zeros(3,2,numcond,20);
for h=1:2 % h=1 no drug, h=2 drug
    for COND=1:numcond
        test=subplot(3,numcond+2,COND);
        if COND==1
            ylabel('aligned to stimulus onset');
        end
        hold on
        prl= squeeze(meandata(1,h,COND,1:length(t_stim)));
        %%%%% this does the smoothing of histograms
        plotdata=gaussfit(30,0,prl');
        
        %%%% assigne the colors related to attend RF attend away, drug no
        %%%% drug....
        if  h==1 && (COND==1 || COND==4)
            col=[1.0 0.0 0.0]; %% red
        elseif h==1 && (COND==2 || COND==5 || COND==3 || COND==6)
            col=[0.0 0.0 1.0]; %% blue
        elseif h==2 && (COND==1 || COND==4 )
            col=[0.0 1.0 0.0]; %% green
        elseif h==2 && (COND==2 || COND==5 || COND==3 || COND==6)
            col=[0.0 0.0 0.0]; %% black
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% this is for time resolved FF
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for jkl=1:20
            prl=squeeze(sps_varT(1,h,COND,jkl, 1:ndir(h,COND)));
            FF_Tres(1,h,COND,jkl)=mean(prl);
            FF_Tres(2,h,COND,jkl)=var(prl);
            prl=squeeze(sps_varT(2,h,COND,jkl, 1:ndir(h,COND)));
            FF_Tres(3,h,COND,jkl)=mean(prl);
            FF_Tres(4,h,COND,jkl)=var(prl);
            prl=squeeze(sps_varT(3,h,COND,jkl, 1:ndir(h,COND)));
            FF_Tres(5,h,COND,jkl)=mean(prl);
            FF_Tres(6,h,COND,jkl)=var(prl);
        end
        
        stim_bins=[-270:60:-270+19*60];
        dim_bins=[-870:60:-870+19*60];
        %%%%%%%%%%%%%% plot the STIM histograms
        plot(t_stim,plotdata','color',col);
        hold on;
%         prl=squeeze(FF_Tres(1,h,COND,:));
%         plot(stim_bins,(prl/3)*100,'--','color',col);
        set(test,'XLim',[-250 1000])
        if plot_rasters==0
            set(test,'YLim',[0 100])
        end
        set(test,'YTick',[0 100],'YTicklabel',[0 round(norm*1000)]);
        
        
        %%%%%%%%%%%%%% plot the CUE histograms
        test=subplot(3,numcond+2,COND+(numcond+2));
        hold on
        if COND==1
            ylabel('aligned to cue onset');
        end
        prl= squeeze(meandata(2,h,COND,1:length(t_stim)));
        plotdata=gaussfit(30,0,prl');
        plot(t_stim,plotdata','color',col);
        hold on;
%         prl=squeeze(FF_Tres(2,h,COND,:));
%         plot(stim_bins,(prl/3)*100,'--','color',col);
        set(test,'XLim',[-250 1000]);
        if plot_rasters==0
            set(test,'YLim',[0 100]);
        end
        set(test,'YTick',[0 100],'YTicklabel',[0 round(norm*1000)]);
        
        %%%%%%%%%%%%%% plot the DIM histograms
        test=subplot(3,numcond+2,COND+2*(numcond+2));
        hold on
        if COND==1
            ylabel('aligned to dim onset');
        end
        prl= squeeze(meandata(3,h,COND,1:length(t_stim)));
        plotdata=gaussfit(30,0,prl');
        plot(t_dim,plotdata','color',col);
        hold on;
%         prl=squeeze(FF_Tres(3,h,COND,:));
%         plot(dim_bins,(prl/3)*100,'--','color',col);
        set(test,'XLim',[-1000 250]);
        if plot_rasters==0
            set(test,'YLim',[0 100]);
        end
        set(test,'YTick',[0 100],'YTicklabel',[0 round(norm*1000)]);
    end
end

%%%%%%%% now write the info relating to the recording, drug,  and cell.

prl1=strfind(NLX_Path,'\');
prl2=strfind(NLX_Path,'\pen');
prl3=find(prl2==prl1);
if ((length(NLX_Path)>prl2+6) && (prl3>2))
    monkeytext = NLX_Path(prl1(prl3-2)+1:prl2+6);
else
    monkeytext = 'No text avalible';
end

test=subplot(3,numcond+2,1);
title(monkeytext);
test=subplot(3,numcond+2,2);
celltext=sprintf('File: %s%d Cell: %d',name,exten(1), CELL(1,1));
title(celltext);
test=subplot(3,numcond+2,3);
drugtext2=sprintf('Drug: %s ',drugtext);
title(drugtext2);

%%%DO the ANOVA analysis
col=['b^';'r^';'k^';'g^'];
col2=['b';'r';'k';'g'];
col3=['b*';'r*';'k*';'g*'];
L_one=['att';'no1';'no1';'att';'no1';'no1'];
L_two=['hold';'drug';'hold';'drug'];
s_dir=['1';'1';'1';'2';'2';'2'];
time_label=['STIM';'CUE ';'DIM '];
for j=1:3 %%%% three different time windows (stim, cue, dim)
    all_data=[];
    all_g1=[];
    all_g2=[];
    g_three=[];
    Adata_all=[];
    fac1_all=[];
    fac2_all=[];
    fac3_all=[];
    ranked_data=[];
    for COND=1:numcond
        Adata=[];
        group=[];
        g_one=[];
        g_two=[];
        count=0;
        for h=1:2
            for n=1:ndir(h,COND)
                Adata=[Adata sps(j,h,COND,n)];
                g_one=[g_one ;L_one(COND,:)];
                g_two=[g_two ;L_two(h,:)];
                fac1_all=[fac1_all;L_one(COND,:)];
                fac2_all=[fac2_all;L_two(h,:)];
                fac3_all=[fac3_all;s_dir(COND,:)];
                Adata_all=[Adata_all sps(j,h,COND,n)];
                count=count+1;
            end
        end
        group={g_one;g_two};
        all_data=[all_data Adata];
        all_g1=[all_g1;g_one];
        all_g2=[all_g2;g_two];
    end
    
    all_groups={all_g1;all_g2;fac3_all};
    txtcell{1,1}=sprintf('att');
    txtcell{2,1}=sprintf('dru');
    txtcell{3,1}=sprintf('dir');
    [P,table,stats,terms]  = anovan(all_data,all_groups,'full',3,txtcell,'off');
    overall_P=P;
    
    
    %%%%%%%%%%%%%%%% write ANOVA results into the plot
    test=subplot(3,4,j*(4));
    set(test,'Xlim',[-2 110],'Ylim',[0 11]);
    axis off
    counter=10;
    %%%% Format the ANOVA table output
    for kl=1:10
        for lm=1:7
            currtxt=table{kl,lm};
            if kl>1 && lm==2
                currtxt=sprintf('%6.1f',currtxt);
            end
            if kl>1 && (lm==3 || lm==4)
                currtxt=sprintf('%d',currtxt);
            end
            if kl>1 && (lm==5)
                currtxt=sprintf('%5.1f',currtxt);
            end
            if kl>1 && (lm==6)
                currtxt=sprintf('%5.2f',currtxt);
            end
            if kl>1 && (lm==7)
                currtxt=sprintf('%4.3f',currtxt);
            end
            if lm==1 && kl==1
                text('Position',[0 counter-(kl-1) ],'Fontsize',[12], 'String',time_label(j,:));
            elseif lm==1 && kl>1
                text('Position',[0 counter-(kl-1) ],'Fontsize',[8], 'String',currtxt);
            elseif kl==1 && lm==4
                text('Position',[9+(lm-1)*17 counter-(kl-1) ],'Fontsize',[8], 'String','Sing.');
            elseif kl==1 && lm==5
                text('Position',[9+(lm-1)*17 counter-(kl-1) ],'Fontsize',[8], 'String','Mean Sq.');
            elseif kl==1 && lm==6
                text('Position',[15+(lm-1)*17 counter-(kl-1) ],'Fontsize',[8], 'String','F');
            else
                text('Position',[12+(lm-1)*17 counter-(kl-1) ],'Fontsize',[8], 'String',currtxt);
                hold on
            end
        end
        
    end
    %text('Position',[1 5],'Fontsize',[4], 'String',table);
end


%%%%%%%%%%%%%%% assign whether attention, drug etc had significant effect
%%%%%%%%%%%%%%%% for further writing to files....
attP=0.9;
drugP=0.9;
intAttDrug=0.9;
if P(1)<0.05 | P(4)<0.05 |P(5)<0.05 | P(7)<0.05
    attP=0.01;
end
if P(2)<0.05 | P(4)<0.05 |P(6)<0.05 | P(7)<0.05
    drugP=0.01;
end
if  P(4)<0.05 | P(7)<0.05
    intAttDrug=0.01;
end


fighandle2=  figure('Color',[1,1,1],'Units', 'Normalized', 'Position',pos); %
set(fighandle2, 'NumberTitle', 'off', 'Name', figinfo);
set(fighandle2, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.0 0.0 27 8]);


%%%%%%%%%%%%%%%%%%%% calculate drug ROCs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ ATTrocs,ARMI, DrugROC, DrugMI, mean_act ] = jones_ROC_MI( sps,ndir, 1, t_win);
%%%% ATTrocs=zeros(2,3,6); (drug,time window (i.e. alignment),attention comparison)
%%%  ARMI=zeros(2,3,6); (drug,time window (i.e. alignment),attention comparison)
%%%% DrugMI=zeros(3,6); (alignment, ROCs no drug-drug for the 6 stimulus/attention conditions)
%%%% DrugROC=zeros(3,6);(alignment, ROCs no drug-drug for the 6 stimulus/attention conditions)
[ Fanos ] = jones_Fano( sps,ndir, 1 ); %%% Fanos=zeros(2,3,6):= Fanos(drug, time window (i.e. alignment),attention/stimulus comparison, i.e. 6 conditions)


if write_histos==1
    currdir=cd;
    jones_write_histo( meandata, t_stim, P, 'd:\data', drugtext);
end
currdir=cd;
penloc=findstr(NLX_Path,'pen');
penloc=NLX_Path(penloc:penloc+5);
cell_text=sprintf('%s %s%d %d',penloc, NAME{1}, exten, CELL(1));

if write_data==1 
    jones_write_data(ATTrocs,ARMI, DrugROC, DrugMI,mean_act, Fanos, P, 'd:\data', cell_text, drugtext, unit_identifier, corr_prctile, ndir);
end
if write_Tres_FF==1
    jones_write_TRes_FF( FF_Tres,P, 'd:\data', cell_text, drugtext, unit_identifier, corr_prctile, ndir,mean_act);
end
if plot_rasters==1
    pause
end
if write_eps
    cell_text=sprintf('%s_%s%d_Cell%d',penloc, NAME{1}, exten, CELL(1));
    figure(fighandle1);
    evaltext=sprintf('%s_Cycles%d_%d',cell_text, cycles(1), cycles(2));
    print('-depsc2','-tiff','-r300',evaltext);
end
eval(evaldir);

