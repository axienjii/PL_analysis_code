function check_NCS_across_chs
%Written by Xing 3/4/2014, to compare .ncs data across channels and verify
%whether voltage values are incorrectly duplicated across channels.

%Set the following variables (monkey, area, session):
monkey = 'blanco';
monkey = 'jack';
area = 'v1';
session=121;
session=123;
session=122;
session=124;
session=125;

%No need to adjust the rest:
FieldSelection(1) = 1;
FieldSelection(2) = 1;
FieldSelection(3) = 1;
FieldSelection(4) = 1;
FieldSelection(5) = 1;
[TimeStamps, ChNumbers, SampFreq, NumValSamp, Samp, Header] = Nlx2MatCSC(['I:\pl_spk_artifactRemoved\jack\MehdiCheetahData\',num2str(session),'\SPK_NCS_artifactRemoved\spk_CSC7_artifactRemoved.ncs'], FieldSelection,  1, 1);
[TimeStamps2, ChNumbers2, SampFreq2, NumValSamp2, Samp2, Header2] = Nlx2MatCSC(['I:\pl_spk_artifactRemoved\jack\MehdiCheetahData\',num2str(session),'\SPK_NCS_artifactRemoved\spk_CSC9_artifactRemoved.ncs'], FieldSelection,  1, 1);
% [TimeStamps, ChNumbers, SampFreq, NumValSamp, Samp, Header] = Nlx2MatCSC(['I:\jack\MehdiCheetahData\',num2str(session),'\spk_CSC7.ncs'], FieldSelection,  1, 1);
% [TimeStamps2, ChNumbers2, SampFreq2, NumValSamp2, Samp2, Header2] = Nlx2MatCSC(['I:\jack\MehdiCheetahData\',num2str(session),'\spk_CSC9.ncs'], FieldSelection,  1, 1);
conditions=[1 2 3 4];
rightCueOnly=0;
plotAcrossConds=0;%draw rasters in same subplot regardless of condition, sum PSTH values across conditions
%bottom-up passive fixation in V4 or V1:
%for blanco: channels=[1 2 3 4 7 12 13 14 15 18 20 22 24 33 34 36 37 38 40 42 49 50 51 52 53 54 55 57 58 59 60];
cdText=['cd ','V:\thielelab\Groups\ThieleGroup\monkey_data\Jack\_jackgrid\j_bu_events_files\',num2str(session)];
eval(cdText);
if strcmp(area,'v4')&&strcmp(monkey,'jack')
channels=[1 2 3 4 5 6 8 10 24 35 37 39 40 41 49 50 52 53 54 56];%V4
switch(session)
case 121
file_name='bu_pasv4.1';
folder='G:\Jack\2012-08-28_14-33-11_j121';
[vals durations]=jack_2target_psycho_EV_v3_mex_bu_pas(file_name,session,conditions,area,monkey);
case 123
rightCueOnly=1;
file_name='bu_rfrig.1';
folder='G:\Jack\2012-09-10_15-23-02_j123';
conditions=[1 2 3 4 5 6];
[vals durations]=jack_2target_psycho_EV_v3_mex_rightcue_pas(file_name,session,conditions,area,monkey);
case 124
file_name='bu_jack.1';
folder='G:\Jack\2012-10-04_15-24-45_j124';
conditions=[1 2 3 4 5 6 7 8];
[vals durations]=jack_2target_psycho_EV_v3_mex_bu_sac(file_name,session,conditions,area,monkey);
end
elseif strcmp(area,'v1')&&strcmp(monkey,'jack')
channels=[7 9 11 12 13 14 15 16 17 18 19 20 21 22 23 25 26 27 28 29 30 31 32 51 55];%V1
switch(session)
case 122
file_name='bu_pasv1.1';
folder='G:\Jack\2012-09-04_10-46-11_j122';
[vals durations]=jack_2target_psycho_EV_v3_mex_bu_pas(file_name,session,conditions,area,monkey);
case 125
file_name='bu_jack.1';
folder='G:\Jack\2012-10-09_16-04-18_j125';
conditions=[1 2 3 4 5 6 7 8];
[vals durations]=jack_2target_psycho_EV_v3_mex_bu_sac(file_name,session,conditions,monkey);
%             [vals durations]=jack_2target_psycho_EV_v3_mex_bu_sac(file_name,session,conditions,area,monkey);
end
elseif strcmp(area,'v1')&&strcmp(monkey,'blanco')
channels=[8 9 10 11 15 17 19 21 23 25 26 27 28 29 31 39 44 45 46 47 48 61 62 63 64];%V1
session=378;
file_name='bu_sac.1';
folder='F:\bottom_up_analysis\blanco_378\2010-09-08_18-40-01_blanco_spike3DgoodCh';
conditions=[1 2 3 4 5 6 7 8];
[vals durations]=jack_2target_psycho_EV_v3_mex_bu_sac(file_name,session,conditions,area,monkey);
end
figure
for trialNum=1:size(vals,1)
    temp1=find(TimeStamps<vals(trialNum,15));
    temp2=find(TimeStamps(temp1)>=vals(trialNum,7));
    Samp(temp2)
    plot(1:length(Samp(temp2)),Samp(temp2),'k-');hold on
    temp1=find(TimeStamps2<vals(trialNum,15));
    temp2=find(TimeStamps2(temp1)>=vals(trialNum,7));
    Samp2(temp2)
    plot(1:length(Samp2(temp2)),Samp2(temp2),'r-');
end
figure
for trialNum=1:size(vals,1)
    temp1=find(TimeStamps<vals(1,15));%spontaneous activity
    temp2=find(TimeStamps(temp1)>=vals(1,10)-58800);
    Samp(temp2)
    plot(1:length(Samp(temp2)),Samp(temp2),'k-');hold on
    temp1=find(TimeStamps2<vals(1,15));%spontaneous activity
    temp2=find(TimeStamps2(temp1)>=vals(1,10)-58800);
    Samp2(temp2)
    plot(1:length(Samp2(temp2)),Samp2(temp2),'r-');
end

diff=TimeStamps2-TimeStamps;
record=[];
for count=1:length(TimeStamps)
    if diff(count)==0
        record=[record count];
    end
end