function[vals]=jack_2target_psycho_EV_v3_mex_2SFs(file_of_int,sampleContrast,testContrast,session,conditions,area,analysisFolderAppend,monkey,roving)
%%jack_2target_psycho_EV_v3_mex_2SFs('24753928.1',30,[5 10 15 20 22 25 28
%%32 35 40 45 50 60 90],191,1:28,[],'jack',0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%This version is used on dual-quad core processor, calls the mex version of
%Nlx2MatEV_v4, Mat2NlxEV.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Written on 13/02/12 by Xing, modified from blanco's version:
%blanco_2target_psycho_EV_v3_mex
%Version 1 was written by Xing on 18/03/10, modified from previous function,
%'blanco_2target_psycho_NLX' dating from 25/01/10, which was modified from
%previous function, 'blanco_2target_psycho' dating from 30/10/09.
%Version 2 written on 22/03/10 to read NLX encoding of NumRandRew, an additional
%parameter which occurs during rewarded trials only. In version 2,
%length(encode_arr)==39 for correct saccade trials, but
%length(encode_arr)==37 for distractor saccade trials.
%In version 1, length(encode_arr) was always 37.
%Version 3 takes a 4th input arg, 'session,' and prints it on the graphs
%generated.


%Usage:
%NAME= Filename to analyse
%For analysing performance with presentation of a fixation point.
%Timing files used for Cortex programmes include gabornlx.stt and any of
%the perceptual learning contrast discrimination tasks, using gabor stimuli,
%e.g. pl23.stt and easy.stt.
%Reads data file from training session and calculates performance
%over time, for each condition.
%Use the variable 'sizebin' to adjust number of trials into which entire
%session is divided, eg. size of 500 gives ~2000/500= ~4 quarters.
%Set sizebin=trial to calculate performance averaged across all trials
%across 1 session.

% % encode(10001+cond_no)
% % encode(25000+sample_time)
% % encode(26000+response_time)
% % encode(START_EYE_DATA) 		100: put_eye_data_in_buf(ON)
% % encode(TURN_FIXSPOT_ON)     fix spot on
% % encode(STIM1_ON)            sample on
% % encode(STIM1_OFF)           sample off
% % encode(STIM2_ON)            test on
% % encode(STIM2_OFF)           test off
% % encode(STIM3_ON)            targets on
% % encode(FIX_COL_CHANGE)      fix spot colour change
% % encode(FIXATION_OCCURS)     fixation begins
% % encode(SACCADE_CORR)        correct saccade
% % encode(FIX_BREAK)           fixation break
% % encode(SACCADE_DIST)        incorrect saccade to distractor
% % encode(SACCADE_TEST)        incorrect saccade to test
% % encode(END_EYE_DATA)        101
% % encode(320+num_rand_rew)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%codes listed in CX_CODES.STT:
PARAMBASE              =10000;
TRIALPARAM_START       =300;
TRIALPARAM_END         =301;
STIMPARAM_START        =302;
STIMPARAM_END          =303;
STIM_SWITCH            =304;
REWARDPARAM_START      =305;
REWARDPARAM_END        =306;
STIM1_ON               =307; %sample
STIM1_OFF              =308;
STIM2_ON               =309; %test
STIM2_OFF              =310;
STIM3_ON               =311; %targets
FIX_COL_CHANGE         =312;
SACCADE_CORR           =313;
SACCADE_DIST           =314;
SACCADE_TEST           =315;
FIX_BREAK              =316;
REW_NUM_BASE           =320;
SACCADE_ONSET          =5006;

NLX_TRIAL_START      =255;
NLX_RECORD_START      =2;
NLX_SUBJECT_START     =4;
NLX_STIM_ON           =8;
NLX_STIM_OFF          =16;
NLX_TEST_ON           =102;
NLX_TEST_OFF          =103;
NLX_SUBJECT_END       =32;
NLX_RECORD_END        =64;
NLX_TRIAL_END        =254;
NLX_FIX_BREAK         =39;
NLX_TARGET_ON         =40;
NLX_TARGET_OFF        =41;
NLX_FIX_COL_CHANGE        =42;
NLX_SACCADE_CORR          =43;
NLX_SACCADE_ERROR         =44;
NLX_SACCADE_TEST          =45;
NLX_REW_NUM               =46;
NLX_FIX_SPOT_ON           =59;
NLX_FIX_START             =60;

NLX_TESTDIMMED        =17;
NLX_DISTDIMMED        =18;
NLX_BARRELEASED       =19;
NLX_CUE_ON             =20;
NLX_CUE_OFF        	=21;
NLX_EVENT_1        =9;
NLX_EVENT_2        =10;
NLX_EVENT_3        =11;
NLX_EVENT_4        =12;
NLX_EVENT_5        =13;
NLX_EVENT_6        =14;
NLX_EVENT_7        =15;

NLX_READ_DATA        =128;

%condition parameter encodes (send as 1 bytes)
NLX_TRIALPARAM_START =253;
NLX_TRIALPARAM_END   =252;
%Stimparameter encodes (send as 2 bytes)
NLX_STIMPARAM_START =251;
NLX_STIMPARAM_END   =250;

%codes listed in ENCODES.H:
REWARD					 					=3;
FIXATION_OCCURS					 		    =8;
START_INTER_TRIAL					 		=9;
END_INTER_TRIAL					 		    =10;
START_WAIT_FIXATION					     	=11;
END_WAIT_FIXATION					 		=12;
START_PRE_TRIAL					 	    	=15;
END_PRE_TRIAL					 			=16;
START_POST_TRIAL					 		=17;
END_POST_TRIAL					 			=18;
TURN_FIXSPOT_ON					 		    =35;
START_EYE_DATA								=100;
END_EYE_DATA								=101;

% [time_arr,event_arr,eog_arr,header,trial]=readcort(NAME);
% sampleContrast=30;
% testContrast=[10 15 20 25 27 28 29 31 32 33 35 40 50 60];
% testContrast=[10 15 20 25 27 29 31 33 35 40 50 60];
numconds=length(testContrast);
SFs=[4 2];
SFCol='kb';
calculateTangent=1;

% conditions=1:numconds;
NAME=file_of_int;

[TimeStamps, EventIDs, Nttls, Extras, EventStrings, Header] =Nlx2MatEV('Events.Nev', [1 1 1 1 1],1,1,1);
onfile=sprintf('%s on',file_of_int);
onfile2=sprintf('on');
on_time=[];
off_time=[];
offfile=sprintf('%s off',file_of_int);
for j=1:length(EventStrings)
    prl=EventStrings(j);
    if strcmp(prl,onfile)%look for keyboard entry which notes start of Cortex file writing
        on_time=TimeStamps(j);
    end
    if strcmp(prl,offfile)%look for keyboard entry which notes end of Cortex file writing
        off_time=TimeStamps(j);
    end
    k=findstr(onfile2, char(prl));
    if ~isempty(k)
        prl;
    end
end
if isempty(on_time)%first value in TimeStamps
    sprintf('exact match for ON time not found,press any key to continue')
    %     pause
    for j=1:length(EventStrings)
        prl=EventStrings(j);
        a=uint8(onfile);b=uint8(char(prl));
        if a(1:3)==b(1:3)
            if a(end-2:end)==b(end-2:end)%long, annoying way of comparing first and last 3 characters
                on_time=TimeStamps(j);
                sprintf('rough match for ON time found,press any key to continue')
                %                 pause
            end
        end
    end
end
if isempty(off_time)
    sprintf('exact match for OFF time not found,press any key to continue')
    %     pause
    for j=1:length(EventStrings)
        prl=EventStrings(j);
        a=uint8(offfile);b=uint8(char(prl));
        if a(1:3)==b(1:3)
            if a(end-2:end)==b(end-2:end)%long, annoying way of comparing first and last 3 characters
                off_time=TimeStamps(j);
                sprintf('rough match for OFF time found,press any key to continue')
                %                 pause
            end
        end
    end
end
if isempty(on_time)%first value in TimeStamps
    on_time=TimeStamps(1);
end
if isempty(off_time)
    off_time=TimeStamps(end);%last value in TimeStamps
end

if session==415
    on_time=TimeStamps(1);
    off_time=TimeStamps(end);
end

on_time
off_time
prl2=find(TimeStamps>on_time & TimeStamps<off_time);%from 2nd to second-last value in TimeStamps
timestamps2=TimeStamps(prl2);
Nttls2=Nttls(prl2);
event_arr=zeros(100,2000);%arbitrarily large empty matrix
time_arr=zeros(100,2000);
count_resp_trials=0;
max_event_length=0;
corr_resp=zeros(1,numconds);
error_resp=zeros(1,numconds);
prl=find(Nttls2==255);%#define NLX_TRIAL_START (Alex previously set to 253, which was #define NLX_TRIALPARAM_START 253)
for SFInd=1:length(SFs)
    for k=1:length(prl)-1%look through all trials
        for j=1:numconds
            if Nttls2(prl(k)+6)==conditions(j)+(numconds*(SFInd-1))%condition number (Alex previously set to 4, as count started from NLX_TRIALPARAM_START instead of NLX_TRIAL_START
                encode_arr=Nttls2(prl(k):prl(k+1));%first value is 255 for start of current trial, ends with 255 of next trial
                corr_resp(1,j)=corr_resp(1,j)+length(find(encode_arr==43));%#define NLX_SACCADE_CORR
                error_resp(1,j)=error_resp(1,j)+length(find(encode_arr==44));%#define NLX_SACCADE_ERROR
                if ~isempty(find(encode_arr==43))&&length(encode_arr)==39||~isempty(find(encode_arr==44))&&length(encode_arr)==37||~isempty(find(encode_arr==43))&&length(encode_arr)==37
                    count_resp_trials=count_resp_trials+1;
                    if length(encode_arr)>max_event_length
                        max_event_length=length(encode_arr);
                    end
                    if ~(length(encode_arr)==37||length(encode_arr)==39)
                        length(encode_arr)
                    end
                    event_arr(1:length(encode_arr)+1,count_resp_trials)=[k encode_arr]';%each trial occupies 1 column, trial number is value in first row
                    time_arr(1:length(encode_arr)+1,count_resp_trials)=[k timestamps2(prl(k):prl(k+1))]';
                end
            end
        end
    end
end
event_arr=event_arr(1:max_event_length,1:count_resp_trials);
time_arr=time_arr(1:max_event_length,1:count_resp_trials);
% prop_corr(1:numconds/2)=error_resp(1:numconds/2)./(error_resp(1:numconds/2)+corr_resp(1:numconds/2));
% prop_corr((numconds/2)+1:numconds)=1-(error_resp((numconds/2)+1:numconds)./(error_resp((numconds/2)+1:numconds)+corr_resp((numconds/2)+1:numconds)));
% set(figure,'Name','whole session');
% plot(testContrast, prop_corr,'ok');
% hold on
% X0=[40 2];
% X= fminsearch('weibull_zero_one',X0,[],testContrast,prop_corr);
% xvals=0:1:65;
% yvals=1-exp(-((xvals/X(1)).^X(2)));
% plot(xvals,yvals,'r');
% hold off

% fig1 =  figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.2, 0.3, 0.5, 0.5]);
% set(fig1, 'NumberTitle', 'off', 'Name', 'fig1');
% set(fig1, 'PaperUnits', 'centimeters', 'PaperType', 'usletter', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 26.65 20.305]);
% clf;
trial=size(event_arr,2);
tic

vals=zeros(trial,40)-1;
%columns from 1 to 32:
% 1	NLX_TRIAL_START
% 2	NLX_TRIALPARAM_START
% 3	BlockNr
% 4	CondNr
% 5	NLX_TRIALPARAM_END
% 6	NLX_RECORD_START
% 7	NLX_FIX_SPOT_ON
% 8	NLX_FIX_START
% 9	NLX_FIX_BREAK
% 10	NLX_STIM_ON
% 11	NLX_STIM_OFF
% 12	NLX_TEST_ON
% 13	NLX_TEST_OFF
% 14	NLX_TARGET_ON
% 15	NLX_FIX_COL_CHANGE
% 16	NLX_SACCADE_CORR
% 17	NLX_SACCADE_ERROR
% 18	NLX_SACCADE_TEST
% 19	NLX_RECORD_END
% 20	NLX_REW_NUM
% 21	NLX_TRIAL_END
% 22	RT saccade to correct target (calculated in Matlab)
% 23	RT saccade to distractor (calculated in Matlab)
% 24	time between fixation and early fix break (calculated in Matlab)
% 25	fixation held continuously (no response- calculated in Matlab)
% 26	fix1
% 27	sample duration
% 28	fix2
% 29	test duration
% 30	fix3
% 31	interval between targets onset & fixspot colour change
% 32	trial duration


sizebin=trial;
% sizebin=floor(trial/4);
bins=zeros(1,floor(trial/sizebin));
for j=1:length(bins)
    bins(j)=j*sizebin;
end;
for SFInd=1:length(SFs)
    correct=zeros(length(testContrast),1);
    incorrect=zeros(length(testContrast),1);
    RT=[];
    RT_distract=[];
    for k=1:length(bins)%eg.2000/50=40 bins of 50 trials, in 1 session
        trial=bins(k);
        for h=1:numconds
            for i=trial-sizebin+1:trial
                %if (header(13,i)==0)
                if event_arr(8,i)==conditions(h)+(numconds*(SFInd-1))% condition number, test has lower contrast than sample
                    vals(i,4)=event_arr(8,i);
                    %             elseif ~isempty(find(event_arr(:,i)==10001+5+h))% condition
                    %             number, test has higher contrast than sample
                    %                 vals(i,1)=h;%this collapses across direction of stim colour change- can do h+5 if want to keep separate
                end
                if vals(i,4)==conditions(h)+(numconds*(SFInd-1))%||(conds(i)==h+5)
                    x1=find(event_arr(2:end,i)==START_PRE_TRIAL)+1;%deal with this! not included in event_arr at present
                    if ~isempty(x1)
                        %                     vals(i,2)=time_arr(x1,i);
                    end
                    x1=find(event_arr(2:end,i)==NLX_TRIAL_START)+1;
                    if ~isempty(x1)
                        vals(i,1)=time_arr(x1(1),i);%value in first cell of x1, as 255 is encoded two times- once for present trial, once for next
                    end
                    x1=find(event_arr(2:end,i)==NLX_TRIALPARAM_START)+1;%add one, because find returns value numbering from 2nd cell in event_rr, not from first
                    x2=find(event_arr(2:end,i)==NLX_TRIALPARAM_END)+1;
                    if ~isempty(x1)&&~isempty(x2)
                        vals(i,2)=time_arr(x1,i);%TRIALPARAM_START
                        vals(i,5)=time_arr(x2,i);%TRIALPARAM_END
                        x3=event_arr(3,i);
                        vals(i,3)=x3;
                        %                     x3=(find(event_arr(2:end,i)>=PARAMBASE));
                        %                     if length(x3)==6
                        %                         vals(i,5:10)=(x3-PARAMBASE)/100;%coordinates for fixspot, RF, size of fix window
                        %                     end
                    end
                    x1=find(event_arr(14:end,i)==NLX_RECORD_START)+13;
                    if ~isempty(x1)
                        vals(i,6)=time_arr(x1,i);
                    end
                    x1=find(event_arr(15:end,i)==NLX_FIX_SPOT_ON)+14; % find fix point onset
                    if ~isempty(x1)
                        vals(i,7)=time_arr(x1,i);
                    end
                    x1=find(event_arr(15:end,i)==NLX_FIX_START)+14; % find fixation start (behavioural)
                    if ~isempty(x1)
                        vals(i,8)=time_arr(x1,i);
                    end
                    x1=find(event_arr(15:end,i)==NLX_FIX_BREAK)+14; % find time of early fix break
                    if ~isempty(x1)
                        vals(i,9)=time_arr(x1,i);
                    end
                    x1=find(event_arr(15:end,i)==NLX_STIM_ON)+14; % find time of sample onset
                    if ~isempty(x1)
                        vals(i,10)=time_arr(x1,i);
                    end
                    x1=find(event_arr(15:end,i)==NLX_STIM_OFF)+14; % find time of sample offset
                    if ~isempty(x1)
                        vals(i,11)=time_arr(x1,i);
                    end
                    x1=find(event_arr(15:end,i)==NLX_TEST_ON)+14; % find time of test onset
                    if ~isempty(x1)
                        vals(i,12)=time_arr(x1,i);
                    end
                    x1=find(event_arr(15:end,i)==NLX_TEST_OFF)+14; % find time of test offset
                    if ~isempty(x1)
                        vals(i,13)=time_arr(x1,i);
                    end
                    x1=find(event_arr(15:end,i)==NLX_TARGET_ON)+14; % find time of target(s) onset
                    if ~isempty(x1)
                        vals(i,14)=time_arr(x1,i);
                    end
                    x1=find(event_arr(15:end,i)==NLX_FIX_COL_CHANGE)+14; % find time of fix spot colour change
                    if ~isempty(x1)
                        vals(i,15)=time_arr(x1,i);
                    end
                    x1=find(event_arr(15:end,i)==NLX_SACCADE_CORR)+14; % find time of saccade to correct target
                    if ~isempty(x1)
                        vals(i,16)=time_arr(x1,i);
                        correct(h)=correct(h)+1;
                        RT=[RT vals(i,16)-vals(i,14)];%from target onset to correct saccade
                    end
                    x1=find(event_arr(15:end,i)==NLX_SACCADE_ERROR)+14; % find time of saccade to incorrect target
                    if ~isempty(x1)
                        vals(i,17)=time_arr(x1,i);
                        incorrect(h)=incorrect(h)+1;
                        RT_distract=[RT_distract vals(i,17)-vals(i,14)];%from target onset to saccade to distractor
                    end
                    x1=find(event_arr(15:end,i)==NLX_SACCADE_TEST)+14; % find time of incorrect saccade to test
                    if ~isempty(x1)
                        vals(i,18)=time_arr(x1,i);
                    end
                    x1=find(event_arr(15:end,i)==NLX_RECORD_END)+14;
                    if ~isempty(x1)
                        vals(i,19)=time_arr(x1,i);
                    end
                    x1=find(event_arr(15:end,i)>=65)+14;
                    if ~isempty(x1)
                        y1=find(event_arr(x1,i)<84);%allows up to 20 reward pulses
                        if ~isempty(y1)
                            vals(i,20)=event_arr(x1(y1),i)-65;%number of reward pulses
                        end
                    end
                    x1=find(event_arr(15:end,i)==NLX_TRIAL_END)+14;
                    if ~isempty(x1)
                        vals(i,21)=time_arr(x1,i);
                    end
                end
            end
        end
    end
    ave_RT=mean(RT);
    std_RT=std(RT);
    ave_RTerror=mean(RT_distract);
    std_RTerror=std(RT_distract);
    [a b ci stats]=ttest2(RT,RT_distract)
    perf_bins=correct./(correct+incorrect);
    valsSF{SFInd}=vals;
    lowerInd=find(testContrast<sampleContrast);
    lowerInd=lowerInd(end);
    report_higher_contrast=[1.-perf_bins(1:lowerInd,:);perf_bins(lowerInd+1:numconds,:)];
    prop_corr=report_higher_contrast';
    plot(testContrast,prop_corr(k,:),'o','Color',SFCol(SFInd));
    hold on
    X0=[2 30 0.2 0.1];
    X=fminsearch(@fit_weibull,X0,[],testContrast,prop_corr(k,:),[],'least_square',[0 0 0 0],[],[0 0 0 0],[]);
    allX(k,:)=X;
    if calculateTangent==1
        allX(k,1)=100*X(1)*X(3)*exp(-(sampleContrast/X(2))^X(1) )*sampleContrast^(X(1)-1)*(1/X(2))^X(1);%multiply by 100 as prop_corr given as fraction, not percentage
    end
    %     X=fminsearch('weibull_zero_one',X0,[],testContrast,prop_corr)
    xvals=0:1:testContrast(end)+10;
    yvals=1-X(4)-X(3).*exp(-(xvals./X(2)).^X(1));
    %     yvals=1-exp(-((xvals/X(1)).^X(2)));
    PSE(k)=X(2).*(-log((0.5-X(4))/X(3))).^(1/X(1));
    allX(k,2)=PSE(k);
    allX(k,:)
    plot(xvals,yvals,'Color',SFCol(SFInd));%,'LineStyle','--');
    set(gca,'FontSize',[6],'YLim',[0,1.01],'XLim',[0,testContrast(end)+10],'YTickMode','manual');%'YTick',[0.1]
    line(PSE(k),0:0.01:1,'Color',SFCol(SFInd));    
end

% figure(fig1);
% for i=1:2
%     i=1;
%     test=subplot(2,1,i);
%     test=subplot(1,1,i);
%     %if i==1
%     titletext=sprintf('%s performance',NAME);
%     set(fig1,'NumberTitle','off','Name',titletext);
%     %end;
%     x=1.5*trial/20:1:2.5*trial/20;
%     cond=[mean(perf_bins(1:numconds/2,:));mean(perf_bins(numconds/2+1:numconds,:))]
%     for j=1:2
%         line(bins,cond(j,:),'Color',[0.5*j 0.5*(2-j) 0.5*j]);hold on
%         line(x,j*0.02,'Color',[0.5*j 0.5*(2-j) 0.5*j]);hold on
%         ptext1=sprintf('%d',j+(j-1)*5);
%         text('Position',[3*trial/20 j*0.02],'FontSize',[7],'String',ptext1);
%     end
%
%     % for j=1:numconds
%     %     line(bins,perf_bins(j,:),'Color',[0.1*(numconds-j) 0.1*j 0.1*j]);hold on
%     %     line(x,j*0.02,'Color',[0.1*(numconds-j) 0.1*j 0.1*j]);hold on
%     %     ptext1=sprintf('%d',j);
%     %     text('Position',[3*trial/20 j*0.02],'FontSize',[7],'String',ptext1);
%     % end
% end;

% pause
close all

