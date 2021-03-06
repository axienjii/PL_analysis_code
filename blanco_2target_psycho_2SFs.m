function[]=blanco_2target_psycho_2SFs(NAME)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Written by Xing on 08/05/13
%Usage:
%NAME= Filename to analyse
%For analysing performance with presentation of a fixation point.

%Reads data file from training session and calculates performance
%over time, for each condition.
%Use the variable 'sizebin' to adjust number of trials into which entire
%session is divided, eg. size of 500 gives ~2000/500= ~4 quarters.
%Set sizebin=trial to calculate performance averaged across all trials
%across 1 session.

% % encode(START_PRE_TRIAL)
% % encode(10001+cond_no)
% % encode(TRIALPARAM_START)
% % encode(PARAMBASE+((int) (fix_x*100.0)));
% % encode(PARAMBASE+((int) (fix_y*100.0)));
% % encode(PARAMBASE+((int) (rf_x*100.0)));
% % encode(PARAMBASE+((int) (rf_y*100.0)));
% % encode(PARAMBASE+ (int) (FixWinX*100));
% % encode(PARAMBASE+ (int) (FixWinY*100));
% % encode(PARAMBASE+((int) (SacWin*100.0)));
% % encode(PARAMBASE+ (int) (SF*100.0));
% % encode(PARAMBASE+ (int) (size*100.0));
% % encode(TRIALPARAM_END);
% % encode(END_PRE_TRIAL);
% % encode(START_EYE_DATA) 		100: put_eye_data_in_buf(ON)
% % encode(TURN_FIXSPOT_ON)     fix spot on
% % encode(FIXATION_OCCURS)     fixation begins
% % encode(STIM1_ON)            sample on
% % encode(STIM1_OFF)           sample off
% % encode(STIM2_ON)            test on
% % encode(STIM2_OFF)           test off
% % encode(STIM3_ON)            targets on
% % encode(FIX_COL_CHANGE)      fix spot colour change
% % encode(313)               correct saccade
% % encode()               fixation break
% % encode(314)               incorrect saccade to distractor
% % encode()               incorrect saccade to test
% % encode(END_EYE_DATA)        101
% % encode(+num_rand_rew)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% path(path,'H:\My Received Files\_work\cortex_files\alex_prg')
[time_arr,event_arr,eog_arr,header,trial]=readcort(NAME);

fig1 =  figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.2, 0.3, 0.5, 0.5]);
set(fig1, 'NumberTitle', 'off', 'Name', 'fig1');
set(fig1, 'PaperUnits', 'centimeters', 'PaperType', 'usletter', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 26.65 20.305]);
clf;

calculateTangent=1;

sampleContrast=30;
testContrast=[5 10 15 20 22 25 28 32 35 40 45 50 60 90];
savePsychoCurves=0;
trial
tic
SFs=[4 2];
SFCol='kb';
numconds=14;
rews=zeros(trial,1);
vals=zeros(trial,16)-1;%columns from 1 to 12:
%1  cond num 2 sample duration 3 max response dur 4 fix point onset 5 fixation begins
%6 sample onset 7 sample offset 8 test onset 9 test offset 10 target(s) onset
%11 fixspot colour change 12 correct target saccade 13 distractor saccade 14 test saccade
%15 early fix break 16 num of reward units

sizebin=trial;
% sizebin=80;%set value to 'trial' to calculate perf across entire session
%and generate graph
bins=zeros(1,floor(trial/sizebin));
for j=1:length(bins)
    bins(j)=j*sizebin;
end;
for SFInd=1:length(SFs)
    RT=[];
    RT_distract=[];
    perf_bins=zeros(numconds,length(bins));%max num of trials per condition manually set to 300- adjust as needed
    ndir=zeros(1,numconds);
    perf=zeros(numconds,4);
    trialnum=zeros(1,32);
    RT=[];RT_early=[];RT_distract=[];
    all_ave_RT=zeros(1,length(bins));
    for k=1:length(bins)%eg.2000/50=40 bins of 50 trials, in 1 session
        trial=bins(k);
        for h=1:numconds
            correct=0;
            wrong=0;
            early=0;
            held=0;heldcheck=0;
            for i=trial-sizebin+1:trial
                if event_arr(2,i)-10000==h+1+(numconds*(SFInd-1))
                    vals(i,1)=event_arr(2,i)-10000;%condition number
                    if event_arr(24,i)==313
                        vals(i,12)=1;%correct saccade
                        correct=correct+1;
                        RT=[RT time_arr(24,i)-time_arr(23,i)];
                    elseif event_arr(24,i)==314
                        vals(i,12)=0;%correct saccade
                        wrong=wrong+1;
                        RT_distract=[RT_distract time_arr(24,i)-time_arr(23,i)];
                    end
                end
            end
            perf(h,1)=correct;
            perf(h,2)=wrong;
            perf(h,3)=early;
            perf(h,4)=held;%change to parameter 'heldcheck' to check if it's same as with 'held'
        end
        ave_RT(SFInd)=mean(RT);
        std_RT(SFInd)=std(RT);
        ave_RTerror(SFInd)=mean(RT_distract);
        std_RTerror(SFInd)=std(RT_distract);
        allRT{SFInd}=RT;
        allRTerror{SFInd}=RT_distract;
        perf
        %perf(:,1:2)
        response_trials=zeros(numconds,1);
        percent_perf=zeros(numconds,1);
        for j=1:numconds
            response_trials(j)=perf(j,1)+perf(j,2);
            percent_perf(j)=perf(j,1)/response_trials(j);
        end
        %     index=find(RT>=180);%optional cut-off point: Cortex code should impose limit at 180 ms
        %     RT=RT(index);
        %     index=find(RT<=700);%optional cut-off point: Cortex code should impose limit at 700 ms
        %     RT=RT(index);
        %RT_early
        size(RT)
        max(RT)
        min(RT)
        ave_RT=mean(RT)
        all_ave_RT(k,SFInd)=ave_RT;
        perf_bins(1:numconds,k)=percent_perf(:,1);
    end
    perf_binsSF{SFInd}=perf_bins%performance for each condition in each time bin
    mean_perf_bins(SFInd)=mean(perf_bins,1)%average across conditions
    all_ave_RT
    
    figure(fig1);
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
    if savePsychoCurves==1&&k==1
        filename1=[rootFolder,'\PL\psycho_data\',monkey,'\psycho_curves\',num2str(session),'_',num2str(sampleContrast)];
        printtext=sprintf('print -dpng %s',filename1);
        eval(printtext);
    end
    xlabel('Test contrast (%)');
    ylabel('Proportion of correct trials');
end
[a b ci stats]=ttest2(allRT{1},allRT{2})
save 'F:\PL\SFs_4_2_jack_control\RTs' allRT allRTerror
% close all

% load 'F:\PL\SFs_4_2_jack_control\RTs_191' allRT allRTerror
% allRTtemp=allRT;
% allRTerrortemp=allRTerror;
% load 'F:\PL\SFs_4_2_jack_control\RTs_190' allRT allRTerror
% allRTreshape=[allRTtemp{1} allRTtemp{2} allRT{1} allRT{2} allRTerrortemp{1} allRTerrortemp{2} allRTerror{1} allRTerror{2}];
% allce=[zeros(1,length(allRTtemp{1})+length(allRTtemp{2})+length(allRT{1})+length(allRT{2}))+1 zeros(1,length(allRTerrortemp{1})+length(allRTerrortemp{2})+length(allRTerror{1})+length(allRTerror{2}))+2];
% allsession=[zeros(1,length(allRTtemp{1})+length(allRTtemp{2})+length(allRTerrortemp{1})+length(allRTerrortemp{2}))+1 zeros(1,length(allRT{1})+length(allRT{2})+length(allRTerror{1})+length(allRTerror{2}))+2];
% allsf=[zeros(1,length(allRTtemp{1})+length(allRTerrortemp{1})+length(allRT{1})+length(allRTerror{1}))+1 zeros(1,length(allRTtemp{2})+length(allRTerrortemp{2})+length(allRT{2})+length(allRTerror{2}))+2];
% [p,table,stats]=anovan(allRTreshape',{allsf',allce',allsession'})
% allRTSF4=[allRTtemp{1} allRT{1} allRTerrortemp{1} allRTerror{1}];
% mean(allRTSF4)
% std(allRTSF4)
% allRTSF2=[allRTtemp{2} allRT{2} allRTerrortemp{2} allRTerror{2}];
% mean(allRTSF2)
% std(allRTSF2)
