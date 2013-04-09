function vals=analyseVals2(vals,testContrast,conditions,sampleContrast,monkey,area,analysisFolderAppend,session,roving,onExternalHD)

if onExternalHD==1
    rootFolder='G:\PL_backup_060413';
else
    rootFolder='F:';
end
writeWeibullPsycho=0;%remember to set correct directory!
savePsychoCurves=0;
writeMeanPerf=0;
writeNumTrials=1;
calculateTangent=1;
numconds=length(testContrast);
%     perf_bins=zeros(numconds,length(bins));%max num of trials per condition manually set to 300- adjust as needed
perf=zeros(numconds,4);
RT=[];RT_early=[];RT_distract=[];corrTrialList=[];RT_conds=[];RTerror_conds=[];RTce=[];
%     all_ave_RT=zeros(1,length(bins));
fix1_dur=[];fix2_dur=[];fix3_dur=[];sample_dur=[];test_dur=[];fixchange_dur=[];trial_dur=[];
trial=size(vals,1);
sessSegments=3;%look at all trials across session
for k=1:sessSegments;%examine all trials, then just first and last 30% of trials
    if k==1
        sizebin=trial;%trial
        trialsRange=1:size(vals,1);
    elseif k==2
        sizebin=floor(size(vals,1)*0.3);%trial
        trialsRange=1:sizebin;
    elseif k==3
        sizebin=floor(size(vals,1)*0.3);%trial
        trialsRange=size(vals,1)-sizebin+1:size(vals,1);
    end
    for h=1:numconds
        correct=0;
        wrong=0;
        early=0;
        held=0;heldcheck=0;
        RT_cond=[];
        RTerror_cond=[];
        for trial=trialsRange
            if vals(trial,4)==conditions(h)
                if vals(trial,16)>-1 %if correct saccade occurred
                    RT=[RT vals(trial,16)-vals(trial,14)];%from target onset to correct saccade
                    RT_cond=[RT_cond vals(trial,16)-vals(trial,14)];
                    RTce=[RTce vals(trial,16)-vals(trial,14)];%combine across correct and error trials
                    vals(trial,22)=vals(trial,16)-vals(trial,14);
                    correct=correct+1;
                    corrTrialList=[corrTrialList trial];%compile list of trials where correct saccade made, to check against Mehdi's
                end
                if vals(trial,17)>-1
                    RT_distract=[RT_distract vals(trial,17)-vals(trial,14)];%from target onset to saccade to distractor
                    RTerror_cond=[RTerror_cond vals(trial,17)-vals(trial,14)];
                    RTce=[RTce vals(trial,17)-vals(trial,14)];%combine across correct and error trials
                    vals(trial,23)=vals(trial,17)-vals(trial,14);
                    wrong=wrong+1;
                end
                if vals(trial,9)>-1
                    RT_early=[RT_early vals(trial,9)-vals(trial,8)];%from start of fixation (behavioural) to early fix break (reference differs from RT for correct responses)
                    vals(trial,24)=vals(trial,9)-vals(trial,8);
                    early=early+1;
                end
                if (vals(trial,8)~=-1)&&(vals(trial,9)==-1)&&(vals(trial,16)==-1)&&(vals(trial,17)==-1)&&(vals(trial,18)==-1)%fix begun but never broken
                    vals(trial,25)=1;
                    held=held+1;
                end
                if (vals(trial,22)>-1)%only examine times for correct trials
                    if vals(trial,10)>-1%sample was turned on after spontaneous period
                        fix1_dur=[fix1_dur vals(trial,10)-vals(trial,8)];%from start of fixation (behavioural)TO sample onset
                        vals(trial,26)=vals(trial,10)-vals(trial,8);
                    end
                    if vals(trial,11)>-1%calculate sample duration
                        sample_dur=[sample_dur vals(trial,11)-vals(trial,10)];
                        vals(trial,27)=vals(trial,11)-vals(trial,10);
                    end
                    if vals(trial,12)>-1%sample-test interval
                        fix2_dur=[fix2_dur vals(trial,12)-vals(trial,11)];
                        vals(trial,28)=vals(trial,12)-vals(trial,11);
                    end
                    if vals(trial,13)>-1%calculate test duration
                        test_dur=[test_dur vals(trial,13)-vals(trial,12)];
                        vals(trial,29)=vals(trial,13)-vals(trial,12);
                    end
                    if vals(trial,14)>-1%test-target interval
                        fix3_dur=[fix3_dur vals(trial,14)-vals(trial,13)];
                        vals(trial,30)=vals(trial,14)-vals(trial,13);
                    end
                    if vals(trial,15)>-1%interval between target onset and fixspot colour change
                        fixchange_dur=[fixchange_dur vals(trial,15)-vals(trial,14)];
                        vals(trial,31)=vals(trial,15)-vals(trial,14);
                    end
                    if vals(trial,21)>-1%interval between trial start and end (not including pre-trial period)
                        trial_dur=[trial_dur vals(trial,21)-vals(trial,1)];
                        vals(trial,32)=vals(trial,21)-vals(trial,1);
                    end
                end
            end
        end
        perf(h,1)=correct;
        perf(h,2)=wrong;
        perf(h,3)=early;
        perf(h,4)=held;
        RT_conds{h}=RT_cond;
        RTerror_conds{h}=RTerror_cond;
        mean_RT_conds(h,1)=mean(RT_cond)/1000;
        mean_RT_conds(h,2)=std(RT_cond)/1000;
        mean_RTerror_conds(h,1)=mean(RTerror_cond)/1000;
        mean_RTerror_conds(h,2)=std(RTerror_cond)/1000;
        if k==1
            perfWholeTrial=perf;
        end
    end
    perf
    rew_pulse_tally=[];
    for l=1:max(vals(:,20))
        rew_pulse_tally=[rew_pulse_tally length(find(vals(:,20)==l))];
    end
    rew_pulse_tally
    %perf(:,1:2)
    percent_perf=perf(:,1)./(perf(:,1)+perf(:,2));
    %     index=find(RT>=180);%optional cut-off point: Cortex code should impose limit at 180 ms
    %     RT=RT(index);
    %     index=find(RT<=700);%optional cut-off point: Cortex code should impose limit at 700 ms
    %     RT=RT(index);
    %RT_early
    size(RT);
    max(RT);
    min(RT);
    ave_RT=mean(RT);
    std_RT=std(RT);
    ave_RTerror=mean(RT_distract);
    std_RTerror=std(RT_distract);
    all_ave_RT(k)=ave_RT;
    all_std_RT(k)=std_RT;
    all_ave_RTerror(k)=ave_RTerror;
    all_std_RTerror(k)=std_RTerror;
    meanRTce=mean(RTce);
    stdRTce=std(RTce);
    perf_bins(1:numconds,k)=percent_perf;
    perf_bins%performance for each condition in each time bin
    mean_perf_bins=mean(perf_bins,1)%average across conditions
    lowerInd=find(testContrast<sampleContrast);
    lowerInd=lowerInd(end);
    report_higher_contrast=[1.-perf_bins(1:lowerInd,:);perf_bins(lowerInd+1:numconds,:)];
    prop_corr=report_higher_contrast';
    figureNames=[{'whole session'} {'first 30%'} {'last 30%'}];
    set(figure,'Name',figureNames{k});
    plot(testContrast,prop_corr(k,:),'ok');
    hold on
    X0=[2 30 0.2 0.1];
    X=fminsearch(@fit_weibull,X0,[],testContrast,prop_corr(k,:),[],'least_square',[0 0 0 0],[],[0 0 0 0],[])
    allX(k,:)=X;
    if calculateTangent==1
        allX(k,1)=100*X(1)*X(3)*exp(-(sampleContrast/X(2))^X(1) )*sampleContrast^(X(1)-1)*(1/X(2))^X(1);%multiply by 100 as prop_corr given as fraction, not percentage
    end
    %     X=fminsearch('weibull_zero_one',X0,[],testContrast,prop_corr)
    xvals=0:1:testContrast(end)+10;
    yvals=1-X(4)-X(3).*exp(-(xvals./X(2)).^X(1));
    %     yvals=1-exp(-((xvals/X(1)).^X(2)));
    PSE(k)=X(2).*(-log((0.5-X(4))/X(3))).^(1/X(1));
    plot(xvals,yvals,'r');
    set(gca,'FontSize',[6],'YLim',[0,1.01],'XLim',[0,testContrast(end)+10],'YTickMode','manual');%'YTick',[0.1]
    line(PSE(k),0:0.01:1,'Color','r');
    if savePsychoCurves==1&&k==1
        filename1=[rootFolder,'\PL\psycho_data\',monkey,'\psycho_curves\',num2str(session),'_',num2str(sampleContrast)];
        printtext=sprintf('print -dpng %s',filename1);
        eval(printtext);
    end
    % write Weibull function constants for psychometric curve to file:
    if writeWeibullPsycho==1
        currdir=cd
        if roving==1
            analysisFolderAppend=[analysisFolderAppend,'\sample_',num2str(sampleContrast)];
        end
        if strcmp(monkey,'jack')
            cdText=['cd ',rootFolder,'\jack\j_',area,'_roc_analysis',analysisFolderAppend];
        elseif strcmp(monkey,'blanco')
            cdText=['cd ',rootFolder,'\blanco\',area,'_roc_analysis',analysisFolderAppend];
        end
        eval(cdText)
        fid=fopen('psycho_constants','a+');
        fprintf(fid,'%d',session);
        for j=1:length(testContrast)
            fprintf(fid,' %f',prop_corr(j));
        end
        fprintf(fid,' %f %f',PSE,X(1),X(3),X(4));
        %         fprintf(fid,' %f %f',X(1),X(2));
        fprintf(fid,'\n');
        fclose(fid);
        chdirtext=sprintf('cd ''%s''',currdir);
        eval(chdirtext);
        cd
%     all_ave_RT(k)=ave_RT;
%     perf_bins(1:numconds,k)=percent_perf;
%     perf_bins%performance for each condition in each time bin
%     mean_perf_bins=mean(perf_bins,1)%average across conditions
%     all_ave_RT
    end
end
durations=zeros(7,2);
durations(1,1)=mean(all_ave_RT);
durations(2,1)=mean(fix1_dur);%ave_fix1
durations(3,1)=mean(fix2_dur);%ave_fix2
durations(4,1)=mean(fix3_dur);%ave_fix3
durations(5,1)=mean(sample_dur);%ave_sample_dur
durations(6,1)=mean(test_dur);%ave_test_dur
durations(7,1)=mean(trial_dur);%ave_trial_dur
durations(1,2)=std(RT);
durations(2,2)=std(fix1_dur);
durations(3,2)=std(fix2_dur);
durations(4,2)=std(fix3_dur);
durations(5,2)=std(sample_dur);
durations(6,2)=std(test_dur);
durations(7,2)=std(trial_dur);
durations=durations./1000;
round(durations)
    
% else for cond=1:size(perf_bins,2)
%         titletext=sprintf('bin %d',cond);
%         set(figure,'Name',titletext);
%         plot(testContrast,prop_corr(cond,:),'ok');
%         %plot(testContrast(1:numconds),report_higher_contrast(:,cond),'LineStyle', 'none','Marker','o','MarkerSize',4);hold on
%         line(sampleContrast,0:0.01:1);
%         hold on
%         %         X0=[40 2];
%         X0=[2 30 0 1];
%         %         X=fminsearch('weibull_zero_one',X0,[],testContrast,prop_corr(cond,:))%X(1) is mid point, X(2) is slope
%         X=fminsearch(@fit_weibull,X0,[],testContrast,prop_corr,[],'least_square',[0 0 0 0],[],[0 0 0 0],[])
%         xvals=0:1:testContrast(end)+10;
%         %         yvals=1-exp(-((xvals/X(1)).^X(2)));
%         yvals=1-X(4)-X(3).*exp(-(xvals./X(2)).^X(1));
%         plot(xvals,yvals,'r');
%         set(gca,'FontSize',[6],'YLim',[0,1.01],'XLim',[0,testContrast(end)+10],'YTickMode','manual');%'YTick',[0.1]
%     end

if writeMeanPerf==1
    matPath=[rootFolder,'\PL\psycho_data\',monkey,'\allMeanPerf_',area,'_',num2str(sampleContrast),'.mat'];
    if ~exist(matPath,'file')
        allMeanPerf=[];
    else
        loadText=['load ',matPath,' allMeanPerf'];
        eval(loadText)
    end
    if ~isempty(allMeanPerf)
        sessionInd=find(allMeanPerf(:,1)==session);
        if ~isempty(sessionInd)
            if sessSegments==1
                allMeanPerf(sessionInd,:)=[session mean_perf_bins(1) prop_corr(1,:) PSE(1) allX(1,1) allX(1,3) allX(1,4) all_ave_RT(1) all_std_RT(1) mean_RT_conds(:,1)' mean_RT_conds(:,2)' durations(:,1)' durations(:,2)' all_ave_RTerror(1) all_std_RTerror(1) mean_RTerror_conds(:,1)' mean_RTerror_conds(:,2)' meanRTce stdRTce];
            elseif sessSegments==3
                allMeanPerf(sessionInd,:)=[session mean_perf_bins(1) prop_corr(1,:) PSE(1) allX(1,1) allX(1,3) allX(1,4) all_ave_RT(1) all_std_RT(1) mean_RT_conds(:,1)' mean_RT_conds(:,2)' durations(:,1)' durations(:,2)' all_ave_RTerror(1) all_std_RTerror(1) mean_RTerror_conds(:,1)' mean_RTerror_conds(:,2)' meanRTce stdRTce mean_perf_bins(2) prop_corr(2,:) PSE(2) allX(2,1) allX(2,3) allX(2,4) all_ave_RT(2) all_std_RT(2) mean_perf_bins(3) prop_corr(3,:) PSE(3) allX(3,1) allX(3,3) allX(3,4) all_ave_RT(3) all_std_RT(3)];
            end
        else
            if sessSegments==1
                allMeanPerf=[allMeanPerf;session mean_perf_bins(1) prop_corr(1,:) PSE(1) allX(1,1) allX(1,3) allX(1,4) all_ave_RT(1) all_std_RT(1) mean_RT_conds(:,1)' mean_RT_conds(:,2)' durations(:,1)' durations(:,2)' all_ave_RTerror(1) all_std_RTerror(1) mean_RTerror_conds(:,1)' mean_RTerror_conds(:,2)' meanRTce stdRTce];
            elseif sessSegments==3
                allMeanPerf=[allMeanPerf;session mean_perf_bins(1) prop_corr(1,:) PSE(1) allX(1,1) allX(1,3) allX(1,4) all_ave_RT(1) all_std_RT(1) mean_RT_conds(:,1)' mean_RT_conds(:,2)' durations(:,1)' durations(:,2)' all_ave_RTerror(1) all_std_RTerror(1) mean_RTerror_conds(:,1)' mean_RTerror_conds(:,2)' meanRTce stdRTce mean_perf_bins(2) prop_corr(2,:) PSE(2) allX(2,1) allX(2,3) allX(2,4) all_ave_RT(2) all_std_RT(2) mean_perf_bins(3) prop_corr(3,:) PSE(3) allX(3,1) allX(3,3) allX(3,4) all_ave_RT(3) all_std_RT(3)];
            end
        end
    else
        allMeanPerf=[];
        if sessSegments==1
            allMeanPerf=[allMeanPerf;session mean_perf_bins(1) prop_corr(1,:) PSE(1) allX(1,1) allX(1,3) allX(1,4) all_ave_RT(1) all_std_RT(1) mean_RT_conds(:,1)' mean_RT_conds(:,2)' durations(:,1)' durations(:,2)' all_ave_RTerror(1) all_std_RTerror(1) mean_RTerror_conds(:,1)' mean_RTerror_conds(:,2)' meanRTce stdRTce];
        elseif sessSegments==3
            allMeanPerf=[allMeanPerf;session mean_perf_bins(1) prop_corr(1,:) PSE(1) allX(1,1) allX(1,3) allX(1,4) all_ave_RT(1) all_std_RT(1) mean_RT_conds(:,1)' mean_RT_conds(:,2)' durations(:,1)' durations(:,2)' all_ave_RTerror(1) all_std_RTerror(1) mean_RTerror_conds(:,1)' mean_RTerror_conds(:,2)' meanRTce stdRTce mean_perf_bins(2) prop_corr(2,:) PSE(2) allX(2,1) allX(2,3) allX(2,4) all_ave_RT(2) all_std_RT(2) mean_perf_bins(3) prop_corr(3,:) PSE(3) allX(3,1) allX(3,3) allX(3,4) all_ave_RT(3) all_std_RT(3)];
        end
    end
    if sessSegments==1
        segmentsText='all_early_late';
    else
        segmentsText=[];
    end
    saveText=['save ',rootFolder,'\PL\psycho_data\',monkey,'\allMeanPerf_',area,'_',num2str(sampleContrast),'.mat allMeanPerf'];
    eval(saveText)
end
if writeNumTrials==1
    matPath=[rootFolder,'\PL\psycho_data\',monkey,'\allNumTrials_',area,'_',num2str(sampleContrast),'.mat'];
    if ~exist(matPath,'file')
        perfAll=[];
    else
        loadText=['load ',rootFolder,'\PL\psycho_data\',monkey,'\allNumTrials_',area,'_',num2str(sampleContrast),'.mat perfAll'];
        eval(loadText)
    end
    if ~isempty(perfAll)
        sessionInd=find(perfAll(:,1)==session);
        if ~isempty(sessionInd)
            perfAll(sessionInd,:)=[session perfWholeTrial(:,1)' sum(perfWholeTrial(:,1))];
        else
            perfAll=[perfAll;session perfWholeTrial(:,1)' sum(perfWholeTrial(:,1))];
        end
    else
        perfAll=[];
        perfAll=[perfAll;session perfWholeTrial(:,1)' sum(perfWholeTrial(:,1))];
    end
    saveText=['save ',rootFolder,'\PL\psycho_data\',monkey,'\allNumTrials_',area,'_',num2str(sampleContrast),'.mat perfAll'];
    eval(saveText)
end
