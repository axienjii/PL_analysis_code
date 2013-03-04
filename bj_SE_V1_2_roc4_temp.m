function[rocvals]=bj_SE_V1_2_roc4(animal,channel,session,verbose)
%Modified 30/01/13 to process artifact-free, spontaneous-activity-matched data.
%Compares spike activity between sample
%and test presentations, for each condition, and calculates ROC value based
%on all trials for each condition. Plots ROC values against contrast using
%Weibull function.
%
%Analyses spike data (from 1 channel at a time) and matching event file.
%Input arguments are:
%1. 'name' (full path name for the spike data file, extension '.nse')
%2. 'file_of_int' (name of Cortex data file, no need for full path- used for
%identification of start/ end of Cortex file recording within Neuralynx events file)
%Must first change directory to folder containing events file (.nev).
%E.g. set directory to R:\monkey_data\blanco\_grid\320\2010-04-09_09-55-16_blanco_raw1
%and key in:
%blanco_SE_smallbins('R:\monkey_data\blanco\_grid\320\2010-04-09_17-29-43_b
%lanco_range\SpikeCh_12.NSE','21613614.2') at Matlab prompt
%Note that the events file is located in raw folder data whereas the spike
%files are in the playback folder- this is because manually-entered events
%are lost during playback processing.
%
%This function examines spike activity during correct trials, plots rasters
%and average response graphs for each condition.
%Note that as sample-test interval duration is set to random(512,1024), and
%as the exact timing of events has some between-trial variation anyway, this
%function analyses spike activity according to the actual timestamp values,
%and then aligns the activity across trials to a standardised timeline.
%
%Divides each trial into 5 periods:
%-------------------------------    --------
%Period                             Duration
%-------------------------------    --------
%1. spontaneous (fix1)              512 ms
%2. sample                          512 ms
%3. sample-test interval (fix2)     512 - 1024 ms
%4. test                            512 ms
%5. test-target interval (fix3)     400 ms
%
%Thus, rasters and average activity graphs obtain data from trials with
%temporally variable events, and align them according to an idealised time sequence.
%The third time period can be visualised in its entirety: 'set(test,'XLim',[512 1536]);'
%OR only for the first 512 ms: 'set(test,'XLim',[512 1024]);,' as desired.

artifactTrialsName=[num2str(session),'_corrtrialartifact.mat'];
artifactTrialsPath=fullfile('F:','PL','pl_corr_art_trials',animal,artifactTrialsName);
loadText=['load ',artifactTrialsPath,' rlist removeTrialsTimestamps'];
eval(loadText);%removeTrialsTimestamps column 1: NLX_TRIAL_START; column 21: NLX_TRIAL_END
% removeTrialsTimestamps=[0 0];
plotRedArtifacts=1;
if ~sum(session==[355.2 405.2 435.2])%for split sessions, run .1 and .2 at the same time- when .1 is processed.
    splitSess=1;
    [file_of_int,testContrasts,sampleContrasts,expt_type,rotated,area]=session_metadata(session,animal);
    if length(sampleContrasts)==36
        sampleContrasts=[30 20 40];
        testContrasts=[testContrasts(1:12);testContrasts(13:24);testContrasts(25:36)];
        allConditions=[13:24;1:12;25:36];
    elseif length(sampleContrasts)==14
        sampleContrasts=30;
        roving=0;
        allConditions=[1:12];
    end    
    
    nseName=['SpikeCh_',num2str(channel),'_.nse'];
    if sum(session==[355.1 405.1 435.1])%combine data from split sessions
        sessionName=[num2str(floor(session)),'.1'];
        nsePath=fullfile('I:','pl_spnorm_nse',animal,sessionName,nseName);
        [SE_TimeStamps,missing]=open_nse_file(nsePath);
        sessionName=[num2str(floor(session)),'.2'];
        nsePath=fullfile('I:','pl_spnorm_nse',animal,sessionName,nseName);
        [SE_TimeStamps2,missing2]=open_nse_file(nsePath);
        missing=max([missing missing2]);
        splitSess=2;
    else
        nsePath=fullfile('I:','pl_spnorm_nse',animal,num2str(session),nseName);
        [SE_TimeStamps,missing]=open_nse_file(nsePath);
    end
    if missing==0
        for k=1:length(sampleContrasts)
            sampleContrast=sampleContrasts(k);
            testContrast=testContrasts(k,:);
            conditions=allConditions(k,:);
            fig=figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.2, 0.04, 0.5, 0.9]);
            set(fig,'PaperUnits','centimeters','PaperType','A4','PaperOrientation', 'portrait', 'PaperPosition', [0.63452 0.63452 21 28.41]);
            numconds=length(testContrast);
            binwidth=10*1000;%10 ms
            lengthBinsOver53=0;lengthBinsOver100=0;
            minBins=zeros(1,5)+56;%arbitrarily large number
            mean_act=zeros(numconds,5);
            if session==451
                a=find(vals(:,:)>=SE_TimeStamps(1));
                vals=vals(a(1):end,:);
            end
            all_spontan=[];all_sample=[];all_fix2=[];all_test=[];all_fix3=[];
            cond_arr=[];
            numTrialsCond=zeros(1,numconds);
            histo1=zeros(numconds,100);histo2=histo1;histo3=zeros(numconds,200);histo4=histo1;histo5=histo1;
            rocvals=[];%for roc analysis
            for h=1:length(conditions)
                spontan=[];sample_act=[];fix2_act=[];test_act=[];fix3_act=[];
                for sessHalf=1:splitSess
                    if splitSess==2
                        sessionName=[num2str(floor(session)),'.',num2str(sessHalf)];
                        nsePath=fullfile('I:','pl_spnorm_nse',animal,sessionName,nseName);
                        [SE_TimeStamps,dummy]=open_nse_file(nsePath);
                        valsFileName=['vals_',num2str(floor(session)),'.',num2str(sessHalf),'.mat'];
                        valsFolder=fullfile('F:','PL','vals_perf',animal); %#ok<NASGU>
                        valsPath=fullfile('F:','PL','vals_perf',animal,valsFileName);
                        load(valsPath);
                    elseif splitSess==1
                        nsePath=fullfile('I:','pl_spnorm_nse',animal,num2str(session),nseName);
                        [SE_TimeStamps,dummy]=open_nse_file(nsePath);
                        valsFileName=['vals_',num2str(session),'.mat'];
                        valsFolder=fullfile('F:','PL','vals_perf',animal); %#ok<NASGU>
                        valsPath=fullfile('F:','PL','vals_perf',animal,valsFileName);
                        load(valsPath);
                    end
                    if length(sampleContrasts)==3
                        if sampleContrast==30
                            vals=vals{1};
                        elseif sampleContrast==20
                            vals=vals{2};
                        elseif sampleContrast==40
                            vals=vals{3};
                        end
                    end
                    for i=1:size(vals,1)
                        not_artifact_trial=isempty(find(vals(i,1)==removeTrialsTimestamps(:,1), 1))&&isempty(find(vals(i,21)==removeTrialsTimestamps(:,2), 1));
                        if not_artifact_trial
                            rasterColour='k';
                            plotTrial=1;
                        elseif ~not_artifact_trial%this trial contains movement-induced artifacts
                            if plotRedArtifacts==1
                                rasterColour='r';
                                plotTrial=1;
                            elseif plotRedArtifacts==0
                                plotTrial=0;
                            end
                        end
                        if plotTrial
                            if vals(i,4)==conditions(h)%check condition number
                                %     for j=1:size(vals,2)
                                %     SE_EV_TimeStamps(j)=find(SE_TimeStamps<vals(i,j+1)&&SE_TimeStamps>=vals(i,j));
                                %     spikes(j)=DataPoints(:,:,SE_EV_TimeStamps(j));
                                %     end
                                if (vals(i,22)>-1)
                                    numTrialsCond(1,h)=numTrialsCond(1,h)+1;%tally number of trials per condition
                                    temp1=find(SE_TimeStamps<vals(i,10));%spontaneous activity
                                    temp2=find(SE_TimeStamps(temp1)>=vals(i,8));
                                    bins=vals(i,8):binwidth:vals(i,10);
                                    if length(bins)<minBins(1)
                                        minBins(1)=length(bins);
                                    end
                                    if length(bins)>minBins(1)
                                        lengthBinsOver53=lengthBinsOver53+1;
                                    end
                                    spikeTimes=SE_TimeStamps(temp2)-vals(i,8)-512000;%calculate times relative to beginning of spontan period, convert to ms
                                    [N X]=hist(SE_TimeStamps(temp2),bins);%returns the distribution of timestamps
                                    %among bins with centers specified by bins. The first bin includes
                                    %data between -inf and the first center and the last bin
                                    %includes data between the last bin and inf.
                                    if size(N,2)==1
                                        N=N';
                                    end
                                    histo1(h,1:length(N))=histo1(h,1:length(N))+N;%tally spikes across trials for each bin
                                    assignedspikes(h,numTrialsCond(h),1:length(temp2))=spikeTimes;
                                    test=subplot(numconds,5,h*5-4);
                                    for m=1:length(spikeTimes)%raster plots
                                        plot([spikeTimes(m)/1000 spikeTimes(m)/1000],[100+numTrialsCond(h)-0.5 100+numTrialsCond(h)+0.5],rasterColour);hold on
                                    end
                                    if numTrialsCond(h)==1
                                        matarray(h,1,:)={spikeTimes/1000};
                                    else
                                        matarray(h,1,:)={[matarray{h,1}(:,:);{spikeTimes/1000}]};
                                    end
                                    spontan=[spontan length(temp2)*1000000/(vals(i,10)-vals(i,8))];%convert from spikes/microsecond to spikes/second
                                    temp1=find(SE_TimeStamps<vals(i,11));%activity during sample presentation
                                    temp2=find(SE_TimeStamps(temp1)>=vals(i,10));
                                    bins=vals(i,10):binwidth:vals(i,11);
                                    if length(bins)<minBins(2)
                                        minBins(2)=length(bins);
                                    end
                                    if length(bins)>minBins(2)
                                        lengthBinsOver53=lengthBinsOver53+1;
                                    end
                                    spikeTimes=SE_TimeStamps(temp2)-vals(i,10);%calculate times relative to sample onset, convert to ms, convert to standardised sample onset time
                                    [N X]=hist(SE_TimeStamps(temp2),bins);
                                    if size(N,2)==1
                                        N=N';
                                    end
                                    histo2(h,1:length(N))=histo2(h,1:length(N))+N;%tally spikes across trials for each bin
                                    assignedspikes(h,numTrialsCond(h),1:length(temp2))=spikeTimes;
                                    test=subplot(numconds,5,h*5-3);
                                    for m=1:length(spikeTimes)%raster plots
                                        plot([spikeTimes(m)/1000 spikeTimes(m)/1000],[100+numTrialsCond(h)-0.5 100+numTrialsCond(h)+0.5],rasterColour);hold on
                                    end
                                    if numTrialsCond(h)==1
                                        matarray(h,2,:)={spikeTimes/1000};
                                    else
                                        matarray(h,2,:)={[matarray{h,2}(:,:);{spikeTimes/1000}]};
                                    end
                                    sample_act=[sample_act length(temp2)*1000000/(vals(i,11)-vals(i,10))];
                                    temp1=find(SE_TimeStamps<vals(i,12));%activity during sample-test interval
                                    temp2=find(SE_TimeStamps(temp1)>=vals(i,11));
                                    bins=vals(i,11):binwidth:vals(i,12);
                                    if length(bins)<minBins(3)
                                        minBins(3)=length(bins);
                                    end
                                    if length(bins)>minBins(3)
                                        lengthBinsOver100=lengthBinsOver100+1;
                                    end
                                    spikeTimes=SE_TimeStamps(temp2)-vals(i,11)+512000;%calculate times relative to sample offset, convert to ms, convert to standardised sample offset time
                                    [N X]=hist(SE_TimeStamps(temp2),bins);
                                    if size(N,2)==1
                                        N=N';
                                    end
                                    histo3(h,1:length(N))=histo3(h,1:length(N))+N;%tally spikes across trials for each bin
                                    assignedspikes(h,numTrialsCond(h),1:length(temp2))=spikeTimes;
                                    test=subplot(numconds,5,h*5-2);
                                    for m=1:length(spikeTimes)%raster plots
                                        plot([spikeTimes(m)/1000 spikeTimes(m)/1000],[100+numTrialsCond(h)-0.5 100+numTrialsCond(h)+0.5],rasterColour);hold on
                                    end
                                    if numTrialsCond(h)==1
                                        matarray(h,3,:)={spikeTimes/1000};
                                    else
                                        matarray(h,3,:)={[matarray{h,3}(:,:);{spikeTimes/1000}]};
                                    end
                                    fix2_act=[fix2_act length(temp2)*1000000/(vals(i,12)-vals(i,11))];
                                    temp1=find(SE_TimeStamps<vals(i,13));%activity during test presentation
                                    temp2=find(SE_TimeStamps(temp1)>=vals(i,12));
                                    bins=vals(i,12):binwidth:vals(i,13);
                                    if length(bins)<minBins(4)
                                        minBins(4)=length(bins);
                                    end
                                    if length(bins)>minBins(4)
                                        lengthBinsOver53=lengthBinsOver53+1;
                                    end
                                    spikeTimes=SE_TimeStamps(temp2)-vals(i,12)+1024000;%calculate times relative to test onset, convert to ms, convert to standardised test onset (the sample-test interval is variable in reality)
                                    [N X]=hist(SE_TimeStamps(temp2),bins);
                                    if size(N,2)==1
                                        N=N';
                                    end
                                    histo4(h,1:length(N))=histo4(h,1:length(N))+N;%tally spikes across trials for each bin
                                    assignedspikes(h,numTrialsCond(h),1:length(temp2))=spikeTimes;
                                    test=subplot(numconds,5,h*5-1);
                                    for m=1:length(spikeTimes)%raster plots
                                        plot([spikeTimes(m)/1000 spikeTimes(m)/1000],[100+numTrialsCond(h)-0.5 100+numTrialsCond(h)+0.5],rasterColour);hold on
                                    end
                                    if numTrialsCond(h)==1
                                        matarray(h,4,:)={spikeTimes/1000};
                                    else
                                        matarray(h,4,:)={[matarray{h,4}(:,:);{spikeTimes/1000}]};
                                    end
                                    test_act=[test_act length(temp2)*1000000/(vals(i,13)-vals(i,12))];
                                    temp1=find(SE_TimeStamps<vals(i,14));%activity during test-targets interval
                                    temp2=find(SE_TimeStamps(temp1)>=vals(i,13));
                                    bins=vals(i,13):binwidth:vals(i,14);
                                    if length(bins)<minBins(5)
                                        minBins(5)=length(bins);
                                    end
                                    if length(bins)>minBins(5)
                                        lengthBinsOver53=lengthBinsOver53+1;
                                    end
                                    spikeTimes=SE_TimeStamps(temp2)-vals(i,13)+1536000;%calculate times relative to test offset, convert to ms, convert to standardised test offset
                                    [N X]=hist(SE_TimeStamps(temp2),bins);
                                    if size(N,2)==1
                                        N=N';
                                    end
                                    histo5(h,1:length(N))=histo5(h,1:length(N))+N;%tally spikes across trials for each bin
                                    assignedspikes(h,numTrialsCond(h),1:length(temp2))=spikeTimes;
                                    test=subplot(numconds,5,h*5);
                                    for m=1:length(spikeTimes)%raster plots
                                        plot([spikeTimes(m)/1000 spikeTimes(m)/1000],[100+numTrialsCond(h)-0.5 100+numTrialsCond(h)+0.5],rasterColour);hold on
                                    end
                                    if numTrialsCond(h)==1
                                        matarray(h,5,:)={spikeTimes/1000};
                                    else
                                        matarray(h,5,:)={[matarray{h,5}(:,:);{spikeTimes/1000}]};
                                    end
                                    fix3_act=[fix3_act length(temp2)*1000000/(vals(i,14)-vals(i,13))];
                                end
                            end
                        end
                    end
                end
                all_spontan=[all_spontan spontan];
                all_sample=[all_sample sample_act];
                all_fix2=[all_fix2 fix2_act];
                all_test=[all_test test_act];
                all_fix3=[all_fix3 fix3_act];
                cond_arr=[cond_arr zeros(1,length(spontan))+h];
                mean_act(h,1)=mean(spontan);
                mean_act(h,2)=mean(sample_act);
                mean_act(h,3)=mean(fix2_act);
                mean_act(h,4)=mean(test_act);
                mean_act(h,5)=mean(fix3_act);
                %roc analysis:
                [roc]=sglroc3(test_act(:)',sample_act(:)');
                rocvals=[rocvals roc];
            end
            %                 rocvals
            %                 size(rocvals)
            %             figure
            %             plot(testContrast,rocvals,'ok');
            
            %     figure
            %     plot(testContrast,rocvals,'ok');
            %     hold on
            %     X0=[40 2];
            %     X=fminsearch('weibull_zero_one',X0,[],testContrast,rocvals)
            %     xvals=0:1:testContrast(end)+10;
            %     yvals=1-exp(-((xvals/X(1)).^X(2)));
            %     plot(xvals,yvals,'r');
            %     line(sampleContrast,0:0.01:1);
            
            for h=1:numconds
                histo1(h,:)=histo1(h,:)*1000000/(binwidth*numTrialsCond(1,h));%average activity per ms. note:*10 is an arbitrary scaling factor
                histo2(h,:)=histo2(h,:)*1000000/(binwidth*numTrialsCond(1,h));%average activity per ms
                histo3(h,:)=histo3(h,:)*1000000/(binwidth*numTrialsCond(1,h));%average activity per ms
                histo4(h,:)=histo4(h,:)*1000000/(binwidth*numTrialsCond(1,h));%average activity per ms
                histo5(h,:)=histo5(h,:)*1000000/(binwidth*numTrialsCond(1,h));%average activity per ms
            end
            maxhisto=[max(histo1) max(histo2) max(histo3) max(histo4) max(histo5)];%normalise to highest firing rates
            maxhisto=max(maxhisto);
            histo1=histo1.*100./maxhisto;histo2=histo2.*100./maxhisto;histo3=histo3.*100./maxhisto;histo4=histo4.*100./maxhisto;histo5=histo5.*100./maxhisto;
            for h=1:numconds
                test=subplot(numconds,5,h*5-4);
                plot((-1*minBins(1)*binwidth:binwidth:-1*binwidth)./1000,histo1(h,1:minBins(1)),'k');hold on%convert x-axis values to ms
                set(test,'YLim',[0 100+1.1*max(numTrialsCond)]);set(test,'XLim',[-512 0]);set(test,'YTickLabel',[0 round(maxhisto)]);
                ylabel([num2str(testContrast(h)),'     ']);
                set(get(gca,'YLabel'),'Rotation',0.0,'FontSize',8);
                if h~=numconds(end)
                    set(test,'XTickLabel','');
                else
                    set(test,'XTick',[-512 0]);
                    set(test,'XTickLabel',[-512 0]);
                end
                if h==1
                    ptext=sprintf('%s   %s   sample %s%%',num2str(floor(session)),file_of_int,num2str(sampleContrast));
                    text('Position',[-600 100+2.6*max(numTrialsCond)],'FontSize',9,'String',ptext);
                end
                test=subplot(numconds,5,h*5-3);
                plot((0:binwidth:minBins(2)*(binwidth-1))./1000,histo2(h,1:minBins(2)),'k');hold on%convert x-axis values to ms
                set(test,'YLim',[0 100+1.1*max(numTrialsCond)]);set(test,'XLim',[0 512]);set(test,'YTickLabel','');
                if h~=numconds(end)
                    set(test,'XTickLabel','');
                else
                    set(test,'XTick',[0 512]);
                    set(test,'XTickLabel',[0 512]);
                end
                test=subplot(numconds,5,h*5-2);
                plot((512000:binwidth:512000+minBins(3)*(binwidth-1))./1000,histo3(h,1:minBins(3)),'k');hold on%convert x-axis values to ms
                set(test,'YLim',[0 100+1.1*max(numTrialsCond)]);set(test,'XLim',[512 1024]);set(test,'YTickLabel','');%set(test,'XLim',[512 1536]);
                if h~=numconds(end)
                    set(test,'XTickLabel','');
                else
                    set(test,'XTick',[512 1024]);
                    set(test,'XTickLabel',[512 1024]);
                end
                test=subplot(numconds,5,h*5-1);
                plot((1024000:binwidth:1024000+minBins(4)*(binwidth-1))./1000,histo4(h,1:minBins(4)),'k');hold on%convert x-axis values to ms
                set(test,'YLim',[0 100+1.1*max(numTrialsCond)]);set(test,'XLim',[1024 1536]);set(test,'YTickLabel','');
                if h~=numconds(end)
                    set(test,'XTickLabel','');
                else
                    set(test,'XTick',[1024 1536]);
                    set(test,'XTickLabel',[1024 1536]);
                end
                test=subplot(numconds,5,h*5);
                plot((1536000:binwidth:1536000+minBins(5)*(binwidth-1))./1000,histo5(h,1:minBins(5)),'k');hold on%convert x-axis values to ms
                set(test,'YLim',[0 100+1.1*max(numTrialsCond)]);set(test,'XLim',[1536 1936]);set(test,'YTickLabel','');
                if h~=numconds(end)
                    set(test,'XTickLabel','');
                else
                    set(test,'XTick',[1536 1936]);
                    set(test,'XTickLabel',[1536 1936]);
                end
            end
            lengthBinsOver53;
            mean_act;
            
            spikeDataName=[num2str(channel),'_',num2str(floor(session)),'_',num2str(sampleContrast)]
            spikeDataFolder=fullfile('F:','PL','spikeData',animal,spikeDataName);
            saveText=['save ',spikeDataFolder,'.mat matarray'];
            eval(saveText);
            
            formats=[{'eps'} {'png'}];
            formats={'fig'};
            for j=1:1%2
                imageFolderName=fullfile('F:','PL','PSTHs',animal,num2str(channel),formats{j});
                if ~exist(imageFolderName,'dir')
                    mkdir(imageFolderName);
                end
                imageName=fullfile(imageFolderName,spikeDataName);
                %export_fig(imageName,formats{j});
                printtext=sprintf('print -d%s %s',formats{j},imageName);
                saveas(fig,imageName,'fig') 
            end
            
            %write ROC values and Weibull function constants to file:
            X=[0 0];
            ROCMatFileName=[num2str(channel),'_',num2str(sampleContrast),'_roc'];
            ROCMatFolder=fullfile('F:','PL','ROC_mat_files',animal);
            ROCMatPath=fullfile(ROCMatFolder,ROCMatFileName);
            matExists=0;
            rocValsMat=[];
            listing=dir(ROCMatFolder);
            if size(listing,1)>2%check whether the mat file exists for that channel and session
                for r=3:size(listing,1)
                    fileName=listing(r,1).name;
                    ind=find(fileName=='.');
                    if ind
                        chFolder=fileName(1:ind(1)-1);
                        if strcmp(ROCMatFileName,chFolder)
                            matExists=1;
                        end
                    end
                end
            end
            if matExists==1
                loadText=['load ',ROCMatPath,' rocValsMat'];
                eval(loadText);
            end
            addVar={[num2str(channel),'_',num2str(sampleContrast)]};
            addVar=[addVar floor(session)];
            addVar=[addVar rocvals];
            addVar=[addVar X];
            row=[];
            for rowNum=1:size(rocValsMat,1)
                if find(floor(rocValsMat{rowNum,2})==floor(session))
                    row=rowNum;
                end
            end
            if isempty(row)
                rocValsMat=[rocValsMat;addVar];
            else
                rocValsMat(row,:)=addVar;
            end
            saveText=['save ',ROCMatPath,'.mat rocValsMat'];
            eval(saveText);
            
            %pause
            close all hidden
        end
    end
end
