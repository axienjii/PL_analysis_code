function[rocvals]=jack_SE_V4_1_roc(folder,channel,file_of_int,vals,sampleContrast,testContrast,session,skip,area,analysisFolderAppend,BUPLanalysisType)
%Modified on 13/02/2012 for Jack data analysis,
%from blanco_SE_V4_1_roc.
%Called by function blanco_SE_V1_batch_roc, checks folder specified by first
%input arg, and if it contains a sorted spike file, analyses spike
%activity from that file, for specified channel.

%Written by Xing on 08/05/10
%Modified from blanco_SE_smallbins, compares spike activity between sample
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

listing=dir(folder);
if size(listing,1)>2
    for q=3:size(listing,1)
        fileName=listing(q,1).name;
        ch=fileName(9:10);
        if ch(end)=='_'
            ch=ch(1);
        end
        if ch==num2str(round(channel))
            name=[folder,'\',fileName];
            [all_SE_TimeStamps,ScNumbers,CellNumbers,Params,DataPoints,NlxHeader]=Nlx2MatSpike(name,[1 1 1 1 1],1,1,1);
%             for k=1:1
                for k=1:max(CellNumbers)
                fig=figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.2, 0.04, 0.5, 0.9]);
                set(fig,'PaperUnits','centimeters','PaperType','A4','PaperOrientation', 'portrait', 'PaperPosition', [0.63452 0.63452 21 28.41]);
                numconds=length(testContrast);
                binwidth=10*1000;%10 ms
                lengthBinsOver53=0;lengthBinsOver100=0;
                minBins=zeros(1,5)+56;%arbitrarily large number
                mean_act=zeros(numconds,5);
                temp1=find(CellNumbers==k);
                SE_TimeStamps=all_SE_TimeStamps(temp1);%time stamp arrays for each of the sorted cells
                all_spontan=[];all_sample=[];all_fix2=[];all_test=[];all_fix3=[];
                cond_arr=[];
                numTrialsCond=zeros(1,numconds);
                histo1=zeros(numconds,100);histo2=histo1;histo3=zeros(numconds,200);histo4=histo1;histo5=histo1;
                rocvals=[];%for roc analysis
                rocvalsOld=[];
                for h=1:numconds
                    spontan=[];sample_act=[];fix2_act=[];test_act=[];fix3_act=[];
                    for i=1:size(vals,1)
                        if isempty(find(skip==i,1))
                            if vals(i,4)==h%check condition number
                                %     for j=1:size(vals,2)
                                %     SE_EV_TimeStamps(j)=find(SE_TimeStamps<vals(i,j+1)&&SE_TimeStamps>=vals(i,j));
                                %     spikes(j)=DataPoints(:,:,SE_EV_TimeStamps(j));
                                %     end
                                if (vals(i,22)>-1)
                                    numTrialsCond(1,h)=numTrialsCond(1,h)+1;%tally number of trials per condition
%                                     if numTrialsCond(1,h)<11%temporary measure to get equal number of BU and PL trials
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
                                        plot([spikeTimes(m)/1000 spikeTimes(m)/1000],[100+numTrialsCond(h)-0.5 100+numTrialsCond(h)+0.5],'k');hold on
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
                                        plot([spikeTimes(m)/1000 spikeTimes(m)/1000],[100+numTrialsCond(h)-0.5 100+numTrialsCond(h)+0.5],'k');hold on
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
                                        plot([spikeTimes(m)/1000 spikeTimes(m)/1000],[100+numTrialsCond(h)-0.5 100+numTrialsCond(h)+0.5],'k');hold on
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
                                        plot([spikeTimes(m)/1000 spikeTimes(m)/1000],[100+numTrialsCond(h)-0.5 100+numTrialsCond(h)+0.5],'k');hold on
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
                                        plot([spikeTimes(m)/1000 spikeTimes(m)/1000],[100+numTrialsCond(h)-0.5 100+numTrialsCond(h)+0.5],'k');hold on
                                    end
                                    if numTrialsCond(h)==1
                                        matarray(h,5,:)={spikeTimes/1000};
                                    else
                                        matarray(h,5,:)={[matarray{h,5}(:,:);{spikeTimes/1000}]};
                                    end
                                    fix3_act=[fix3_act length(temp2)*1000000/(vals(i,14)-vals(i,13))];
%                                     end
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
                    [roc]=sglroc3(test_act(:)',sample_act(:)');%old method
                    rocvalsOld=[rocvalsOld roc];
                    roc=sum(test_act>sample_act)/(sum(test_act>sample_act)+sum(test_act<sample_act));%new trial-wise method
                    rocvals=[rocvals roc];
                end
                rocvals
                size(rocvals)
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
                    set(test,'YLim',[0 100+1.1*max(numTrialsCond)]);set(test,'XLim',[-512 0]);set(test,'YTick',[0 100],'YTickLabel',[0 round(maxhisto)]);
                    if h~=numconds(end)
                        set(test,'XTickLabel','');
                    else
                        set(test,'XTick',[-512 0]);
                        set(test,'XTickLabel',[-512 0]);
                    end
                    if h==1
                        ptext=sprintf('%s   %s',name,file_of_int);
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
                %     %set(test,'YTickLabel',[0 20]);
                lengthBinsOver53
                mean_act
                name=[folder,'_',fileName];
                if max(CellNumbers)>1
                    filename1=[name(1:end-4),'_',num2str(k),'_act',analysisFolderAppend];
                else
                    filename1=[name(1:end-4),'_act',analysisFolderAppend];
                end
                currdir=cd;
                cd 'F:\PL\BUPL\j_mat_files'
               if round(channel)~=channel
                    filename2=[num2str(ch),'_',num2str(10*(channel-round(channel)))];
                else
                    filename2=[num2str(channel)];
                end  
                savetext=sprintf('save %s_%s_%s.mat matarray',filename2,num2str(session),BUPLanalysisType)
                eval(savetext);
                chdirtext=sprintf('cd ''%s''',currdir);
                eval(chdirtext);
                cd
                %             if max(CellNumbers)>1
                %                 filename1=[ch,'_',num2str(k),'_act'];
                %             else
                %                 filename1=[ch,'_act'];
                %             end
                %             if ch<10
                %                 filename1=['0',filename1];
                %             end
                %             printtext=sprintf('print -dpdf %s',filename1);
                printtext=sprintf('print -dpng %s',filename1);
                eval(printtext);
                %     fig=figure('Color',[1,1,1]);
                %     for h=1:numconds
                %         plot(1:1:5,mean_act(h,:),'Color',[1-h/numconds,1-h/numconds,1-h/numconds]);hold on
                %     end
                %     max_mean_act=max(max(mean_act))
                %     x=1.4:0.005:1.7;
                %     line(x,max_mean_act/8,'Color',[1-1/numconds,1-1/numconds,1-1/numconds]);hold on
                %     line(x,max_mean_act/12,'Color',[0 0 0]);hold on
                %     ptext1=sprintf('lowest contrast: %d%%',testContrast(1));
                %     text('Position',[1.8 max_mean_act/8],'FontSize',9,'String',ptext1);
                %     ptext2=sprintf('highest contrast: %d%%',testContrast(end));
                %     text('Position',[1.8 max_mean_act/12],'FontSize',9,'String',ptext2);
                %     titletext=sprintf('activity during 5 periods')
                %     set(fig,'Name',titletext);
                %     xlabel('time period');ylabel('spikes/s');
                %     text('Position',[1 max_mean_act*1.2],'FontSize',9,'String',ptext);
                %     [p1(1) table stats]=anovan(all_spontan,{cond_arr});%,'model','interaction');
                %     [p1(2) table stats]=anovan(all_sample,{cond_arr});%,'model','interaction');
                %     [p1(3) table stats]=anovan(all_fix2,{cond_arr});%,'model','interaction');
                %     [p1(4) table stats]=anovan(all_test,{cond_arr});%,'model','interaction');
                %     [p1(5) table stats]=anovan(all_fix3,{cond_arr});%,'model','interaction');
                %     orient landscape
                %     text('Position',[5.05 max_mean_act/8],'FontSize',7,'String',p1);
                %     p1
                %     filename1=[filename1(1:end-4),'_bin'];
                %     printtext=sprintf('print -dpdf %s',filename1);
                %     %eval(printtext);
                
                
                %write ROC values and Weibull function constants to file:
                currdir=cd
                matFolderText=['F:\PL\BUPL\',area,'_roc_analysis',analysisFolderAppend]; 
                matName=['ROCvals',analysisFolderAppend,'_first10.mat'];
                matName=['ROCvals',analysisFolderAppend,'_all.mat'];
                matPath=fullfile(matFolderText,matName);
                if exist(matPath,'file')
                    loadText=['load ',matPath];
                    eval(loadText)
                    if ~isempty(ROCvals)
                        ROCvals=[ROCvals;channel session rocvals];
                    else
                        ROCvals=[channel,session,rocvals];
                    end
                else
                    ROCvals=[channel,session,rocvals];
                end
                saveText=['save ',matPath,' ROCvals'];
                eval(saveText)
                
                chdirtext=sprintf('cd ''%s''',currdir);
                eval(chdirtext);
                cd
                
                %pause
                close all hidden
            end
        end
    end
end

