function bj_plot_adaptation_act(exampleFig,cutoff,animals,useISI,areas,excludeSuppressed,normalize)
%Written by Xing 26/07/13
%Reads activity across all channels, sessions and trials for a particular
%test contrast condition (the one with a test contrast immediately above
%the sample contrast), plots grand PSTH.
%Set useISI to 1: based on pre-test vs test, not on sample vs test.
%Set useISI to 0: sample vs test.
%Set excludeSuppressed to 1 to exclude channels with stimulus-evoked
%suppression, i.e. blanco 13, 24, 42 and jack 49.
smoothing=1;
plotSeparateSessions=1;%set to 0 to calculate average act across sessions; set to 1 to keep sessions distinct
onExternalHD=0;
if onExternalHD==1
    rootFolder='G:\PL_backup_060413';
else
    rootFolder='F:';
end
if nargin<3||isempty(animals)
    animals=[{'blanco'} {'jack'}];
    animalTexts=[{'Monkey 1'} {'Monkey 2'}];
    % animals={'blanco'};
end
if nargin<4||isempty(areas)
    areas=[{'v4_1'} {'v4_2'} {'v1_1'} {'v1_2'}];
    areas=[{'v4_1'} {'v1_1'} {'v1_2_1'} {'v1_2_2'} {'v1_2_3'}];
    areas=[{'v4_1'} {'v1_1'}];
end
sigma=3;
mu=0;
binWidth=10;
numbins=ceil(512/binWidth);
bins=0:binWidth:numbins*binWidth;
figROCconds=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
set(figROCconds, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        [sampleContrasts testContrasts]=area_metadata(area);
        channels=main_channels(animal,area);
        sessionNums=main_raw_sessions_final(animal,area,[],0);
        if strcmp(area,'v4_1')
            conds=9;%test contrast of 31%
        elseif strcmp(area,'v1_1')
            conds=9;%test contrast of 32%
        elseif strcmp(area,'v1_2_1')||strcmp(area,'v1_2_2')||strcmp(area,'v1_2_3')
            conds=[4;9;8];%test contrast of 22, 35, 42% for respective sample contrast conditions
        end
        subplot(2,2,animalInd+(areaInd-1)*2);
        for sampleContrastsInd=1:length(sampleContrasts)
            sampleContrast=sampleContrasts(sampleContrastsInd);
            testContrast=testContrasts(sampleContrastsInd,:); 
            %read in list of included channels
            if cutoff~=1
                matname=['good_SNR_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10),'.mat'];
                pathname=fullfile(rootFolder,'PL','SNR',animal,matname);
            else
                matname=['good_SNR_',area,'_',num2str(sampleContrast),'.mat'];
                pathname=fullfile(rootFolder,'PL','SNR',animal,'cutoff_SNR_1',matname);
            end
            loadText=['load ',pathname,' includeSessionsAll'];
            eval(loadText);
            colmapText=colormap(jet(size(testContrast,2)));
            colmapText=[colmapText(1,:);132/255 22/255 216/255;202/255 65/255 223/255;colmapText(3:7,:);157/255 212/255 61/255;colmapText(10:12,:);178/255 111/255 12/255;colmapText(end,:)];
            allSpikeHist1=zeros(1,numbins);
            allSpikeHist2=zeros(1,numbins);
            allSpikeHist3=zeros(1,numbins);
            allSpikeHist4=zeros(1,numbins);
            allSpikeHist5=zeros(1,numbins);
            count=0;%keep track of divisor
            for subplotMultiple=1:length(conds)
                condInd=conds(sampleContrastsInd,subplotMultiple);
                if plotSeparateSessions==1
                    sessHist=zeros(length(sessionNums),numbins*5);
                end
                for i=1:length(sessionNums)
                    sessCount=0;%keep track of number of trials per session
                    sessSpikeHist1=zeros(1,numbins);
                    sessSpikeHist2=zeros(1,numbins);
                    sessSpikeHist3=zeros(1,numbins);
                    sessSpikeHist4=zeros(1,numbins);
                    sessSpikeHist5=zeros(1,numbins);
                    matFolder=['F:\PL\spikeData\',animal];
                    if useISI==0
                        allEpoch2=[];%array to store summed up activity across channels
                    elseif useISI==1
                        allEpoch3=[];%array to store summed up activity across channels
                    end
                    allEpoch4=[];
                    for chInd=1:length(channels)
                        includeCh=1;
                        if excludeSuppressed==1
                            if strcmp(animal,'blanco')&&channels(chInd)==13||channels(chInd)==24||channels(chInd)==42
                                includeCh=0;
                            elseif strcmp(animal,'jack')&&channels(chInd)==49
                                includeCh=0;
                            end
                        end
                        chStr=[num2str(channels(chInd)),'_',num2str(sessionNums(i)),'_',num2str(sampleContrast),'.mat'];
                        matPath=fullfile(matFolder,chStr);
                        matExists=0;
                        if exist(matPath,'file')
                            matExists=1;
                        end
                        includeRows=includeSessionsAll(find(includeSessionsAll(:,1)==channels(chInd)),2);%include this session in analysis
                        includeRow=find(includeRows==sessionNums(i));
                        if matExists==1&&~isempty(includeRow)&&includeCh
                            valsText=['load ',matPath,' matarray'];
                            eval(valsText);
                            if useISI==0
                                if size(matarray{condInd,2},1)~=size(matarray{condInd,4},1)
                                    pause%if number of trials are not equal
                                end
                                if isempty(allEpoch2)
                                    allEpoch2=zeros(size(matarray{condInd,4}));
                                end
                            elseif useISI==1
                                if size(matarray{condInd,3},1)~=size(matarray{condInd,4},1)
                                    pause%if number of trials are not equal
                                end
                                if isempty(allEpoch3)
                                    allEpoch3=zeros(size(matarray{condInd,4}));
                                end
                            end
                            if normalize==1
                                maxAll=[];%find the highest firing rate across all conditions and trials, across both the sample and test presentation periods
                                for n=1:size(matarray{condInd,4})
                                    if useISI==0
                                        maxAll=[maxAll length(matarray{condInd,2}{n})*1000/512 length(matarray{condInd,4}{n})*1000/512];
                                    elseif useISI==1
                                        temp3=matarray{condInd,3}{n}>512*2-256;%activity during second half of ISI
                                        spikes=matarray{condInd,3}{n}(temp3);
                                        temp3=spikes<512*2;
                                        spikes=spikes(temp3);
                                        maxAll=[maxAll length(spikes)*1000/256 length(matarray{condInd,4}{n})*1000/512];
                                    end
                                end
                                maxval=max(maxAll)/100;
                            else
                                maxval=1;
                            end
                            if isempty(allEpoch4)
                                allEpoch4=zeros(size(matarray{condInd,4}));
                            end
                            for n=1:size(matarray{condInd,4})           
                                spikeHist=histc(matarray{condInd,1}{n},bins-512);
                                allSpikeHist1=allSpikeHist1+spikeHist(1:end-1);%sum up spontan activity within bins
                                spikeHist=histc(matarray{condInd,2}{n},bins);
                                allSpikeHist2=allSpikeHist2+spikeHist(1:end-1);%sum up sample activity within bins
                                spikeHist=histc(matarray{condInd,3}{n},bins+512);
                                allSpikeHist3=allSpikeHist3+spikeHist(1:end-1);%sum up ISI activity
                                spikeHist=histc(matarray{condInd,4}{n},bins+1024);
                                allSpikeHist4=allSpikeHist4+spikeHist(1:end-1);%sum up test activity
                                spikeHist=histc(matarray{condInd,5}{n},bins+1536);
                                allSpikeHist5=allSpikeHist5+spikeHist(1:end-1);%sum up post-test activity
                                count=count+1;    
                                spikeHist=histc(matarray{condInd,1}{n},bins-512);
                                sessSpikeHist1=sessSpikeHist1+spikeHist(1:end-1);%sum up spontan activity within bins
                                spikeHist=histc(matarray{condInd,2}{n},bins);
                                sessSpikeHist2=sessSpikeHist2+spikeHist(1:end-1);%sum up sample activity within bins
                                spikeHist=histc(matarray{condInd,3}{n},bins+512);
                                sessSpikeHist3=sessSpikeHist3+spikeHist(1:end-1);%sum up ISI activity
                                spikeHist=histc(matarray{condInd,4}{n},bins+1024);
                                sessSpikeHist4=sessSpikeHist4+spikeHist(1:end-1);%sum up test activity
                                spikeHist=histc(matarray{condInd,5}{n},bins+1536);
                                sessSpikeHist5=sessSpikeHist5+spikeHist(1:end-1);%sum up post-test activity
                                sessCount=sessCount+1;
                            end
                        end
                    end
                    if plotSeparateSessions==1
                        sessHist(i,:)=[sessSpikeHist1 sessSpikeHist2 sessSpikeHist3 sessSpikeHist4 sessSpikeHist5];
                        sessHist(i,:)=1000/binWidth*sessHist(i,:)/sessCount;%divide by total number of channels, trial, and sessions, convert to spikes/s
                    end
                end
                allTrialHist=[allSpikeHist1 allSpikeHist2 allSpikeHist3 allSpikeHist4 allSpikeHist5];%concatenate across whole trial; get rid of last 'bin' returned from histc as it contains remainder after last specified bin edge, i.e. zero
                allTrialHist=1000/binWidth*allTrialHist/count;%divide by total number of channels, trial, and sessions, convert to spikes/s
                if smoothing==1
                    fitted_vector=gaussfit(sigma,mu,allTrialHist);
                    for i=1:length(sessionNums)
                        sessHist(i,:)=gaussfit(sigma,mu,sessHist(i,:));
                    end
                else
                    fitted_vector=allTrialHist;
                end
                midBins=binWidth/2-512;
                if plotSeparateSessions==0
                    for binInd=1:length(fitted_vector)-2
                        midBins=[midBins midBins(binInd)+binWidth];
                    end
                    plot(midBins,fitted_vector(1:end-1),'Color','r','LineWidth',2);hold on%last bin contains number of spikes that equal to the final bin edge- can exclude it
                elseif plotSeparateSessions==1
                    for binInd=1:length(fitted_vector)-1
                        midBins=[midBins midBins(binInd)+binWidth];
                    end
                    for i=1:length(sessionNums)
                        plot(midBins,sessHist(i,:),'Color',[1-i/length(sessionNums) 0 i/length(sessionNums)],'LineWidth',0.75);hold on%last bin contains number of spikes that equal to the final bin edge- can exclude it
                    end
                end
                hold on
                xlim([-400 1936])
                yLimVals=get(gca,'YLim');
                ylim([0 yLimVals(2)]);
                if areaInd==1
                    title(animalTexts{animalInd});
                    if animalInd==1
                        xlabel('time (ms)')
                        ylabel('firing rate (spikes/s)');
                    end
                end
            end
        end
    end
end
