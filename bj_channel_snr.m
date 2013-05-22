function bj_channel_snr(cutoff)
%Written by Xing 20/05/13
%Calculate signal-to-noise ratios based on mean peak test-evoked response
%over SD of spontaneous activity.
%Generates list of which sessions to include for which channels.
analysisType='SNR';
onExternalHD=0;
if onExternalHD==1
    rootFolder='G:\PL_backup_060413';
else
    rootFolder='F:';
end
animals=[{'blanco'} {'jack'}];
areas=[{'v4_1'} {'v4_2'} {'v1_1'} {'v1_2'}];
areas=[{'v4_1'} {'v1_1'} {'v1_2_1'} {'v1_2_2'} {'v1_2_3'}];
binwidth=10;%10 ms
bins=512*2:binwidth:512*3;
sigma=3;
calcSNRs=1;
if calcSNRs==1
    for animalInd=1:length(animals)
        animal=animals{animalInd};
        for areaInd=1:length(areas)
            area=areas{areaInd};
            [sampleContrasts testContrasts]=area_metadata(area);
            channels=main_channels(animal,area);
            sessionNums=main_raw_sessions_final(animal,area,[],0);
            for chInd=1:length(channels)
                for sampleContrastsInd=1:length(sampleContrasts)
                    sampleContrast=sampleContrasts(sampleContrastsInd);
                    testContrast=testContrasts(sampleContrastsInd,:);
                    includeSessions=[];
                    excludeSessions=[];
                    for i=1:length(sessionNums)
                        matFolder=['F:\PL\spikeData\',animal];
                        allSponSD=[];
                        allHistoTest=[];
                        snr=[];
                        for condInd=1:length(testContrast)
                            condSpon=[];%list of spontaneous activity levels across trials for each condition
                            histoTest=[];%matrix of test PSTH activity across trials for each condition
                            chStr=[num2str(channels(chInd)),'_',num2str(sessionNums(i)),'_',num2str(sampleContrast),'.mat'];
                            matPath=fullfile(matFolder,chStr);
                            matExists=0;
                            if exist(matPath,'file')
                                matExists=1;
                            end
                            if matExists==1
                                valsText=['load ',matPath,' matarray'];
                                eval(valsText);
                                if size(matarray{condInd,2},1)~=size(matarray{condInd,4},1)
                                    pause%if number of trials are not equal
                                end
                                for n=1:size(matarray{condInd,2})
                                    condSpon(n)=length(matarray{condInd,1}{n})/512*1000;%rate of spontaneous activity in spikes/s
                                    [N X]=hist(matarray{condInd,4}{n},bins);%test-evoked activity
                                    if size(N,2)==1
                                        N=N';
                                    end
                                    histoTest(n,1:length(N))=N;%tally spikes across trials for each bin
                                end
                                allSponSD(condInd)=std(condSpon);
                                allHistoTest(condInd,:)=mean(histoTest,1);
                                allHistoTest(condInd,:)=gaussfit(sigma,0,allHistoTest(condInd,:)*1000/binwidth);%average activity per ms. note:*10 is an arbitrary scaling factor
                                snr(condInd)=max(allHistoTest(condInd,:))/allSponSD(condInd);
                            end
                        end
                        sessionSNR(i)=max(snr);
                        if sessionSNR(i)>=cutoff
                            includeSessions=[includeSessions;channels(chInd) sessionNums(i) sessionSNR(i)];
                        else
                            excludeSessions=[excludeSessions;channels(chInd) sessionNums(i) sessionSNR(i)];
                        end
                    end
                    if size(includeSessions,1)>=0.8*length(sessionNums)%if at least 80% of sessions have good SNR for that channel
                        matname=['good_SNR_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10),'.mat'];
                        pathname=fullfile(rootFolder,'PL',analysisType,animal,matname);
                        if exist(pathname,'file')
                            loadText=['load ',pathname,' includeSessionsAll'];
                            eval(loadText);
                            includeSessionsAll=[includeSessionsAll;includeSessions];
                        else
                            includeSessionsAll=includeSessions;
                        end
                        saveText=['save ',pathname,' includeSessionsAll'];
                        eval(saveText);
                        if ~isempty(excludeSessions)
                            %save info for excluded session as well:
                            matname=['poor_SNR_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10),'.mat'];
                            pathname=fullfile(rootFolder,'PL',analysisType,animal,matname);
                            if exist(pathname,'file')
                                loadText=['load ',pathname,' excludeSessionsAll'];
                                eval(loadText);
                                excludeSessionsAll=[excludeSessionsAll;excludeSessions];
                            else
                                excludeSessionsAll=excludeSessions;
                            end
                            saveText=['save ',pathname,' excludeSessionsAll'];
                            eval(saveText);
                        end
                    else
                        %save info for excluded channel as well:
                        matname=['poor_SNR_chs',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10),'.mat'];
                        pathname=fullfile(rootFolder,'PL',analysisType,animal,matname);
                        if exist(pathname,'file')
                            loadText=['load ',pathname,' excludeSessionsAll'];
                            eval(loadText);
                            excludeSessionsAll=[excludeSessionsAll;excludeSessions;includeSessions];
                        else
                            excludeSessionsAll=[excludeSessions;includeSessions];
                        end
                        saveText=['save ',pathname,' excludeSessionsAll'];
                        eval(saveText);
                    end
                end
            end
        end
    end
end

for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        [sampleContrasts testContrasts]=area_metadata(area);
        channels=main_channels(animal,area);
        for sampleContrastsInd=1:length(sampleContrasts)
            sampleContrast=sampleContrasts(sampleContrastsInd);
            meanSNR=[];
            matname=['good_SNR_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10),'.mat'];
            pathname=fullfile(rootFolder,'PL',analysisType,animal,matname);
            loadText=['load ',pathname,' includeSessionsAll'];
            eval(loadText);
            for chInd=1:length(channels)
                includeRows=includeSessionsAll(find(includeSessionsAll(:,1)==channels(chInd)),3);%include this session in analysis
                meanSNR(chInd,1:3)=[channels(chInd) mean(includeRows) std(includeRows)];
            end
            matname=['mean_SNR_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10),'.mat'];
            pathname=fullfile(rootFolder,'PL',analysisType,animal,matname);
            saveText=['save ',pathname,' meanSNR'];
            eval(saveText);            
        end
    end
end

                    