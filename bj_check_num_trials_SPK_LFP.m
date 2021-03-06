function bj_check_num_trials_SPK_LFP
%Written by Xing 11/04/14. Check that the number of trials matches between
%LFP and spike data.
errorSessions=cell(1,4);
animals=[{'blanco'} {'jack'}];
areas=[{'v4_1'} {'v1_1'}];
sessionsSkip=[311 322 355];
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        channels=main_channels(animal,area);
        sessions=main_raw_sessions_final(animal,area,[],0);
        for j=1:length(sessions)
            if isempty(find(sessions(j)==sessionsSkip))
                %             [dummy,test_contrasts,sampleContrasts] = session_metadata_alex(sessions(j), animal);
                sampleContrast=30;
                for i=1:length(channels)
                    loadText=['load M:\Xing\pl_LFP\',animal,'\',num2str(sessions(j)),'\SGL_trial_LFP_-512_1536_ch',num2str(channels(i)),'.mat event_arr'];
                    eval(loadText)
                    loadText=['load F:\PL\pl_corr_art_trials\',animal,'\',num2str(sessions(j)),'_corrtrialartifact.mat removeTrialsTimestamps'];
                    eval(loadText);%removeTrialsTimestamps column 1: NLX_TRIAL_START; column 21: NLX_TRIAL_END
                    spikeDataName=[num2str(channels(i)),'_',num2str(floor(sessions(j))),'_',num2str(sampleContrast)];
                    spikeDataFolder=fullfile('F:','PL','spikeData',animal,spikeDataName);
                    loadText=['load ',spikeDataFolder,'.mat matarray'];
                    eval(loadText);
                    event_arrSPK=[];
                    for h=1:size(matarray,1)
                        event_arrSPK=[event_arrSPK length(matarray{h,1})];
                        if event_arr(h)<length(matarray{h,1})
                            errorSessions=[errorSessions {animal} {sessions(j)} {channels(i)} {h}];
                        end
                    end
                    diffTrialNums=event_arr-event_arrSPK;
                    if sum(diffTrialNums)~=length(removeTrialsTimestamps)
                        errorSessions=[errorSessions;{animal} {sessions(j)} {channels(i)} {h}];
                    end
                end
            end
        end
    end
end
errorSessions