function bj_check_numtrials_batch(animal,area)
%created on 30/01/13 for GitHub version control.
%modified from bj_SE_V1_2_batch_roc
%on 23/06/12 to analyse data from either jack or blanco.
%Written by Xing 13/01/11
%Analyses MUA (non-SUA) for sessions 380 and after. 380 to 384: 14 conds
%385 to 387: 12 conditions
%388 onwards: 36 conditions (roving paradigm)
%423 onwards: 6 conditions (roving)
%431 onwards: 36 conditions (roving with flankers)
%Remember to adjust folder to which ROC text files are written in function
%blanco_SE_V1_2_roc!!!
notChopped=0;
chopSessions=311;%311 & 333
comparisonType=1;
% animal='blanco';
% animal='jack';
% area='v1_2';
% animals=[{'blanco'} {'jack'}];
% animals=[{'jack'}];
% areas=[{'v4_1'} {'v4_2'} {'v1_1'} {'v1_2'}];
areas=[{'v4_1'}];
animals=[{'blanco'}];
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        channels=main_channels(animal,area);
        sessionNums = main_raw_sessions_final(animal,area,[],0);
        [sampleContrasts testContrasts]=area_metadata(area);
        for i=1:length(sessionNums)
            for sampleContrastsInd=1:length(sampleContrasts)
                sampleContrast=sampleContrasts(sampleContrastsInd);
                testContrast=testContrasts(sampleContrastsInd,:);
                for h=1:length(channels)
                    channel=channels(h);
                    if comparisonType==1
                        matFolder=['F:\PL\spikeData\',animal];
                    elseif comparisonType==2
                        matFolder=['F:\PL\spikeData\',animal,'\plotRedArtifacts\correct_trials_only'];
                    end
                    chStr=[num2str(channel),'_',num2str(sessionNums(i)),'_',num2str(sampleContrast),'.mat'];
                    matPath=fullfile(matFolder,chStr);
                    valsText=['load ',matPath,' matarray'];
                    if exist(matPath,'file')
                        eval(valsText);
                        for cond=1:length(testContrast)
                            for epoch=1:5
                                sizeAll(cond,epoch)=length(matarray{cond,epoch});
                            end
                        end
                        sizeAllChs{h}=sizeAll;%compile across channels
                    end
                end
                template=sizeAllChs{1};%arbitrarily use first channel as baseline for comparison
                
                for cond=1:length(testContrast)
                    allNumTrials=[];
                    for h=1:length(channels)
                        for epoch=1:5
                            if sessionNums(i)==chopSessions&&notChopped%find mimimum number of trials per condition across channels
                                allNumTrials=[allNumTrials sizeAllChs{h}(cond,epoch)];
                            else
                                numCheck=sizeAllChs{h}(cond,epoch);
                                if numCheck~=template(cond,epoch)
%                                     if channels(h)~=4&&sessionNums(i)~=322
                                    channels(h)
                                    cond
                                    epoch
                                    sessionNums(i)
%                                     end
                                end
                            end
                        end
                    end
                    if sessionNums(i)==chopSessions&&notChopped
                        minNumTrials(cond)=min(allNumTrials);
                    end
                end
            end%for sampleContrastsInd=1:length(sampleContrasts)
            if sessionNums(i)==chopSessions&&notChopped
                for h=1:length(channels)
                    channel=channels(h);
                    if comparisonType==1
                        matFolder=['F:\PL\spikeData\',animal];
                    elseif comparisonType==2
                        matFolder=['F:\PL\spikeData\',animal,'\plotRedArtifacts\correct_trials_only'];
                    end
                    chStr=[num2str(channel),'_',num2str(sessionNums(i)),'_',num2str(sampleContrast),'.mat'];
                    matPath=fullfile(matFolder,chStr);
                    valsText=['load ',matPath,' matarray'];
                    eval(valsText);
                    for cond=1:length(testContrast)
                        for epoch=1:5
                            matarray{cond,epoch}=matarray{cond,epoch}(1:minNumTrials(cond),1);%include only trials that are common across channels, based on minimum number available
                        end
                    end
                    valsText=['save ',matPath,' matarray'];
                    eval(valsText);
                end
            end
        end
    end
end
