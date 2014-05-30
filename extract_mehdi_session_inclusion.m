function extract_mehdi_session_inclusion
%Written by Xing 16/05/14 to read Mehdi's SNR arrays (based on MUA_E, extracted MUA) and
%identify which sessions are to be included in the final analysis for each
%channel. Channels are included if at least 80% of sessions had SNRs of at
%least 1. If a given channel met this criterion, then only sessions for
%which the SNR was at least 1 were included, for that channel.
roving=0;
animals=[{'blanco'} {'jack'}];
if roving==0
    areaTexts=[{'V4'} {'V1'}];
    areas=[{'v4_1'} {'v4_2'} {'v1_1'}];
    areas=[{'v4_1'} {'v1_1'}];
elseif roving==1
    areaTexts={'20' '30' '40'};
end
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        loadText=['load F:\PL\mehdi_snr\',animal,'_SNR_',area,'_N.mat'];
        eval(loadText)
        goodSNR=SNR(:,:)>=1;
        if strcmp(animal,'blanco')%exclude channel 15 from all of Blanco's data as recording site was ambiguous
           channels=V_Channels~=15;
           goodSNR=goodSNR(:,channels);
        end
        if strcmp(area,'v4_1')%exclude control session with horizontal stimulus
            if strcmp(animal,'blanco')
                vertSess=REC_PEN~=342;
            elseif strcmp(animal,'jack')
                vertSess=REC_PEN~=50;
            end
           goodSNR=goodSNR(vertSess,:);
        end
        goodChannels=[];
        goodSessions=[];
        for i=1:size(goodSNR,2)%for each channel            
           if sum(goodSNR(:,i))/size(goodSNR,1)>=0.8%at least 80% of sessions have SNRs of >=1
               goodChannels=[goodChannels;V_Channels(i)];
               goodSessions=[goodSessions;{REC_PEN(goodSNR(:,i)==1)}];
           end
        end
        setdiff(V_Channels,goodChannels);
        saveText=['save F:\PL\SNR\',animal,'_',area,'_includeSessions.mat goodChannels goodSessions'];
        eval(saveText)
    end
end