function bj_SE_xcorr_batch(animal,area)
%Written by Xing 09/09/10
%Modified from blanco_SE_read_epochs_batch_roc.
%Batch file for running blanco_SE_xcorr (itself modified from use_time_periods).
%To check that the identities of individual cells recorded at each
%electrode remained consistent throughout the experiment, correlation
%analysis used to provide a quantitative measure of recording stability. 

sessionNums=main_raw_sessions_final(animal,area,[],0);
channels=main_channels(animal,area);
sigma=8;%in ms
[sampleContrasts allTestContrasts]=area_metadata(area);
for i=1:length(channels)
    for j=1:length(sessionNums)
        for sampleContrastInd=1:length(sampleContrasts)
            sampleContrast=sampleContrasts(sampleContrastInd);
            testContrasts=allTestContrasts(sampleContrastInd,:);
            try
            bj_SE_xcorr(animal,area,channels(i),sessionNums(j),sigma,sampleContrast,testContrasts);%note that epochTimes(ind,2) is currently unused
            %         all_onset2=[all_onset2 onset2];
            catch ME
                disp(ME)
                load F:\PL\xcorr\missingFigSessions.mat missingSessions
                missingSessions=[missingSessions;{animal} {channels(i)} {sessionNums(j)} {ME}];
                save F:\PL\xcorr\missingSessions.mat missingSessions                
            end
        end
    end
end
% all_onset2

