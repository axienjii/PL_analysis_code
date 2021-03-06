function bj_SE_xcorr_batch(animal,area,channels,sessionNums)
%Written by Xing 09/09/10
%Modified from blanco_SE_read_epochs_batch_roc.
%Batch file for running blanco_SE_xcorr (itself modified from use_time_periods).
%To check that the identities of individual cells recorded at each
%electrode remained consistent throughout the experiment, correlation
%analysis used to provide a quantitative measure of recording stability. 

if nargin<3 || isempty(channels)
    channels = main_channels(animal,area);
end
if nargin<4 || isempty(sessionNums)
    sessionNums = main_raw_sessions_final(animal,area,[],0);
end
sigma=8;%in ms
for i=1:length(channels)
    for j=1:length(sessionNums)
        try
            bj_SE_xcorr(animal,area,channels(i),sessionNums(j),sigma);%note that epochTimes(ind,2) is currently unused
            %         all_onset2=[all_onset2 onset2];
        catch ME
            disp(ME)
            load F:\PL\xcorr\missingSessions.mat missingSessions
            missingSessions=[missingSessions;{animal} {channels(i)} {sessionNums(j)} {ME}];
            save F:\PL\xcorr\missingSessions.mat missingSessions
        end
    end
end
% all_onset2

