function bj_SE_batch_roc_temp(animal,area,channels,sessions,icopy,ncopies,istart)
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

% animal='blanco';
% animal='jack';
% area='v1_2';
verbose=0;
rewriteSessions=1;
if nargin<3 || isempty(channels)
    channels = main_channels(animal,area);
end
if nargin<4 || isempty(sessions)
    sessions = main_raw_sessions(animal,area);
end
if nargin<5
    icopy = 1;
end
if nargin<6
    ncopies = 1;
end
if nargin<7
    istart = 1;
end

% for iter = icopy:ncopies:(length(channels)*length(sessions))    
parfor iter = icopy:(length(channels)*length(sessions)) 
    if iter<istart
        continue;
    end    
    isess = ceil(iter/length(channels));
    ichan = mod(iter-1,length(channels))+1;
    
    session = sessions(isess);
    channel = channels(ichan);
    
    if verbose; fprintf('%d(%d): %d/%d Processing session %.1f (%d/%d), channel %d (%d/%d)\n',...
            icopy, ncopies, iter, length(channels)*length(sessions), ...
            session, isess, length(sessions), ...
            channel, ichan, length(channels));
    end
%     
%     try
        % Check if file with ROC values already exists
        if strcmp(area,'v1_2')
            sampleText=40;
        else
            sampleText=30;
        end
        output_fname = spikeData_finder(animal, session, channel,sampleText);
        % If it does, skip this one
        if ~exist(output_fname, 'file')
            bj_SE_V1_2_roc4_temp(animal, channel, session, max(0,verbose-1))% if not, create it
        elseif rewriteSessions==1 
            bj_SE_V1_2_roc4_temp(animal, channel, session, max(0,verbose-1))% or overwrite it
        end        
%     catch ME
        % If an error occured, output the details about the error
%         animal
%         session
%         channel
%         
%         disp(ME);
%         for i=1:length(ME.stack);
%             last_error_made = ME.stack(i);
%             fprintf('In ==> %s at %d.\n', last_error_made.file, last_error_made.line);
%         end
%     end
end
