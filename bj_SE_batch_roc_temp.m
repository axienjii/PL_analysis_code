function bj_SE_batch_roc_temp(animal,area,channels,sessions,icopy,ncopies,istart,drawFigs,plotRedArtifacts,generateSpikes)
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
if plotRedArtifacts==1
    plotRedArtifactsText='plotRedArtifacts\';
else
    plotRedArtifactsText=[];
end
verbose=0;
rewriteSessions=1;
checkFileDateMod=0;
if nargin<3 || isempty(channels)
    channels = main_channels(animal,area);
end
if nargin<4 || isempty(sessions)
    sessions = main_raw_sessions_final(animal,area,1,0);
end
if nargin<5 || isempty(icopy)
    icopy = 1;
end
if nargin<6 || isempty(ncopies)
    ncopies = 1;
end
if nargin<7 || isempty(istart)
    istart = 1;
end

if generateSpikes==1
    % for iter = icopy:ncopies:(length(channels)*length(sessions))
    for iter = icopy:(length(channels)*length(sessions))
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
%         try
            % Check if file with ROC values already exists
            if strcmp(area,'v1_2')
                sampleText=40;
            else
                sampleText=30;
            end
            output_fname = spikeData_finder(animal, session, channel,sampleText,plotRedArtifactsText);
            % If it does, skip this one
            if ~exist(output_fname, 'file')
                bj_SE_V1_2_roc4_temp(animal, channel, session, max(0,verbose-1),plotRedArtifacts)% if not, create it
            elseif rewriteSessions==1
                if checkFileDateMod==0
                    bj_SE_V1_2_roc4_temp(animal, channel, session, max(0,verbose-1),plotRedArtifacts)% overwrite it regardless of last date of modification
                else
                    fileinfo=dir(output_fname);
                    modDate=fileinfo.date;
                    if ~strcmp(modDate(4:6),'Apr')||str2double(modDate(1:2))<18||str2double(modDate(1:2))==18&&str2double(modDate(13:14))<19&&str2double(modDate(16:17))<08
                        bj_SE_V1_2_roc4_temp(animal, channel, session, max(0,verbose-1),plotRedArtifacts)% or overwrite it
                    end
                end
            end
%         catch ME
%             %If an error occured, output the details about the error
%             animal
%             session
%             channel
%             
%             disp(ME);
%             for i=1:length(ME.stack);
%                 last_error_made = ME.stack(i);
%                 fprintf('In ==> %s at %d.\n', last_error_made.file, last_error_made.line);
%             end
%         end
    end
end
if drawFigs==1
    formats=[{'eps'} {'png'}];
    formatsCommand=[{'epsc'} {'png'}];
    for iter = icopy:(length(channels)*length(sessions))
        if iter<istart
            continue;
        end
        isess = ceil(iter/length(channels));
        ichan = mod(iter-1,length(channels))+1;
        
        session = sessions(isess);
        channel = channels(ichan);
        if strcmp(area,'v1_2')
            sampleContrasts=[20 30 40];
        else
            sampleContrasts=30;
        end
        for i=1:length(sampleContrasts)
            %if ~exist(printFileName,'file')
            %         info=imfinfo(printFileName);
            %         modDate=info.FileModDate;
            %         if ~strcmp(modDate(4:6),'Mar')||str2double(modDate(1:2))<17||str2double(modDate(1:2))==17&&str2double(modDate(13:14))<0&&str2double(modDate(16:17))<24
            figFileName=['F:\PL\PSTHs\',plotRedArtifactsText,animal,'\',num2str(channel),'\fig\',num2str(channel),'_',num2str(session),'_',num2str(sampleContrasts(i)),'.fig'];
            try
                uiopen(figFileName,1)
                for j=1:2
                    folderName=['F:\PL\PSTHs\',plotRedArtifactsText,animal,'\',num2str(channel),'\',formats{j}];
                    if ~exist(folderName,'dir')
                        mkdir(folderName)
                    end
                    printFileName=['F:\PL\PSTHs\',plotRedArtifactsText,animal,'\',num2str(channel),'\',formats{j},'\',num2str(channel),'_',num2str(session),'_',num2str(sampleContrasts(i)),'.',formats{j}];
                    print(sprintf('-d%s',formatsCommand{j}),'-r300',printFileName)
                end
            catch ME
                disp(ME)
                loadText=['load F:\PL\PSTHs\',plotRedArtifactsText,'missingFigSessions.mat missingSessions'];
                eval(loadText)
                missingSessions=[missingSessions;{animal} {channel} {session} {ME}];
                saveText=['save F:\PL\PSTHs\',plotRedArtifactsText,'missingFigSessions.mat missingSessions'];
                eval(saveText)
            end
            %         end
            close all
        end
    end
end

checkImages=0;
if checkImages==1
    %check that new images have been drawn:
    for iter = icopy:(length(channels)*length(sessions))
        if iter<istart
            continue;
        end
        isess = ceil(iter/length(channels));
        ichan = mod(iter-1,length(channels))+1;
        
        session = sessions(isess);
        channel = channels(ichan);
        if strcmp(area,'v1_2')
            sampleContrasts=[20 30 40];
        else
            sampleContrasts=30;
        end
        for i=1:length(sampleContrasts)
            printFileName=['F:\PL\PSTHs\',animal,'\',num2str(channel),'\png\',num2str(channel),'_',num2str(session),'_',num2str(sampleContrasts(i)),'.png'];
            %if ~exist(printFileName,'file')
            info=imfinfo(printFileName);
            modDate=info.FileModDate;
            if ~strcmp(modDate(4:6),'Mar')||str2double(modDate(1:2))<17||str2double(modDate(1:2))==17&&str2double(modDate(13:14))<0&&str2double(modDate(16:17))<24
                printFileName
            end
        end
    end
end