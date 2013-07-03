function bj_SE_batch_write_trial_resp(animals,areas,channels,sessions,icopy,ncopies,istart,drawFigs,plotRedArtifacts)
%Written by Xing 01/07/2013
%Generate array containing values of 1 (correct) and -1 (incorrect) for
%each trial from which spike data is subsequently read and analysed. Check
%number of trials to ensure that procedure is carried out correctly.

plotRedArtifacts=0;
verbose=0;
writeResp=0;
checkLengths=1;
if nargin<1||isempty(animals)
    animals=[{'blanco'} {'jack'}];
end
if nargin<2||isempty(areas)
    areas=[{'v4_1'} {'v4_2'} {'v1_1'} {'v1_2'}];
    areas=[{'v4_1'} {'v1_1'} {'v1_2_1'} {'v1_2_2'} {'v1_2_3'}];
end
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        channels = main_channels(animal,area);
        channels=channels(1);%just run once, as the correctness of the response made on each trial is the same across all channels
        sessions = main_raw_sessions_final(animal,area,[],0);
        if strcmp(animal,'blanco')&&strcmp(area,'v1_1')
            ind=find(sessions==355);
            sessions(ind)=355.1;
            sessions=[sessions(1:ind) 355.2 sessions(ind+1:end)];
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
        
        if writeResp==1
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
%                 try
                    bj_write_trial_resp(animal, channel, session, max(0,verbose-1),plotRedArtifacts,area)% if not, create it
%                 catch ME
                    %If an error occured, output the details about the error
%                     animal
%                     session
%                     channel
%                     
%                     disp(ME);
%                     for i=1:length(ME.stack);
%                         last_error_made = ME.stack(i);
%                         fprintf('In ==> %s at %d.\n', last_error_made.file, last_error_made.line);
%                     end
%                 end
            end
        end
    end
end

if nargin<2||isempty(areas)
    areas=[{'v4_1'} {'v4_2'} {'v1_1'} {'v1_2'}];
    areas=[{'v4_1'} {'v1_1'} {'v1_2_1'} {'v1_2_2'} {'v1_2_3'}];
end
if checkLengths==1%compare number of trials between activity and response arrays, truncate arrays for sessions with large numbers of empty trials
    notEqual=[];
    adjustSessions=[311 318 333 352];
    adjustSessions=[333 352 398 451];
    for animalInd=1:length(animals)
        animal=animals{animalInd};
        for areaInd=1:length(areas)
            area=areas{areaInd};
            channels = main_channels(animal,area);
            sessionNums = main_raw_sessions_final(animal,area,[],0);
            [sampleContrasts testContrasts]=area_metadata(area);
            for i=1:length(sessionNums)
                if sessionNums(i)~=355&&sessionNums(i)~=405&&sessionNums(i)~=435
                    for sampleContrastsInd=1:length(sampleContrasts)
                        sampleContrast=sampleContrasts(sampleContrastsInd);
                        testContrast=testContrasts(sampleContrastsInd,:);
                        respName=[area,'_',num2str(sessionNums(i)),'_',num2str(sampleContrast)];
                        respDataFolder=fullfile('F:','PL','TAS',animal);
                        respDataFolder=fullfile(respDataFolder,respName);
                        loadBehavText=['load ',respDataFolder,'.mat respArray'];
                        eval(loadBehavText);%read list of incorrect and correct trials
                        for h=1:length(channels)
                            if strncmp(area,'v1_2',4)
                                loadText=['load F:\PL\sample_test_activity\',animal,'_v1_2\ch',num2str(channels(h)),'_',num2str(sessionNums(i)),'_example_sample_test_act.mat epoch2 epoch4'];
                            else
                                loadText=['load F:\PL\sample_test_activity\',animal,'_',area,'\ch',num2str(channels(h)),'_',num2str(sessionNums(i)),'_example_sample_test_act.mat epoch2 epoch4'];
                            end
                            eval(loadText);%load channel activity
                            for condInd=1:length(testContrast)
                                if size(epoch2{condInd},2)~=size(epoch4{condInd},2)||size(epoch2{condInd},2)~=size(respArray{condInd},1)
                                    if find(sessionNums(i)==adjustSessions)
                                        if size(epoch2{condInd},2)==size(epoch4{condInd},2)&&size(epoch2{condInd},2)<size(respArray{condInd},1)%if the numbers of trials in the matarrays are fine, truncate just the array containing responses
                                            respArray{condInd}=respArray{condInd}(1:size(epoch4{condInd},2));
                                            saveBehavText=['save ',respDataFolder,'.mat respArray'];
                                            eval(saveBehavText);%save new list of incorrect and correct trials
                                        else
                                            notEqual=[notEqual;animalInd areaInd sampleContrast channels(h) sessionNums(i) size(epoch2{condInd},2) size(epoch4{condInd},2) size(respArray{condInd},1)];
                                        end
                                    else
                                        notEqual=[notEqual;animalInd areaInd sampleContrast channels(h) sessionNums(i) size(epoch2{condInd},2) size(epoch4{condInd},2) size(respArray{condInd},1)];
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
% if drawFigs==1
%     formats=[{'eps'} {'png'}];
%     formatsCommand=[{'epsc'} {'png'}];
%     for iter = icopy:(length(channels)*length(sessions))
%         if iter<istart
%             continue;
%         end
%         isess = ceil(iter/length(channels));
%         ichan = mod(iter-1,length(channels))+1;
%         
%         session = sessions(isess);
%         channel = channels(ichan);
%         if strcmp(area,'v1_2')
%             sampleContrasts=[20 30 40];
%         else
%             sampleContrasts=30;
%         end
%         for i=1:length(sampleContrasts)
%             %if ~exist(printFileName,'file')
%             %         info=imfinfo(printFileName);
%             %         modDate=info.FileModDate;
%             %         if ~strcmp(modDate(4:6),'Mar')||str2double(modDate(1:2))<17||str2double(modDate(1:2))==17&&str2double(modDate(13:14))<0&&str2double(modDate(16:17))<24
%             figFileName=['F:\PL\PSTHs\',plotRedArtifactsText,animal,'\',num2str(channel),'\fig\',num2str(channel),'_',num2str(session),'_',num2str(sampleContrasts(i)),'.fig'];
%             try
%                 uiopen(figFileName,1)
%                 for j=1:2
%                     folderName=['F:\PL\PSTHs\',plotRedArtifactsText,animal,'\',num2str(channel),'\',formats{j}];
%                     if ~exist(folderName,'dir')
%                         mkdir(folderName)
%                     end
%                     printFileName=['F:\PL\PSTHs\',plotRedArtifactsText,animal,'\',num2str(channel),'\',formats{j},'\',num2str(channel),'_',num2str(session),'_',num2str(sampleContrasts(i)),'.',formats{j}];
%                     print(sprintf('-d%s',formatsCommand{j}),'-r300',printFileName)
%                 end
%             catch ME
%                 disp(ME)
%                 loadText=['load F:\PL\PSTHs\',plotRedArtifactsText,'missingFigSessions.mat missingSessions'];
%                 eval(loadText)
%                 missingSessions=[missingSessions;{animal} {channel} {session} {ME}];
%                 saveText=['save F:\PL\PSTHs\',plotRedArtifactsText,'missingFigSessions.mat missingSessions'];
%                 eval(saveText)
%             end
%             %         end
%             close all
%         end
%     end
% end
% 
% checkImages=0;
% if checkImages==1
%     %check that new images have been drawn:
%     for iter = icopy:(length(channels)*length(sessions))
%         if iter<istart
%             continue;
%         end
%         isess = ceil(iter/length(channels));
%         ichan = mod(iter-1,length(channels))+1;
%         
%         session = sessions(isess);
%         channel = channels(ichan);
%         if strcmp(area,'v1_2')
%             sampleContrasts=[20 30 40];
%         else
%             sampleContrasts=30;
%         end
%         for i=1:length(sampleContrasts)
%             printFileName=['F:\PL\PSTHs\',animal,'\',num2str(channel),'\png\',num2str(channel),'_',num2str(session),'_',num2str(sampleContrasts(i)),'.png'];
%             %if ~exist(printFileName,'file')
%             info=imfinfo(printFileName);
%             modDate=info.FileModDate;
%             if ~strcmp(modDate(4:6),'Mar')||str2double(modDate(1:2))<17||str2double(modDate(1:2))==17&&str2double(modDate(13:14))<0&&str2double(modDate(16:17))<24
%                 printFileName
%             end
%         end
%     end
% end