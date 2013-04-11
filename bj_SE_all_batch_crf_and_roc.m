function bj_SE_all_batch_crf_and_roc(animal,area)
%Written by Xing on 27/09/10
%Modified from blanco_SE_all_batch_roc
%
%Batch file for V1, V4_1 and V4_2 sessions.
%Reads list of best epoch division times from representative sessions from file 'epoch_times_best3' in F:\blanco. 
%The .mat array of values for the session being analysed is loaded within this batch file.
%Calls function blanco_SE_crf and feeds in channel #, session#, and
%time points, as input args, finds response to stimuli of each contrast
%level for generation of contrast response function by 
%read_blanco_V4_1_crf_batch('F:\blanco\v4_1_crf_analysis'),
%read_blanco_V1_crf_batch('F:\blanco\v1_crf_analysis'), and
%read_blanco_V4_2_crf_batch('F:\blanco\v4_2_crf_analysis') later.
%Takes 'multipleTimes' as input arg 
%(0:average over whole epoch, 1:subdivide into time windows), 
%3rd input arg is minusSpon (0:don't subtract, 1:subtract).

% animal='blanco';
% area='v4_1';
multipleTimes=0;
minusSpon=0;
if multipleTimes==1
    bestSessionName='F:\blanco\epoch_times_best3';%contains list of times from best representative session for each channel
    count=1;
    fid=fopen(bestSessionName,'r');
    while ~feof(fid)
        [A]=fscanf(fid,'%s ', 1);%read channel
        if count==1
            epochTimes(count,1)={A};
        elseif count>1
            A=[num2str(values{1}(end)),A];
            epochTimes(count,1)={A};
        end
        values=textscan(fid,'%d','\n');
        a=values{1,1};
        if mod(length(values{1}),2)==1
            epochTimes(count,2)={values{1}(1:(length(values{1})-1)/2)};%epoch times from test onset to onset latency to test offset
            epochTimes(count,3)={values{1}((length(values{1})+1)/2:length(values{1})-1)};%epoch times from samp onset to onset latency to samp offset
        else
            epochTimes(count,2)={values{1}(1:(length(values{1}))/2)};%epoch times from test onset to onset latency to test offset
            epochTimes(count,3)={values{1}(length(values{1})/2+1:length(values{1}))};%epoch times from samp onset to onset latency to samp offset
        end
        count=count+1;
    end;
    fclose(fid);
    
    %get rid of session # and just keep ch # from epochTimes array:
    for i=1:size(epochTimes,1)
        index=find(epochTimes{i,1}=='_');
        if length(index)>1
            chNum=str2num(epochTimes{i,1}(1:index(1)-1))+0.1*str2num(epochTimes{i,1}(index(1)+1:index(2)-1))
        else
            chNum=str2num(epochTimes{i,1}(1:index(1)-1));
        end
        epochTimes(i,1)={chNum};
    end
end

% %if average activity levels are to be calculated across 512 ms in epoch 4
% %and 400 ms in epoch 5, i.e. without subdividing analysis periods into
% %smaller time windows:

channels = main_channels(animal,area);
sessionNums = main_raw_sessions_final(animal,area,[],0);
if multipleTimes==0
    for i=1:length(channels)
        epochTimes(i,1)={channels(i)};
    end
    epochTimes(:,3)={[0;512;512*2;512*3]};
%     for i=1:size(epochTimes,1)
%         epochTimes{i,3}(:)=[];
%         epochTimes(i,3)={[1024;1536;1936]};
%     end
end

cellEpochTimes={0 512 512*2 512*3};%{[0 40 300] 529 [529*2 529*2+40 529*2+300] 529*3}
[sampleContrasts testContrasts]=area_metadata(area);
for h=1:length(channels)
    channel=channels(h);
    for sampleContrastsInd=1:length(sampleContrasts)
        sampleContrast=sampleContrasts(sampleContrastsInd);
        testContrast=testContrasts(sampleContrastsInd,:);
        for i=1:length(sessionNums)
            matFolder=['F:\PL\spikeData\',animal];
            chStr=[num2str(channel),'_',num2str(sessionNums(i)),'_',num2str(sampleContrast),'.mat'];
            matPath=fullfile(matFolder,chStr);
            matExists=0;
            if exist(matPath,'file')
                matExists=1;
            end
            if matExists==1
                for row=1:size(epochTimes,1)
                    if epochTimes{row,1}==channels(h)
%                         cellEpochTimes=epochTimes(row,3);
                    end
                end
                valsText=['load ',matPath,' matarray'];
                eval(valsText);
                bj_SE_crf_and_roc(channels(h),sessionNums(i),cellEpochTimes,minusSpon,matarray,animal,area,sampleContrast,testContrast)
            end
        end
    end
end


