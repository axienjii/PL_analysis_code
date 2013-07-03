function bj_SE_batch_cp(animals,areas,sessionNums)
%Written by Xing on 27/06/13
%Batch file for V1, V4_1 and V1_2 sessions.
%Calculate choice probability.

% animal='blanco';
% area='v4_1';
if nargin<1||isempty(animals)
    animals=[{'blanco'} {'jack'}];
    % animals={'blanco'};
end
if nargin<2||isempty(areas)
    areas=[{'v4_1'} {'v4_2'} {'v1_1'} {'v1_2'}];
    areas=[{'v4_1'} {'v1_1'} {'v1_2_1'} {'v1_2_2'} {'v1_2_3'}];
end
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
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        channels = main_channels(animal,area);
        if isempty(sessionNums)||nargin<4
            sessionNums = main_raw_sessions_final(animal,area,[],0);
        end
        excludeFewTrials=1;
        if excludeFewTrials==1
            if find(sessionNums==398)
                sessionNums=sessionNums(~(sessionNums==398));
            end
            if find(sessionNums==451)
                sessionNums=sessionNums(~(sessionNums==451));
            end
        end
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
                if strncmp(area,'v1',2)
                    closeCond=find(abs(testContrast-sampleContrast)<=5);
                elseif strncmp(area,'v4',2)
                    closeCond=find(abs(testContrast-sampleContrast)<5);
                end
                closeTestContrast=testContrast(closeCond);
                if strcmp(CPorTrialsAcrossSess,'TAS')
                    if strcmp(area,'v4_1')||strcmp(area,'v1_1')
                        figTASconds=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
                        set(figTASconds, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                    elseif strcmp(area,'v1_2_1')||strcmp(area,'v1_2_2')||strcmp(area,'v1_2_3')
                        for animalInd=1:2
                            for areaInd=1:2
                                figTASconds(animalInd+(areaInd-1)*2)=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
                                set(figTASconds(animalInd+(areaInd-1)*2), 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                            end
                        end
                    end
                end
                for i=1:length(sessionNums)
                    matFolder=['F:\PL\spikeData\',animal];
                    chStr=[num2str(channel),'_',num2str(sessionNums(i)),'_',num2str(sampleContrast),'.mat'];
                    matPath=fullfile(matFolder,chStr);
                    mat1Exists=0;
                    if exist(matPath,'file')
                        mat1Exists=1;
                    end
                    matFolder2=['F:\PL\spikeData\',animal,'\plotRedArtifacts\correct_trials_only'];
                    matPath2=fullfile(matFolder2,chStr);
                    mat2Exists=0;
                    if exist(matPath2,'file')
                        mat2Exists=1;
                    end
                    if mat1Exists==1&&mat2Exists==1
                        for row=1:size(epochTimes,1)
                            if epochTimes{row,1}==channels(h)
                                %                         cellEpochTimes=epochTimes(row,3);
                            end
                        end
                        valsText=['load ',matPath,' matarray'];
                        eval(valsText);
                        allTrialsArr=matarray;%all trials, regardless of correctness of response
                        valsText=['load ',matPath2,' matarray'];
                        eval(valsText);
                        corrTrialsArr=matarray;%correct response trials only
                        bj_SE_cp(channels(h),sessionNums(i),cellEpochTimes,minusSpon,allTrialsArr,corrTrialsArr,animal,area,sampleContrast,testContrast,closeCond)
                        end
                    end
                end
            end
        end
    end
end


