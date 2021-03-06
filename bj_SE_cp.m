function bj_SE_cp(ch,session,test_epochs,minusSpon,allTrialsArr,corrTrialsArr,animal,area,sampleContrast,testContrast,closeCond)
%Writes CP values to mat file:
%session number, time epochs, CP in spikes/s for each condition during
%spontaneous period, during sample, during ISI, and during test.
%
%Set the parameter 'minusSpon' to 1 to subtract spontaneous firing rates:
%spontaneous levels are calculated in units of spikes/s during -150 to 0 ms
%before sample onset. Average firing rates (also in spikes/s) are
%calculated for epochs 4 and 5, and spontaneous rates are subtracted.

durSpon=150;%length of period prior to sample onset from which spontaneous rates are calculated. Can take on a value of up to 512 ms.
minTrials=10;%set value of minumum number of trials for inclusion of session

%             if length(test_epochs)>2
%                 samp_epochs=[0 test_epochs(2)-1024 512];
%             else
%                 samp_epochs=[0 512];
%             end
numTrials=zeros(1,length(testContrast));
for cond=closeCond
    numTrials(cond)=min([length(corrTrialsArr{cond,1}) length(corrTrialsArr{cond,2}) length(corrTrialsArr{cond,3}) length(corrTrialsArr{cond,4}) length(corrTrialsArr{cond,5})]);
end
if min(numTrials(closeCond))>=minTrials
    for epoch=4:4%only examine test presentation period
        if epoch==1
            periods=[-durSpon 0];
        else
            periods=[test_epochs{epoch-1} test_epochs{epoch}(1)];
        end
        for subPeriod=1:length(periods)-1
            startTime=periods(subPeriod);
            endTime=periods(subPeriod+1);
            corrTrialsAct=cell(length(testContrast),1);
            incorrTrialsAct=cell(length(testContrast),1);
            incorrTrialsArr=cell(length(testContrast),5);
            CPvals=[];
            for cond=closeCond
                allTrialRowInd=1;
                for n=1:numTrials(cond)
                    corrTrialAct=corrTrialsArr{cond,4}{n};%test act
                    correctTrialExists=0;
                    %check that trial exists in new array
                    for allrowInd=1:size(allTrialsArr{cond,4},1)
                        if length(allTrialsArr{cond,4}{allrowInd})==length(corrTrialAct)
                            if allTrialsArr{cond,4}{allrowInd}==corrTrialAct
                                correctTrialExists=1;
                            end
                        end
                    end
                    %if that correct trial exists, get the spike act details, otherwise, skip it
                    if correctTrialExists==1
                        foundMatch=0;
                        while foundMatch==0&&allTrialRowInd<=size(allTrialsArr{cond,4},1)
                            if length(allTrialsArr{cond,4}{allTrialRowInd})==length(corrTrialAct)
                                if allTrialsArr{cond,4}{allTrialRowInd}==corrTrialAct
                                    corrTrialsAct{cond,1}=[corrTrialsAct{cond,1};length(corrTrialAct)];%store number of spikes during test epoch for correct trials
                                    allTrialRowInd=allTrialRowInd+1;
                                    foundMatch=1;
                                else
                                    incorrTrialsArr{cond,4}=[incorrTrialsArr{cond,4};allTrialsArr{cond,4}(allTrialRowInd)];%store spike times in new matarray for incorrect trials
                                    incorrTrialsAct{cond,1}=[incorrTrialsAct{cond,1};length(allTrialsArr{cond,4}{allTrialRowInd})];%store number of spikes during test epoch for incorrect trials
                                    allTrialRowInd=allTrialRowInd+1;
                                end
                            else
                                incorrTrialsArr{cond,4}=[incorrTrialsArr{cond,4};allTrialsArr{cond,4}(allTrialRowInd)];%store spike times
                                incorrTrialsAct{cond,1}=[incorrTrialsAct{cond,1};length(allTrialsArr{cond,4}{allTrialRowInd})];%store number of spikes during test epoch for incorrect trials
                                allTrialRowInd=allTrialRowInd+1;
                            end
                        end
                    end
                end
%                 if testContrast(cond)<sampleContrast
%                     [roc]=sglroc3(incorrTrialsAct{cond}',corrTrialsAct{cond}');
%                 elseif testContrast(cond)>sampleContrast
                    if isempty(incorrTrialsAct{cond})||isempty(corrTrialsAct{cond})
                        roc=NaN;
                    else
                        [roc]=sglroc3(corrTrialsAct{cond}',incorrTrialsAct{cond}');%if test contrast is below sample contrast, expect CP value to be less than 0.5, vice versa for higher test contrast
                    end
%                 end
                CPvals=[CPvals roc];%20 vs 30; 30 vs 40; 20 vs 40
            end
            startEndTime=['_',num2str(periods(subPeriod)),'_to_',num2str(periods(subPeriod+1))];
            if round(ch)~=ch
                CPmatName=['CP_Ch',num2str(round(ch)),'_',num2str(10*(ch-round(ch))),'_',num2str(sampleContrast),startEndTime,'.mat'];
            else
                CPmatName=['CP_Ch',num2str(ch),'_',num2str(sampleContrast),startEndTime,'.mat'];
            end
            CPmatFolder=fullfile('F:','PL','CP',animal,area);
            if ~exist(CPmatFolder,'dir')
                mkdir(CPmatFolder);
            end
            CPmatPath=fullfile('F:','PL','CP',animal,area,CPmatName);
            CPmatTemp=[{session} {test_epochs} {CPvals}];
            if ~exist(CPmatPath,'file')
                CPmat=CPmatTemp;
            elseif exist(CPmatPath,'file')
                loadText=['load ',CPmatPath,' CPmat'];
                eval(loadText);
                replace=[];
                for rowInd=1:size(CPmat,1)
                    if CPmat{rowInd,1}==session
                        replace=rowInd;
                    end
                end
                if ~isempty(replace)
                    CPmat(replace,:)=CPmatTemp;
                else
                    CPmat=[CPmat;CPmatTemp];
                end
            end
            saveText=['save ',CPmatPath,' CPmat'];
            eval(saveText);
        end
    end
end

