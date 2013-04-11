function bj_calc_durations
%Written by Xing 11/04/13 
%Check if performance differs depending on the duration of the
%sampel-test-interval.
%Durations of the ISI ony vary for Blanco V4_1 sessions (from around 529 to
%1035 ms).
% 11	NLX_STIM_OFF
% 12	NLX_TEST_ON
% 13	NLX_TEST_OFF
onExternalHD=0;
if onExternalHD==1
    rootFolder='G:\PL_backup_060413';
else
    rootFolder='F:';
end
numBins=6;
animals=[{'blanco'} {'jack'}];
areas=[{'v4_1'} {'v4_2'} {'v1_1'} {'v1_2'}];
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        [sampleContrasts testContrasts]=area_metadata(area);
        sessions = main_raw_sessions_final(animal,area,[],0);
        for sampleInd=1:length(sampleContrasts)
            durations=[];
            allP123=[]
            allP13=[];
            for i=1:length(sessions)
                loadText=['load F:\PL\vals_perf\',animal,'\vals_',num2str(sessions(i)),'.mat'];
                eval(loadText)
                if strcmp(area,'v1_2')
                    vals=[vals{1};vals{2};vals{3}];
                end
                for cond=1:size(testContrasts,2)
                    indCond=find(vals(:,4)==cond);
                    valsCond=vals(indCond,:);
                    duration=[];
                    for trial=1:size(valsCond,1)
                        duration(trial,1)=valsCond(trial,12)-valsCond(trial,11);
                        duration(trial,2)=valsCond(trial,13)-valsCond(trial,12);
                        if valsCond(trial,16)~=-1
                            duration(trial,3)=1000;%correct response
                        elseif valsCond(trial,17)~=-1
                            duration(trial,3)=0;%incorrect response
                        elseif valsCond(trial,16)==-1&&valsCond(trial,17)==-1
                            duration(trial,3)=-1000;%no response
                        end
                    end
                    duration=duration/1000;
                    durations{i,cond}=duration;
                    if strcmp(area,'v4_1')&&strcmp(animal,'blanco')
                        maxDuration=max(duration(:,1));%maximum possible ISI duration
                        minDuration=min(duration(:,1));%maximum possible ISI duration
                        durationsBins=minDuration:(maxDuration-minDuration)/numBins:maxDuration;%divide into a certain number of bins
                        for durationsBinsInd=1:numBins
                            binInd=duration(:,1)>=durationsBins(durationsBinsInd);
                            binInd2=binInd+(duration(:,1)<durationsBins(durationsBinsInd)+(maxDuration-minDuration)/numBins);
                            correct(durationsBinsInd,cond)=sum(duration(binInd2==2,3)==1)/(sum(duration(binInd2==2,3)==1)+sum(duration(binInd2==2,3)==0));%exclude trials where no response given
                        end
                    end
                end
                if strcmp(area,'v4_1')&&strcmp(animal,'blanco')
                    allCorrect{i}=correct;
                    durList=[];
                    condList=[];
                    for rowCount=1:numBins
                        durList=[durList;zeros(1,size(testContrasts,2))+rowCount];
                        condList=[condList;1:1:size(testContrasts,2)];
                    end
                    perfReshape=reshape(correct,1,size(correct,1)*size(correct,2));
                    durReshape=reshape(durList,1,size(durList,1)*size(durList,2));
                    condReshape=reshape(condList,1,size(condList,1)*size(condList,2));
                    [pDurBinsANOVA,table,stats]=anovan(perfReshape,{durReshape,condReshape})
                    allPBins(i,:)=pDurBinsANOVA;
                    perfReshape=reshape(correct([1 numBins],:),1,size(correct([1 numBins],:),1)*size(correct([1 numBins],:),2));
                    durReshape=reshape(durList([1 numBins],:),1,size(durList([1 numBins],:),1)*size(durList([1 numBins],:),2));
                    condReshape=reshape(condList([1 numBins],:),1,size(condList([1 numBins],:),1)*size(condList([1 numBins],:),2));
                    [pDurFirstlastANOVA,table,stats]=anovan(perfReshape,{durReshape,condReshape})
                    %[cDurANOVA(animalInd+2*(areaInd-1),:),m,h]=multcompare(stats)
                    if pDurFirstlastANOVA(1)<0.05
                        multcompare(stats)
                        pause
                    end
                    allPFirstlast(i,:)=pDurFirstlastANOVA;
                    close all hidden
                end
            end
            sessCorrect=[];
            sessDurList=[];
            sessCondList=[];
            for sessCount=1:length(allCorrect)
                sessCorrect=[sessCorrect;allCorrect{sessCount}];
                sessDurList=[sessDurList;durList];
                sessCondList=[sessCondList;condList];
            end
            sessCorrectReshape=reshape(sessCorrect,1,size(sessCorrect,1)*size(sessCorrect,2));
            durReshape=reshape(sessDurList,1,size(sessDurList,1)*size(sessDurList,2));
            condReshape=reshape(sessCondList,1,size(sessCondList,1)*size(sessCondList,2));
            [pDurSessANOVA,table,stats]=anovan(sessCorrectReshape,{durReshape,condReshape})
            %[cDurANOVA(animalInd+2*(areaInd-1),:),m,h]=multcompare(stats)
            close all hidden
            saveMatName=['ISI_test_durations_',area,'_',num2str(sampleContrasts(sampleInd)),'_',num2str(numBins),'bins'];
            saveMatFolder=fullfile(rootFolder,'PL','vals_perf',animal);
            saveMatPath=fullfile(saveMatFolder,saveMatName);
            if strcmp(area,'v4_1')&&strcmp(animal,'blanco')
                saveText=['save ',saveMatPath,' durations allCorrect allPBins allPFirstlast'];
            else
                saveText=['save ',saveMatPath,' durations'];
            end
            eval(saveText);
        end
    end
end
