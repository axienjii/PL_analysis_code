function noise_correlations
%Written by Xing 05/02/14.
%Calculates noise correlations across channels, for each trial

figure
animals=[{'blanco'} {'jack'}];
areas=[{'v4_1'} {'v1_1'}];
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        allSessR=[];
        folderName=['F:\PL\sample_test_activity\',animal,'_',area];
        channels = main_channels(animal,area);
        sessionNums = main_raw_sessions_final(animal,area,[],0);
        [sampleContrasts testContrasts]=area_metadata(area);
        
        for sessionInd=1:length(sessionNums)
            if sessionNums(sessionInd)~=322
                zChAct=cell(1,length(channels));
                for condInd=1:length(testContrasts)
                    allChSampleAct=[];
                    for chInd=1:length(channels)
                        fileName=['ch',num2str(channels(chInd)),'_',num2str(sessionNums(sessionInd)),'_example_sample_test_act'];
                        matName=fullfile(folderName,fileName);
                        load(matName);
                        allChSampleAct=[allChSampleAct epoch2{condInd,1}];
                    end
                    zAllChSampleAct=zscore(allChSampleAct);%z-score for each condition
                    for chInd=1:length(channels)%compile z-scored values across conditions
                        zChAct{1,chInd}=[zChAct{1,chInd} zAllChSampleAct(1,1+(chInd-1)*length(epoch2{condInd,1}):chInd*length(epoch2{condInd,1}))];
                    end
                end
                rvals=[];
                for chPair1=1:length(channels)-1
                    for chPair2=chPair1+1:length(channels)
                        tempR=corrcoef(zChAct{1,chPair1},zChAct{1,chPair2});
                        rvals=[rvals tempR(2)];%for each pairwise comparison between channels
                    end
                end
                allSessR(sessionInd,:)=rvals;%compile across sessions
            end
        end
        numCombinedSess=5;
        early=[];
        late=[];
        for sessInd=1:numCombinedSess
            early=[early allSessR(sessionInd,:)];
            late=[late allSessR(length(sessionNums)-sessionInd+1,:)];
        end
        mean(early)
        mean(late)
        subplot(2,2,animalInd+(areaInd-1)*2);
        hist(early,50);
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','r','facealpha',0.5)
        hold on
        hist(late,50);
        h = findobj(gca,'Type','patch');
        set(h,'facealpha',0.5)
    end
end