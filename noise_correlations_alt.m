function noise_correlations_alt
%Written by Xing 06/02/14.
%Calculates noise correlations across channels, for each trial, with a
%different (and I think incorrect) method- calculates z-scores based on
%each channel, rather than across channels.

originalSessList=0;
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
        if originalSessList==1
            if strcmp(animal,'blanco')&&strcmp(area,'v4_1')
                sessionNums=[307 308 311 313 314 318 320 321 329 330 331:1:341];% blanco V4
            end
        end
        
        for sessionInd=1:length(sessionNums)
            if sessionNums(sessionInd)~=322
                zChAct=cell(1,length(channels));
                rvals=[];
                for chPair1=1:length(channels)-1
                    fileName=['ch',num2str(channels(chPair1)),'_',num2str(sessionNums(sessionInd)),'_example_sample_test_act'];
                    matName=fullfile(folderName,fileName);
                    load(matName);
                    ch1Act=epoch2;
                    for chPair2=chPair1+1:length(channels)
                        fileName=['ch',num2str(channels(chPair2)),'_',num2str(sessionNums(sessionInd)),'_example_sample_test_act'];
                        matName=fullfile(folderName,fileName);
                        load(matName);
                        ch2Act=epoch2;
                        zCh1Act=[];
                        zCh2Act=[];
                        for condInd=1:length(testContrasts)
                            zCh1Act=[zCh1Act zscore(ch1Act{condInd,1})];
                            zCh2Act=[zCh2Act zscore(ch2Act{condInd,1})];
                        end
                        tempR=corrcoef(zCh1Act,zCh2Act);
                        rvals=[rvals tempR(2)];%for each pairwise comparison between channels
                    end
                end
                allSessR(sessionInd,:)=rvals;%compile across sessions
            end
        end
        numCombinedSess=5;
%         numCombinedSess=floor(length(sessionNums)/2);
        early1=[];
        late1=[];
        for sessInd=1:numCombinedSess
            early1=[early1 allSessR(sessInd,:)];
            late1=[late1 allSessR(length(sessionNums)-sessInd+1,:)];
        end
        mean(early1)
        mean(late1)
        bins=[-0.95:0.01:0.95];
        subplot(2,2,animalInd+(areaInd-1)*2);
        hist(early1(:),bins);
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','r','facealpha',0.5)
        hold on
        hist(late1(:),bins);
        h = findobj(gca,'Type','patch');
        set(h,'facealpha',0.5)
        xlim([-0.3 0.5]);
        allSessR=[sessionNums' allSessR];
        saveText=['save F:\PL\noise_trial_corr\',animal,'_',area,'_noise_R.mat allSessR'];
        eval(saveText);
    end
end