function CRF_light(roving,animals,areas)
cutoff=1;
excludeSuppressed=0;
normaliseCh=0;
normaliseSpontan=0;
excludeNonmonotonic=0;
useISI=0;
analysisType='CRF_new';
onExternalHD=0;
if onExternalHD==1
    rootFolder='G:\PL_backup_060413';
else
    rootFolder='F:';
end
slopeNeuroNew=[];PNENew=[];diffPNENew=[];minRateNew=[];maxRateNew=[];chSSENew=[];
plotDiffC50_30=1;
calculateTangent=1;
startEndTime='_1024_to_1536';
if nargin<2||isempty(animals)
    animals=[{'blanco'} {'jack'}];
end
if nargin<3||isempty(areas)
    if roving==0
        areas=[{'v4_1'} {'v1_1'}];
    elseif roving==1
        areas=[{'v1_2_1'} {'v1_2_2'} {'v1_2_3'}];
    end
end
figure
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        [sampleContrasts testContrasts]=area_metadata(area);
        for sampleInd=1:length(sampleContrasts)
            sampleContrast=sampleContrasts(sampleInd);
            testContrast=testContrasts(sampleInd,:);
            channels=main_channels(animal,area);
            sessions=main_raw_sessions_final(animal,area,[],0);
            meanEpoch4=zeros(length(sessions),length(testContrasts));%array to store activity across channels for each condition and session
            for sessionInd=1:length(sessions)
                for cond=1:length(testContrasts)
                    for chInd=1:length(channels)
                        loadText=['load F:\PL\spikeData\',animal,'\',num2str(channels(chInd)),'_',num2str(sessions(sessionInd)),'_',num2str(sampleContrast),'.mat'];
                        eval(loadText);
                        epoch3=[];
                        epoch4=[];
                        normAct=[];
                        if chInd==1
                            allEpoch4=zeros(1,length(matarray{cond,3}));%initialise array when number of trials known
                        end
                        for n=1:length(matarray{cond,3})
                            temp=matarray{cond,3}{n}(matarray{cond,3}{n}<1024);
                            temp=temp>1024-300;
                            epoch3(n)=sum(temp)/300*1000;
                            epoch4(n)=length(matarray{cond,4}{n})/512*1000;
                            if useISI==1
                                normAct(n)=epoch4(n)-epoch3(n);
                            elseif useISI==0
                                normAct(n)=epoch4(n);
                            end
                        end
                        allEpoch4=allEpoch4+normAct;%sum across channels, subtract spontan
                    end
                    meanEpoch4(sessionInd,cond)=mean(allEpoch4/(length(channels)));%find mean for condition across trials
                end
                [slopeNeuroNew,PNENew,diffPNENew,minRateNew,maxRateNew,chSSENew,xvals,yvals]=nr_fitting(meanEpoch4(sessionInd,:),sampleContrast,testContrast,sessionInd,slopeNeuroNew,chSSENew,PNENew,minRateNew,maxRateNew,diffPNENew,plotDiffC50_30,calculateTangent,startEndTime,animal,area);
                plot(xvals,yvals,'k.');hold on
            end
            subFolder='crf_meanchannels';
            if normaliseCh==1
                subFolder=[subFolder,'_normCh'];
            end
            if normaliseSpontan==1
                if useISI==0
                    subFolder=[subFolder,'_normSpontan'];
                elseif useISI==1
                    subFolder=[subFolder,'_normISI'];
                end
            end
            if excludeSuppressed==1
                subFolder=[subFolder,'_excludeSuppressed'];
            end
            if excludeNonmonotonic==1
                subFolder=[subFolder,'_excludeNonmonotonic'];
            end
%             if useISI==0
%                 imagename=['cumulative_CRFs_old_new_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10)];
%             elseif useISI==1
%                 imagename=['cumulative_CRFs_new_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10)];
%             end
            folderpathname=fullfile(rootFolder,'PL',analysisType,animal,subFolder);
            if ~exist(folderpathname,'dir')
                mkdir(folderpathname);
            end
%             pathname=fullfile(rootFolder,'PL',analysisType,animal,subFolder,imagename);
%             printtext=sprintf('print -dpng %s.png',pathname);
%             set(gcf,'PaperPositionMode','auto')
%             eval(printtext);
            matname=['cumulative_CRFs_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10)];
            pathname=fullfile(rootFolder,'PL',analysisType,animal,subFolder,matname);
            saveText=['save ',pathname,'.mat meanEpoch4 slopeNeuroNew PNENew diffPNENew minRateNew maxRateNew chSSENew'];
            eval(saveText);
        end
    end
end

