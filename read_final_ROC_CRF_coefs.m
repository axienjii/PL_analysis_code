function read_final_ROC_CRF_coefs(analysisType,roving,animals,areas,excludeSuppressed,useISI)
%Written by Xing 04/06/14.
%Analyse CRF and neurometric function parameters (slope, PNE/C50, etc)
%across sessions, including v4_0 sessions (304, 305, 22, 23). Session
%inclusion based on Mehdi's SNR>=1 criterion, includes sessions for which
%data were restored from Bluray disc.
rootFolder='F:';
sglroc3IndividualChs=1;
normalize=0;
normaliseCh=0;
normaliseSpontan=0;
excludeNonmonotonic=0;
cutoff=1;
if nargin<3||isempty(animals)
    animals=[{'blanco'} {'jack'}];
    % animals={'blanco'};
end
animalTexts=[{'Monkey 1'} {'Monkey 2'}];
animalCol=[{[0.4 0.4 0.4]} {'r'}];
figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.6]); %AUROC & DICAF vals against within-trial Rs
for animalInd=1:length(animals)
    animal=animals{animalInd};
    allSlopes=[];
    allMidpoints=[];%to cover either CRF or ROC analyses
    allMin=[];
    allMax=[];
    if animalInd>1||(nargin<4||isempty(areas))
        areas=[{'v4_1'} {'v4_2'} {'v1_1'} {'v1_2'}];
        if strcmp(animal,'jack')
            areas=[{'v4_0_1'} {'v4_0_2'} {'v4_1'}];%jack
        elseif strcmp(animal,'blanco')
            areas=[{'v4_0_1'} {'v4_0_2'} {'v4_0_3'} {'v4_1'}];%blanco
        end
        %         areas={'v1_1'};%v1
        if roving==1
            areas=[{'v1_2_1'} {'v1_2_2'} {'v1_2_3'}];
        end
    end
    for areaInd=1:length(areas)
        area=areas{areaInd};
        [sampleContrasts testContrasts]=area_metadata(area);
        if strcmp(analysisType,'CRF')
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
        elseif strcmp(analysisType,'ROC')
            if sglroc3IndividualChs==1
                subFolder='new_vs_old_sglroc3acrosschannels';
            elseif sglroc3IndividualChs==0
                if useISI==0
                    subFolder='new_vs_old_sglrocmeanchannels';
                elseif useISI==1
                    subFolder='new_ROC_useISI_meanchannels';
                end
            end
            if excludeSuppressed
                subFolder=[subFolder,'_excludeSuppressed'];
            end
            if normalize
                subFolder=[subFolder,'_normalised'];
            end
        end
        for sampleContrastsInd=1:length(sampleContrasts)
            sampleContrast=sampleContrasts(sampleContrastsInd);
            testContrast=testContrasts(sampleContrastsInd,:);
            if useISI==0
                matname=['cumulative_',analysisType,'s_old_new_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10)];
            elseif useISI==1
                matname=['cumulative_',analysisType,'s_new_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10)];
            end
            if strcmp(analysisType,'CRF')
                pathname=fullfile(rootFolder,'PL',analysisType,animal,subFolder,matname);
                loadText=['load ',pathname,'.mat allMeanEpoch4AcrossTrials slopeNeuroNew C50New diffC50New minRateNew maxRateNew chSSENew'];
            elseif strcmp(analysisType,'ROC')
                pathname=fullfile(rootFolder,'PL','ROC_zero_one',animal,subFolder,matname);
                if useISI==1
                    pathname=fullfile(rootFolder,'PL','ROC',animal,subFolder,matname);
                end
                loadText=['load ',pathname,'.mat all_rocvals_sglroc3 all_rocvals slopeNeuroNew PNENew diffPNENew minRateNew maxRateNew chSSENew threshold82higher threshold82lower'];
            end
            eval(loadText);
            allSlopes=[allSlopes slopeNeuroNew];
            if strcmp(analysisType,'CRF')
                allMidpoints=[allMidpoints C50New];
            elseif strcmp(analysisType,'ROC')
                allMidpoints=[allMidpoints PNENew];
            end
            allMin=[allMin minRateNew];
            allMax=[allMax maxRateNew];
        end
    end
    [h p]=corrcoef(allSlopes',1:length(allSlopes));
    slopeR(animalInd)=h(2);
    slopep(animalInd)=p(2);
    tableR(animalInd,1)=h(2);
    tablep(animalInd,1)=p(2);
    [h p]=corrcoef(allMidpoints',1:length(allMidpoints));
    midpointR(animalInd)=h(2);
    midpointp(animalInd)=p(2);
    tableR(animalInd,2)=h(2);
    tablep(animalInd,2)=p(2);
    [h p]=corrcoef(allMin',1:length(allMin));
    minR(animalInd)=h(2);
    minp(animalInd)=p(2);
    tableR(animalInd,3)=h(2);
    tablep(animalInd,3)=p(2);
    [h p]=corrcoef(allMax',1:length(allMax));
    maxR(animalInd)=h(2);
    maxp(animalInd)=p(2);
    tableR(animalInd,4)=h(2);
    tablep(animalInd,4)=p(2);
    subplot(2,4,1+4*(animalInd-1));
    plot(1:length(allSlopes),allSlopes,'ko','MarkerSize',4,'MarkerFaceColor',animalCol{animalInd});hold on
    subplot(2,4,2+4*(animalInd-1));
    plot(1:length(allMidpoints),allMidpoints,'ko','MarkerSize',4,'MarkerFaceColor',animalCol{animalInd});hold on
    subplot(2,4,3+4*(animalInd-1));
    plot(1:length(allMax),allMax,'ko','MarkerSize',4,'MarkerFaceColor',animalCol{animalInd});hold on
    subplot(2,4,4+4*(animalInd-1));
    plot(1:length(allMin),allMin,'ko','MarkerSize',4,'MarkerFaceColor',animalCol{animalInd});hold on
end
tableR
tablep
subplot(2,4,1);
title('Slope vs session');
ylabel('Monkey 1');
subplot(2,4,5);
ylabel('Monkey 2');
subplot(2,4,2);
title('Midpoint vs session');
subplot(2,4,3);
title('Maximum vs session');
subplot(2,4,4);
title('Minimum vs session');
imageName=['bj_across_channels_',analysisType,'_',area];
if useISI==1
    imageName=['bj_across_channels_',analysisType,'_',area,'_useISI'];
end
pathname=fullfile(rootFolder,'PL',analysisType,imageName);
printtext=sprintf('print -dpng %s.png',pathname);
set(gcf,'PaperPositionMode','auto')
eval(printtext);