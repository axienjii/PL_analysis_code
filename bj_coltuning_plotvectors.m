function bj_coltuning_plotvectors
%Written by Xing 110511 for feature-selective attention analysis in MT,
%adapted for PL project with Blanco & Jack in V1 & V4.
%Reads data from .mat files and plots vectors of colour tuning on 2 plots:
%1. for all DS cells from teir and indi
%2. for all DS cells from ti that are found to have significant differences
%due to colour of stimulus, as calculated in an ANOVA in
%feat_sel_ti_col_1epoch.
%Note that ANOVA results in the function feat_sel_ti_col_1epoch are
%'contaminated' by direction of motion tuning preferences, as ANOVA did not
%take direction of stimulus motion into consideration. Thus, 2nd figure is
%just to give a rough visualisation of data that is not particularly
%meaningful.
animalTexts=[{'subject B'} {'subject J'}];
animals=[{'blanco'} {'jack'}];
areaTexts=[{'V4'} {'V1'}];
areas=[{'v4_1'} {'v1_2'}];
allV1sess=1;
figure
for animalInd=1:2
    animal=animals{animalInd};
    for areaInd=1:2
        area=areas{areaInd};
        folder=fullfile('F:','PL','STRF',animal);
        if allV1sess==1&&strcmp(area,'v1_2')
            allPrefOriMat='v1_1_and_v1_2_allChOri.mat';
        else
            allPrefOriMat=[area,'_allChOri.mat'];
        end
        allPrefOriPath=fullfile(folder,allPrefOriMat);
        loadText=['load ',allPrefOriPath,' allChOri'];
        eval(loadText)
        
        subplotInd=subplot(2,2,animalInd+2*(areaInd-1));
        ori=allChOri(:,2)*pi/180;
        rho=ones(size(allChOri,1),1);
        [x,y] = pol2cart(ori,rho);
        compass(x,y); hold on
        ylim([0 1]);
        if animalInd==1
            ylabel(areaTexts{areaInd});
        end
        if areaInd==1
            title(animalTexts{animalInd});
        else
            xlabel('orientation');
        end
    end
end

tuningName='STRF_polar_plots';
if allV1sess==1
    tuningName='STRF_polar_plots_combined_nonroving_&_roving';
end
tuningPath=fullfile('F:','PL','STRF',tuningName);
printtext=sprintf('print -dpng %s',tuningPath);
eval(printtext);

sumX=0;sumY=0;ind=[];colPref=[];
for i=1:size(colVectArrAll,1)
    sumX=sumX+colVectArrAll{i,4}/colVectArrAll{i,6};
    sumY=sumY+colVectArrAll{i,5}/colVectArrAll{i,6};
    if colVectArrAll{i,7}~=0
        ind=[ind i];
        colPref=[colPref colVectArrAll(i,7)];
    end
end
meanX=sumX/size(colVectArrAll,1);
meanY=sumY/size(colVectArrAll,1);

rayleighR=sqrt(meanX^2+meanY^2);
rayleighZ=size(colVectArrAll,1)*rayleighR^2;

% load 'H:\My Received Files\_work\cortex_files\alex_prg\cortex\feat_sel\xingtry\checks\teirplusindi\col_analysis\maineffect\col_preferences.mat' colVectArr
% figure
% for i=1:size(colVectArr,1)
%     xVect=colVectArr{i,6};
%     yVect=colVectArr{i,7};
%     plot([0 xVect],[0 yVect]);hold on
% end
% xlabel('G                                            R');ylabel('Y                                            B');
% axis square
% 
load 'H:\My Received Files\_work\cortex_files\alex_prg\cortex\feat_sel\xingtry\checks\teirplusindi\col_analysis\maineffect\col_preferences_rough.mat' colVectArrRough
figure
for i=1:size(colVectArrRough,1)
    xVect=colVectArrRough{i,6};
    yVect=colVectArrRough{i,7};
    plot([0 xVect],[0 yVect],'k');hold on
end
xlabel('G                                            R');ylabel('Y                                            B');
axis square

% load 'H:\My Received Files\_work\cortex_files\alex_prg\cortex\feat_sel\xingtry\checks\teirplusindi\col_analysis\maineffect\col_preferences_strict.mat' colVectArrSig

figure
for i=1:size(colVectArrAll,1)
    xVect=colVectArrAll{i,4};
    yVect=colVectArrAll{i,5};
    plot([0 xVect],[0 yVect],'k');hold on
end

for i=1:length(ind)
    xVect=colVectArrAll{ind(i),4};
    yVect=colVectArrAll{ind(i),5};
    plot([0 xVect],[0 yVect],'r');hold on
end
% for i=1:length(ind)
%     xVect=colVectArrAll{ind(i),4};
%     yVect=colVectArrAll{ind(i),5};
%     if colVectArrAll{ind(i),7}==1
%         lineCol='r';
%     elseif colVectArrAll{ind(i),7}==2
%         lineCol='b';
%     elseif colVectArrAll{ind(i),7}==3
%         lineCol='g';
%     elseif colVectArrAll{ind(i),7}==4
%         lineCol='y';
%     end
%     plot([0 xVect],[0 yVect],lineCol);hold on
% end
xlabel('G                                            R');ylabel('Y                                            B');
axis square




sumX=0;sumY=0;sumInterDir=0;%interaction between col and direction
for i=1:size(colVectArrRough,1)
    sumInterDir=sumInterDir+colVectArrRough{i,4};
    sumX=sumX+colVectArrRough{i,6}/colVectArrRough{i,8};
    sumY=sumY+colVectArrRough{i,7}/colVectArrRough{i,8};
end
meanX=sumX/size(colVectArrRough,1);
meanY=sumY/size(colVectArrRough,1);

rayleighR=sqrt(meanX^2+meanY^2);
rayleighZ=size(colVectArrRough,1)*rayleighR^2;
