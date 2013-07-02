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
analysisType='ROC_zero_one';
calcPartial=0;
plotAbsChange=1;
animalTexts=[{'subject B'} {'subject J'}];
animals=[{'blanco'} {'jack'}];
areaTexts=[{'V4'} {'V1'}];
areas=[{'v4_1'} {'v1_2'}];
allV1sess=1;
fig1=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.4, 0.35]); %
set(fig1, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
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
        compass(x,y,'k'); hold on
        ylim([0 1]);
        if animalInd==1
            ylabel(areaTexts{areaInd});
        end
        if areaInd==1
            title(animalTexts{animalInd});
        else
            xlabel('orientation');
        end
        
        %mark channels with significant changes
        if strcmp(area,'v4_1')
            areaText=area;
        else
            areaText='v1_1';
        end
        tallyMatname=['tallySigCoefs_',areaText];
        if strcmp(analysisType,'CRF')
            %                     if calcPartial==0
            subFolderName=[analysisType,'_coef_mat'];
            %                     elseif calcPartial==1
            %                         subFolderName=[analysisType,'_coef_mat_partialcorr'];
            %                     end
            coefMatFolder=fullfile('F:','PL',analysisType,animal,subFolderName);
            coefMatPathname=fullfile(coefMatFolder,tallyMatname);
        elseif strcmp(analysisType,'ROC_zero_one')
            if calcPartial==0
                coefMatPathname=['F:\PL\ROC_zero_one\',animal,'\ROC_zero_one_coef_mat\',tallyMatname];
            elseif calcPartial==1
                coefMatPathname=['F:\PL\ROC_zero_one\partial_corr_4_factors\',animal,'\ROC_zero_one_coef_mat_partialcorr_4factors\',tallyMatname];
            end
        end
        loadText=['load ',coefMatPathname,'.mat sigChs sigDirection dir30 sigCoefs sigChNamesDir sigChDir30 formattedTable formattedChsTable'];
        eval(loadText)
        
        edges=linspace(0,180,4);%create a certain number of bins (specify number of edges)
        expectedCounts=zeros(1,length(edges)-1)+size(allChOri,1)/(length(edges)-1);
        [h,p,st]=chi2gof(allChOri(:,2),'edges',edges,'expected',expectedCounts)
        pBin60(animalInd+2*(areaInd-1))=p;
        statsBin60{animalInd+2*(areaInd-1)}=st;
        
        edges=[-45 45 135];%create 2 bins (within 45 degrees of vertical or within 45 degrees of orthogonal/horizontal)
        expectedCounts=zeros(1,length(edges)-1)+size(allChOri,1)/(length(edges)-1);
        convertInd=find(allChOri(:,2)>135);%identify channels with PO of more than 135
        convertedOri=allChOri(:,2);
        convertedOri(convertInd)=convertedOri(convertInd)-180;%and convert them to lie adjacent (CW) to the value of 0 degrees
        [h2,p2,st2]=chi2gof(convertedOri,'edges',edges,'expected',expectedCounts)
        pBin90(animalInd+2*(areaInd-1))=p2;
        statsBin90{animalInd+2*(areaInd-1)}=st2;
        
        convertedOri2=2*allChOri(:,2);%Rayleigh test, convert to range of 0 to 360 instead of 0 to 180, can identify unimodal deviation from uniformity
        convertedOri2=convertedOri2(~isnan(convertedOri2))/180*pi;
        pRayleigh(animalInd+2*(areaInd-1))=circ_rtest(convertedOri2);
        
        %Omnibus test- handles multimodal deviations from uniformity
        pOmnibus(animalInd+2*(areaInd-1))=circ_otest(convertedOri2);
        
        paramRows=[1 2 6 7;3 4 8 9];
        rowCols='rmbc';
        rowCols=[1 0 0;1 180/255 200/255;0 0 1;180/255 220/255 1];
        if plotAbsChange==1
            for absOrRel=1:2%plot sig channels based on absolute direction of change (value increase/decrease) or relative (steepness increase/decrease; PNE shift to or away from 30%)
                for paramInd=1:size(paramRows,2)%two directions of change for each of two params (slope and PNE)
                    %sig changes in slope
                    sigChInd=[];
                    for chInd=1:length(formattedChsTable{paramRows(absOrRel,paramInd),:})%slope increases in value
                        [sigChInd(chInd) dummy]=find(formattedChsTable{paramRows(absOrRel,paramInd)}(chInd)==allChOri(:,1));
                    end
                    if ~isempty(sigChInd)
                        sigChOri=ori(sigChInd);
                        rho=ones(size(sigChOri,1),1)-0.25*absOrRel;%slightly shorter stem
                        [x,y] = pol2cart(sigChOri,rho);
                        handle=compass(x,y); hold on
                        set(handle, {'Color'},num2cell(rowCols(paramInd,:),2), 'LineWidth',2)
                        
                        %Watson-Williams test to compare means of two
                        %samples. N.A. if not enough data 
                        convertedSigChOri=sigChOri*2;
                        convertedSigChOri=convertedSigChOri(~isnan(convertedSigChOri))/180*pi;
                        [p stats]=circ_wwtest(convertedOri2,convertedSigChOri);
                        pWW(animalInd+2*(areaInd-1),paramInd)=p;
                        statsWW{animalInd+2*(areaInd-1),paramInd}=stats;
                    end
                end
            end
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



