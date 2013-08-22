function plot_all_psychometric_neurometric3(cutoff,excludeSuppressed,normalised,V1only,V4only)
%Written by Xing 21/08/13
%Compare neurometric and psychometric function parameters. Convert data
%into z-scores.
%Set V1only to 1 to just combine data across V1 sites, set to 0 to include
%V4 data.
onExternalHD=0;
if onExternalHD==1
    rootFolder='K:\PL_backup_190713';
else
    rootFolder='F:';
end
newColCode=1;
fitCurves=1;
separateAnimalFigs=1;
figGauss=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.45, 0.65]); %
set(figGauss, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 3.305 3.325/0.4*0.45 3.305/5*65]);
animals=[{'blanco'} {'jack'}];
if V1only==0&&V4only==0
    numAreas=5;
elseif V1only==0&&V4only==1
    numAreas=1;
elseif V1only==1&&V4only==0
    numAreas=4;
end
allCols=cell(2,numAreas);
PROBMATslopeNeuroNew=cell(2,numAreas);%store all PROBMAT values
PROBMATPNENew=cell(2,numAreas);
PROBMATminRateNew=cell(2,numAreas);
PROBMATmaxRateNew=cell(2,numAreas);
CRFslopeNeuroNew=cell(2,numAreas);%store all CRF values
CRFC50New=cell(2,numAreas);
CRFminRateNew=cell(2,numAreas);
CRFmaxRateNew=cell(2,numAreas);
perfAll=cell(2,numAreas);%store all psychometric perf values
slopeAll=cell(2,numAreas);%store all psychometric slope values
PSEAll=cell(2,numAreas);%store all psychometric PSE values
areaCount=0;
for roving=0:1    
    if roving==0
        if V1only==1
            areas=[{'v1_1'}];
        elseif V4only==1
            areas=[{'v4_1'}];
        else
            areas=[{'v4_1'} {'v1_1'}];
        end
        allTableStats=[];
        allParTableStats=[];
        mpallTableStats=[];
    elseif roving==1
        areas=[{'v1_2_1'} {'v1_2_2'} {'v1_2_3'}];
        if V4only==1
            areas=[];
        end
        rovingCol=[0.5 0 0;0.75 0 0;1 0 0];
        rovingCol=[1 0 0;10/255 170/255 60/255;0 0 1];
        if newColCode==1
            taskCols=[204/256 153/256 256/256;242/256 174/256 15/256;154/256 205/256 50/256];%purple;orange;green
            taskColsDark=[153/256 50/256 204/256;179/256 129/256 14/256;34/256 139/256 34/256];%purple;orange;green
            taskColsVDark=[84/256 3/256 163/256;135/256 97/256 9/256;79/256 117/256 2/256];%purple;orange;green, for fitted line
        end
        allTableStats=cell(1,3);
        allParTableStats=cell(1,3);
        mpallTableStats=cell(1,3);
    end
    for areaInd=1:length(areas)
        areaCount=areaCount+1;
        area=areas{areaInd};
        [sampleContrasts testContrasts]=area_metadata(area);
        if roving==1&&newColCode==1
            taskCol=taskCols(areaInd,:);
            taskColDark=taskColsDark(areaInd,:);
            markerFaceTypes=[{'none'} {taskCol} {taskColDark}];
        end
        for animalInd=1:length(animals)
            animal=animals{animalInd};
            for sampleContrastsInd=1:length(sampleContrasts)
                sampleContrast=sampleContrasts(sampleContrastsInd);
                if roving==1&&newColCode==1
                    if sampleContrast==30
                        markerFaceType=markerFaceTypes{2};
                        fittedLineCol=taskColsDark(areaInd,:);
                    elseif sampleContrast==20
                        markerFaceType=markerFaceTypes{1};
                        fittedLineCol=taskCols(areaInd,:);
                    elseif sampleContrast==40
                        markerFaceType=markerFaceTypes{3};
                        fittedLineCol=taskColsVDark(areaInd,:);
                    end
                end
                if roving==0
                    if strcmp(area,'v4_1')
                        colMarker=[1 0 0];
                    elseif strcmp(area,'v1_1')
                        colMarker=[0 0 1];
                    end
                elseif roving==1
                    colMarker=fittedLineCol;
                end
                sessionNums=main_raw_sessions_final(animal,area,[],0);
                
                %read behavioural perf values
                loadText=['load ',rootFolder,'\PL\psycho_data\',animal,'\allMeanPerf\allMeanPerf_',area,'_',num2str(sampleContrast),'.mat allMeanPerf'];
                eval(loadText)
                perf=[];
                slope=[];
                PSE=[];
                for i=1:length(sessionNums)
                    rowInd=find(allMeanPerf(:,1)==sessionNums(i));
                    perf=[perf allMeanPerf(rowInd,2)];
                    if roving==1
                        slope=[slope allMeanPerf(rowInd,16)];
                        PSE=[PSE allMeanPerf(rowInd,15)];
                        %PSE=[PSE
                        %allMeanPerf(rowInd,15)-sampleContrast];%incorrect- do not use!
                    elseif roving==0
                        slope=[slope allMeanPerf(rowInd,18)];
                        PSE=[PSE allMeanPerf(rowInd,17)];
                        %PSE=[PSE allMeanPerf(rowInd,17)-sampleContrast];%incorrect
                    end
                end
                perfAll{animalInd,areaCount}=[perfAll{animalInd,areaCount} perf];
                slopeAll{animalInd,areaCount}=[slopeAll{animalInd,areaCount} slope];
                PSEAll{animalInd,areaCount}=[PSEAll{animalInd,areaCount} PSE];
                for colCount=1:length(perf)
                    allCols{animalInd,areaCount}=[allCols{animalInd,areaCount};colMarker];%keep track of what colour each data point should be
                end
                
                %read PROBMAT param values (slope, PNE, min, max):
                analysisType='ROC_zero_one';
                subFolder='new_vs_old_sglroc3acrosschannels';
                if excludeSuppressed==1
                    subFolder=[subFolder,'_excludeSuppressed'];
                end
                if normalised
                    subFolder=[subFolder,'_normalised'];
                end
                loadText=['load ',rootFolder,'\PL\',analysisType,'\',animal,'\',subFolder,'\cumulative_ROCs_old_new_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10),'.mat'];
                eval(loadText)
                PROBMATslopeNeuroNew{animalInd,areaCount}=[PROBMATslopeNeuroNew{animalInd,areaCount} slopeNeuroNew];
                PROBMATPNENew{animalInd,areaCount}=[PROBMATPNENew{animalInd,areaCount} PNENew];
                PROBMATminRateNew{animalInd,areaCount}=[PROBMATminRateNew{animalInd,areaCount} minRateNew];
                PROBMATmaxRateNew{animalInd,areaCount}=[PROBMATmaxRateNew{animalInd,areaCount} maxRateNew];
                if length(allMeanPerf(:,2))~=length(slopeNeuroNew)
                    unequal=1;
                end
                
                %read CRF param values (slope, C50, min, max):
                analysisType='CRF';
                subFolder='crf_meanchannels';
                loadText=['load ',rootFolder,'\PL\',analysisType,'\',animal,'\',subFolder,'\cumulative_CRFs_old_new_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10),'.mat'];
                eval(loadText)
                CRFslopeNeuroNew{animalInd,areaCount}=[CRFslopeNeuroNew{animalInd,areaCount} slopeNeuroNew];
                CRFC50New{animalInd,areaCount}=[CRFC50New{animalInd,areaCount} PNENew];
                CRFminRateNew{animalInd,areaCount}=[CRFminRateNew{animalInd,areaCount} minRateNew];
                CRFmaxRateNew{animalInd,areaCount}=[CRFmaxRateNew{animalInd,areaCount} maxRateNew];
                if length(allMeanPerf(:,2))~=length(slopeNeuroNew)
                    unequal=1;
                end
            end
        end
    end
end
%convert data into z-scores:
allzPROBMATslopeNeuroNew=cell(1,2);
allzPROBMATPNENew=cell(1,2);
allzCRFslopeNeuroNew=cell(1,2);
allzCRFC50New=cell(1,2);
allzperfAll=cell(1,2);
allzslopeAll=cell(1,2);
allzPSEAll=cell(1,2);
allzPROBMATminRateNew=cell(1,2);
allzPROBMATmaxRateNew=cell(1,2);
allzCRFminRateNew=cell(1,2);
allzCRFmaxRateNew=cell(1,2);
for animalInd=1:2   
    for areaCount=1:size(PROBMATslopeNeuroNew,2)
        zPROBMATslopeNeuroNew{animalInd,areaCount}=(PROBMATslopeNeuroNew{animalInd,areaCount}-mean(PROBMATslopeNeuroNew{animalInd,areaCount}))/std(PROBMATslopeNeuroNew{animalInd,areaCount});
        zPROBMATPNENew{animalInd,areaCount}=(PROBMATPNENew{animalInd,areaCount}-mean(PROBMATPNENew{animalInd,areaCount}))/std(PROBMATPNENew{animalInd,areaCount});
        zPROBMATminRateNew{animalInd,areaCount}=(PROBMATminRateNew{animalInd,areaCount}-mean(PROBMATminRateNew{animalInd,areaCount}))/std(PROBMATminRateNew{animalInd,areaCount});
        zPROBMATmaxRateNew{animalInd,areaCount}=(PROBMATmaxRateNew{animalInd,areaCount}-mean(PROBMATmaxRateNew{animalInd,areaCount}))/std(PROBMATmaxRateNew{animalInd,areaCount});
        
        zCRFslopeNeuroNew{animalInd,areaCount}=(CRFslopeNeuroNew{animalInd,areaCount}-mean(CRFslopeNeuroNew{animalInd,areaCount}))/std(CRFslopeNeuroNew{animalInd,areaCount});
        zCRFC50New{animalInd,areaCount}=(CRFC50New{animalInd,areaCount}-mean(CRFC50New{animalInd,areaCount}))/std(CRFC50New{animalInd,areaCount});
        zCRFminRateNew{animalInd,areaCount}=(CRFminRateNew{animalInd,areaCount}-mean(CRFminRateNew{animalInd,areaCount}))/std(CRFminRateNew{animalInd,areaCount});
        zCRFmaxRateNew{animalInd,areaCount}=(CRFmaxRateNew{animalInd,areaCount}-mean(CRFmaxRateNew{animalInd,areaCount}))/std(CRFmaxRateNew{animalInd,areaCount});
        
        zperfAll{animalInd,areaCount}=(perfAll{animalInd,areaCount}-mean(perfAll{animalInd,areaCount}))/std(perfAll{animalInd,areaCount});
        zslopeAll{animalInd,areaCount}=(slopeAll{animalInd,areaCount}-mean(slopeAll{animalInd,areaCount}))/std(slopeAll{animalInd,areaCount});
        zPSEAll{animalInd,areaCount}=(PSEAll{animalInd,areaCount}-mean(PSEAll{animalInd,areaCount}))/std(PSEAll{animalInd,areaCount});
        for sessionCount=1:length(zPROBMATslopeNeuroNew{animalInd,areaCount})
            %plot PROBMAT param values:
            subplot(3,4,animalInd+2);
            plot(zslopeAll{animalInd,areaCount}(sessionCount),zPROBMATslopeNeuroNew{animalInd,areaCount}(sessionCount),'Marker','o','Color',allCols{animalInd,areaCount}(sessionCount,:),'LineStyle','none');hold on
            subplot(3,4,4+animalInd+2);
            plot(zperfAll{animalInd,areaCount}(sessionCount),zPROBMATslopeNeuroNew{animalInd,areaCount}(sessionCount),'Marker','o','Color',allCols{animalInd,areaCount}(sessionCount,:),'LineStyle','none');hold on
            subplot(3,4,8+animalInd+2);
            plot(zPSEAll{animalInd,areaCount}(sessionCount),zPROBMATPNENew{animalInd,areaCount}(sessionCount),'Marker','o','Color',allCols{animalInd,areaCount}(sessionCount,:),'LineStyle','none');hold on
            
            %plot CRF param values:
            subplot(3,4,animalInd);
            plot(zslopeAll{animalInd,areaCount}(sessionCount),zCRFslopeNeuroNew{animalInd,areaCount}(sessionCount),'Marker','o','Color',allCols{animalInd,areaCount}(sessionCount,:),'LineStyle','none');hold on
            subplot(3,4,4+animalInd);
            plot(zperfAll{animalInd,areaCount}(sessionCount),zCRFslopeNeuroNew{animalInd,areaCount}(sessionCount),'Marker','o','Color',allCols{animalInd,areaCount}(sessionCount,:),'LineStyle','none');hold on
            subplot(3,4,8+animalInd);
            plot(zPSEAll{animalInd,areaCount}(sessionCount),zCRFC50New{animalInd,areaCount}(sessionCount),'Marker','o','Color',allCols{animalInd,areaCount}(sessionCount,:),'LineStyle','none');hold on
        end
        allzPROBMATslopeNeuroNew{animalInd}=[allzPROBMATslopeNeuroNew{animalInd} zPROBMATslopeNeuroNew{animalInd,areaCount}];
        allzPROBMATPNENew{animalInd}=[allzPROBMATPNENew{animalInd} zPROBMATPNENew{animalInd,areaCount}];
        allzCRFslopeNeuroNew{animalInd}=[allzCRFslopeNeuroNew{animalInd} zCRFslopeNeuroNew{animalInd,areaCount}];
        allzCRFC50New{animalInd}=[allzCRFC50New{animalInd} zCRFC50New{animalInd,areaCount}];
        allzperfAll{animalInd}=[allzperfAll{animalInd} zperfAll{animalInd,areaCount}];
        allzslopeAll{animalInd}=[allzslopeAll{animalInd} zslopeAll{animalInd,areaCount}];
        allzPSEAll{animalInd}=[allzPSEAll{animalInd} zPSEAll{animalInd,areaCount}];
        allzPROBMATminRateNew{animalInd}=[allzPROBMATminRateNew{animalInd} zPROBMATminRateNew{animalInd,areaCount}];
        allzPROBMATmaxRateNew{animalInd}=[allzPROBMATmaxRateNew{animalInd} zPROBMATmaxRateNew{animalInd,areaCount}];
        allzCRFminRateNew{animalInd}=[allzCRFminRateNew{animalInd} zCRFminRateNew{animalInd,areaCount}];
        allzCRFmaxRateNew{animalInd}=[allzCRFmaxRateNew{animalInd} zCRFmaxRateNew{animalInd,areaCount}];
    end
end

%calculate regular correlations
for animalInd=1:2
    [rho p]=corr(allzCRFslopeNeuroNew{animalInd}',allzslopeAll{animalInd}','type','Spearman')
    tableStats(1,[1:3]+3*(animalInd-1))=[length(allzCRFslopeNeuroNew{animalInd})-2 rho p];
    [rho p]=corr(allzPROBMATslopeNeuroNew{animalInd}',allzslopeAll{animalInd}','type','Spearman')
    tableStats(2,[1:3]+3*(animalInd-1))=[length(allzPROBMATslopeNeuroNew{animalInd})-2 rho p];
    [rho p]=corr(allzCRFslopeNeuroNew{animalInd}',allzperfAll{animalInd}','type','Spearman')
    tableStats(3,[1:3]+3*(animalInd-1))=[length(allzCRFslopeNeuroNew{animalInd})-2 rho p];
    [rho p]=corr(allzPROBMATslopeNeuroNew{animalInd}',allzperfAll{animalInd}','type','Spearman')
    tableStats(4,[1:3]+3*(animalInd-1))=[length(allzPROBMATslopeNeuroNew{animalInd})-2 rho p];
    [rho p]=corr(allzCRFC50New{animalInd}',allzPSEAll{animalInd}','type','Spearman')
    tableStats(5,[1:3]+3*(animalInd-1))=[length(allzCRFC50New{animalInd})-2 rho p];
    [rho p]=corr(allzPROBMATPNENew{animalInd}',allzPSEAll{animalInd}','type','Spearman')
    tableStats(6,[1:3]+3*(animalInd-1))=[length(allzPROBMATPNENew{animalInd})-2 rho p];
    [rho p]=corr(allzCRFminRateNew{animalInd}',allzperfAll{animalInd}','type','Spearman')
    tableStats(7,[1:3]+3*(animalInd-1))=[length(allzCRFminRateNew{animalInd})-2 rho p];
    [rho p]=corr(allzPROBMATminRateNew{animalInd}',allzperfAll{animalInd}','type','Spearman')
    tableStats(8,[1:3]+3*(animalInd-1))=[length(allzPROBMATminRateNew{animalInd})-2 rho p];
    [rho p]=corr(allzCRFmaxRateNew{animalInd}',allzperfAll{animalInd}','type','Spearman')
    tableStats(9,[1:3]+3*(animalInd-1))=[length(allzCRFmaxRateNew{animalInd})-2 rho p];
    [rho p]=corr(allzPROBMATmaxRateNew{animalInd}',allzperfAll{animalInd}','type','Spearman')
    tableStats(10,[1:3]+3*(animalInd-1))=[length(allzPROBMATmaxRateNew{animalInd})-2 rho p];
end
subplot(3,4,1);
axis square
xlabel('Psychometric slope');
ylabel('CRF slope');
title('Monkey 1');
subplot(3,4,2);
axis square
title('Monkey 2');
subplot(3,4,3);
axis square
xlabel('Psychometric slope');
ylabel('PROBMAT slope');
title('Monkey 1');
subplot(3,4,4);
axis square
title('Monkey 2');
subplot(3,4,5);
axis square
xlabel('P_c_o_r_r_e_c_t');
ylabel('CRF slope');
subplot(3,4,6);
axis square
subplot(3,4,7);
axis square
xlabel('P_c_o_r_r_e_c_t');
ylabel('PROBMAT slope');
subplot(3,4,8);
axis square
subplot(3,4,9);
axis square
xlabel('PSE');
ylabel('C_5_0');
% xlabel('DI_P_S_E');
% ylabel('DI_C_5_0');
subplot(3,4,10);
axis square
subplot(3,4,11);
axis square
xlabel('PSE');
ylabel('PNE');
% xlabel('DI_P_S_E');
% ylabel('DI_P_N_E');
subplot(3,4,12);
axis square

subplot(3,4,1);
xlim([-3 4]);
ylim([-3 5]);
plot([-5 5],[0 0],'k--');
plot([0 0],[-5 5],'k--');
subplot(3,4,2);
xlim([-3 4]);
ylim([-3 3]);
plot([-5 5],[0 0],'k--');
plot([0 0],[-5 5],'k--');
subplot(3,4,3);
xlim([-3 3]);
ylim([-4 4]);
plot([-5 5],[0 0],'k--');
plot([0 0],[-5 5],'k--');
subplot(3,4,4);
xlim([-3 4]);
ylim([-2 3]);
plot([-5 5],[0 0],'k--');
plot([0 0],[-5 5],'k--');
subplot(3,4,5);
xlim([-4 4]);
ylim([-2 3]);
plot([-5 5],[0 0],'k--');
plot([0 0],[-5 5],'k--');
subplot(3,4,6);
xlim([-4 4]);
ylim([-3 3]);
plot([-5 5],[0 0],'k--');
plot([0 0],[-5 5],'k--');
subplot(3,4,7);
xlim([-4 4]);
ylim([-3 4]);
plot([-5 5],[0 0],'k--');
plot([0 0],[-5 5],'k--');
subplot(3,4,8);
xlim([-3 3]);
ylim([-2 3]);
plot([-5 5],[0 0],'k--');
plot([0 0],[-5 5],'k--');
subplot(3,4,9);
xlim([-4 3]);
ylim([-3 2]);
plot([-5 5],[0 0],'k--');
plot([0 0],[-5 5],'k--');
subplot(3,4,10);
xlim([-3 4]);
ylim([-3 3]);
plot([-5 5],[0 0],'k--');
plot([0 0],[-5 5],'k--');
subplot(3,4,11);
xlim([-4 3]);
ylim([-4 3]);
plot([-5 5],[0 0],'k--');
plot([0 0],[-5 5],'k--');
subplot(3,4,12);
xlim([-4 4]);
ylim([-2 3]);
plot([-5 5],[0 0],'k--');
plot([0 0],[-5 5],'k--');

imagename=['bj_all_z_neuro_psycho_ROC_4_params_sessions_cutoff',num2str(cutoff*10)];
if excludeSuppressed==1
    imagename=['bj_all_neuro_psycho_ROC_4_params_sessions_cutoff',num2str(cutoff*10),'_excludeSuppressed'];
end
pathname=fullfile(rootFolder,'PL',imagename);
printtext=sprintf('print -dpng %s.png',pathname);
set(gcf,'PaperPositionMode','auto')
eval(printtext);
cutoff
allTableStats