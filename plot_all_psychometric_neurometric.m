function plot_all_psychometric_neurometric(cutoff,excludeSuppressed,normalised,V1only)
%Written by Xing 04/08/13
%Compare neurometric and psychometric function parameters.
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
figGauss=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.45, 0.65]); %
set(figGauss, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 3.305 3.325/0.4*0.45 3.305/5*65]);
animals=[{'blanco'} {'jack'}];
PROBMATslopeNeuroNew=cell(1,2);%store all PROBMAT values
PROBMATPNENew=cell(1,2);
PROBMATminRateNew=cell(1,2);
PROBMATmaxRateNew=cell(1,2);
CRFslopeNeuroNew=cell(1,2);%store all CRF values
CRFC50New=cell(1,2);
CRFminRateNew=cell(1,2);
CRFmaxRateNew=cell(1,2);
perfAll=cell(1,2);%store all psychometric perf values
slopeAll=cell(1,2);%store all psychometric slope values
PSEAll=cell(1,2);%store all psychometric PSE values
for roving=0:1    
    if roving==0
        if V1only==1
            areas=[{'v1_1'}];
        else
            areas=[{'v4_1'} {'v1_1'}];
        end
        allTableStats=[];
        allParTableStats=[];
        mpallTableStats=[];
    elseif roving==1
        areas=[{'v1_2_1'} {'v1_2_2'}];
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
                        colMarker='r';
                    elseif strcmp(area,'v1_1')
                        colMarker='b';
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
                        PSE=[PSE allMeanPerf(rowInd,15)-sampleContrast];
                    elseif roving==0
                        slope=[slope allMeanPerf(rowInd,18)];
                        PSE=[PSE allMeanPerf(rowInd,17)-sampleContrast];
                    end
                end
                perfAll{animalInd}=[perfAll{animalInd} perf];
                slopeAll{animalInd}=[slopeAll{animalInd} slope];
                PSEAll{animalInd}=[PSEAll{animalInd} PSE];
                
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
                PROBMATslopeNeuroNew{animalInd}=[PROBMATslopeNeuroNew{animalInd} slopeNeuroNew];
                PNENew=PNENew-sampleContrast;
                PROBMATPNENew{animalInd}=[PROBMATPNENew{animalInd} PNENew];
                PROBMATminRateNew{animalInd}=[PROBMATminRateNew{animalInd} minRateNew];
                PROBMATmaxRateNew{animalInd}=[PROBMATmaxRateNew{animalInd} maxRateNew];
                if length(allMeanPerf(:,2))~=length(slopeNeuroNew)
                    unequal=1;
                end
                subplot(4,4,animalInd+2);
                plot(perf,slopeNeuroNew,'Marker','o','Color',colMarker,'LineStyle','none');hold on
                subplot(4,4,4+animalInd+2);
                plot(PSE,PNENew,'Marker','o','Color',colMarker,'LineStyle','none');hold on
                subplot(4,4,8+animalInd+2);
                plot(perf,minRateNew,'Marker','o','Color',colMarker,'LineStyle','none');hold on
                subplot(4,4,12+animalInd+2);
                plot(perf,maxRateNew,'Marker','o','Color',colMarker,'LineStyle','none');hold on
                
                %read CRF param values (slope, C50, min, max):
                analysisType='CRF';
                subFolder='crf_meanchannels';
                loadText=['load ',rootFolder,'\PL\',analysisType,'\',animal,'\',subFolder,'\cumulative_CRFs_old_new_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10),'.mat'];
                eval(loadText)
                CRFslopeNeuroNew{animalInd}=[CRFslopeNeuroNew{animalInd} slopeNeuroNew];
                PNENew=PNENew-sampleContrast;
                CRFC50New{animalInd}=[CRFC50New{animalInd} PNENew];
                CRFminRateNew{animalInd}=[CRFminRateNew{animalInd} minRateNew];
                CRFmaxRateNew{animalInd}=[CRFmaxRateNew{animalInd} maxRateNew];
                if length(allMeanPerf(:,2))~=length(slopeNeuroNew)
                    unequal=1;
                end
                subplot(4,4,animalInd);
                plot(perf,slopeNeuroNew,'Marker','o','Color',colMarker,'LineStyle','none');hold on
                subplot(4,4,4+animalInd);
                plot(PSE,PNENew,'Marker','o','Color',colMarker,'LineStyle','none');hold on
                subplot(4,4,8+animalInd);
                plot(perf,minRateNew,'Marker','o','Color',colMarker,'LineStyle','none');hold on
                subplot(4,4,12+animalInd);
                plot(perf,maxRateNew,'Marker','o','Color',colMarker,'LineStyle','none');hold on
                
            end
        end
    end
end
%calculate regular correlations
for animalInd=1:2
    [rho p]=corr(CRFslopeNeuroNew{animalInd}',slopeAll{animalInd}','type','Spearman')
    tableStats(1+4*(animalInd-1),1:3)=[length(CRFslopeNeuroNew{animalInd})-2 rho p];
    [rho p]=corr(PROBMATslopeNeuroNew{animalInd}',slopeAll{animalInd}','type','Spearman')
    tableStats(1+4*(animalInd-1),4:6)=[length(PROBMATslopeNeuroNew{animalInd})-2 rho p];
    [rho p]=corr(CRFslopeNeuroNew{animalInd}',perfAll{animalInd}','type','Spearman')
    tableStats(1+4*(animalInd-1),7:9)=[length(CRFslopeNeuroNew{animalInd})-2 rho p];
    [rho p]=corr(PROBMATslopeNeuroNew{animalInd}',perfAll{animalInd}','type','Spearman')
    tableStats(1+4*(animalInd-1),10:12)=[length(PROBMATslopeNeuroNew{animalInd})-2 rho p];
    [rho p]=corr(CRFC50New{animalInd}',PSEAll{animalInd}','type','Spearman')
    tableStats(2+4*(animalInd-1),1:3)=[length(CRFC50New{animalInd})-2 rho p];
    [rho p]=corr(PROBMATPNENew{animalInd}',PSEAll{animalInd}','type','Spearman')
    tableStats(2+4*(animalInd-1),4:6)=[length(PROBMATPNENew{animalInd})-2 rho p];
    [rho p]=corr(CRFminRateNew{animalInd}',perfAll{animalInd}','type','Spearman')
    tableStats(3+4*(animalInd-1),1:3)=[length(CRFminRateNew{animalInd})-2 rho p];
    [rho p]=corr(PROBMATminRateNew{animalInd}',perfAll{animalInd}','type','Spearman')
    tableStats(3+4*(animalInd-1),4:6)=[length(PROBMATminRateNew{animalInd})-2 rho p];
    [rho p]=corr(CRFmaxRateNew{animalInd}',perfAll{animalInd}','type','Spearman')
    tableStats(4+4*(animalInd-1),1:3)=[length(CRFmaxRateNew{animalInd})-2 rho p];
    [rho p]=corr(PROBMATmaxRateNew{animalInd}',perfAll{animalInd}','type','Spearman')
    tableStats(4+4*(animalInd-1),4:6)=[length(PROBMATmaxRateNew{animalInd})-2 rho p];
end
subplot(4,4,1);
xlim([0.68 0.94]);
ylim([0 0.012]);
subplot(4,4,3);
xlim([0.68 0.94]);
ylim([0 0.05]);
subplot(4,4,5);
xlim([-12 18]);
ylim([-20 50]);
plot([-12 18],[0 0],'k--');
plot([0 0],[-20 50],'k--');
% xlim([-0.4 1]);
% ylim([-0.5 2.5]);
% plot([-0.4 1],[0 0],'k--');
% plot([0 0],[-0.5 2.5],'k--');
subplot(4,4,6);
xlim([-20 40]);
ylim([-10 40]);
plot([-20 40],[0 0],'k--');
plot([0 0],[-10 40],'k--');
% xlim([-0.5 1.5]);
% ylim([-0.5 2]);
% plot([-0.5 1.5],[0 0],'k--');
% plot([0 0],[-0.5 2],'k--');
subplot(4,4,7);
xlim([-15 18]);
ylim([-40 30]);
plot([-15 18],[0 0],'k--');
plot([0 0],[-40 30],'k--');
% xlim([-0.4 0.8]);
% ylim([-0.5 1.5]);
% plot([-0.4 0.8],[0 0],'k--');
% plot([0 0],[-0.5 1.5],'k--');
subplot(4,4,8);
xlim([-20 40]);
ylim([-5 20]);
plot([-20 40],[0 0],'k--');
plot([0 0],[-5 20],'k--');
% xlim([-0.5 2]);
% ylim([-0.2 1]);
% plot([-0.5 2],[0 0],'k--');
% plot([0 0],[-0.2 1],'k--');
subplot(4,4,9);
xlim([0.65 0.95]);
ylim([0 18]);
subplot(4,4,10);
ylim([0 15]);
subplot(4,4,11);
ylim([0 0.6]);
subplot(4,4,13);
xlim([0.65 0.95]);
subplot(4,4,14);
xlim([0.5 0.95]);
ylim([0 50]);
% if normalised==0
subplot(4,4,1);
xlabel('Psychometric slope');
ylabel('CRF slope');
title('Monkey 1');
subplot(4,4,2);
title('Monkey 2');
subplot(4,4,3);
xlabel('Psychometric slope');
ylabel('PROBMAT slope');
title('Monkey 1');
subplot(4,4,4);
title('Monkey 2');
subplot(4,4,5);
xlabel('DI_P_S_E');
ylabel('DI_C_5_0');
subplot(4,4,7);
xlabel('DI_P_S_E');
ylabel('DI_P_N_E');
subplot(4,4,9);
ylabel('CRF minimum');
subplot(4,4,11);
ylabel('PROBMAT minimum');
subplot(4,4,13);
ylabel('CRF maximum');
subplot(4,4,15);
ylabel('PROBMAT maximum');

imagename=['bj_all_neuro_psycho_ROC_4_params_sessions_cutoff',num2str(cutoff*10)];
if excludeSuppressed==1
    imagename=['bj_all_neuro_psycho_ROC_4_params_sessions_cutoff',num2str(cutoff*10),'_excludeSuppressed'];
end
pathname=fullfile(rootFolder,'PL',imagename);
printtext=sprintf('print -dpng %s.png',pathname);
set(gcf,'PaperPositionMode','auto')
eval(printtext);
cutoff
allTableStats