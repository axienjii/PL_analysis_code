function bj_conflict_conds
%To compare performance across sample contrasts during stimulus roving, for
%conditions where the identity of the target stimulus was critically
%dependent on the contrast level of the sample stimulus. I.e. when the
%sample contrast was 20%, the conditions that induced this conflict were
%those where the test contrasts were 22, 25 and 28%. When the sample contrast
%was 40%, the conditions that induced this conflict were those where the
%test contrasts were 32, 35 and 38%.      
onExternalHD=0;
plotRTerr=1;
roving=[0 1];
plotStage3=1;
animals=[{'blanco'} {'jack'}];
areas=[{'v1_1'};{'v1_2_1'}];
rows=4;
if plotRTerr==1
    rows=5;
end
if onExternalHD==1
    rootFolder='G:\PL_backup_060413';
else
    rootFolder='F:';
end
% fig=figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.2, 0.04, 0.75/(3-length(areas)), 0.5]);
% set(fig,'PaperUnits','centimeters','PaperType','A4','PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 31.5/(3-length(areas)) 15.78]);
sampleConflictConds=[{6:8} {4:9} {5:7}];20% sample:22 25 28; 30% sample: 22 25 28 32 35 38; 40% sample: 32 35 8
allConflictConds=[22 25	28 32 35 38];%union of conditions which are to be plotted
% colmapText=colormap(jet(size(allConflictConds,2)));
taskCols=[0.7 0.7 0.7;173/256 216/255 230/256;204/256 153/256 256/256;207/256 253/256 162/256;255/256 222/256 173/256;255/256 153/256 153/256];%lighter versions black; blue;purple;green;orange;red
taskColsDark=[0 0 0;13/256 4/256 137/256;153/256 50/256 204/256;0/256 100/256 0/256;255/256 128/256 0/256;204/256 0 0];%black; blue;purple;green;orange;red
% colmapText=[colmapText(1,:);132/255 22/255 216/255;202/255 65/255 223/255;colmapText(3:7,:);157/255 212/255 61/255;colmapText(10:12,:);178/255 111/255 12/255;colmapText(end,:)];
for animalInd=1:length(animals)
    animal=animals{animalInd};
    fig(animalInd)=figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.2, 0.04, 0.5, 0.5]);
    set(fig,'PaperUnits','centimeters','PaperType','A4','PaperOrientation', 'landscape');
    for rovingInd=1:2%do non-roving first, then roving
        %V1_2_1
        area=areas{rovingInd,1};
        [sampleContrasts testContrasts]=area_metadata(area);
        for sampleContrastInd=1:length(sampleContrasts)
            sampleContrast=sampleContrasts(sampleContrastInd);
            if rovingInd==1
                conflictConds=5:9;% for non-roving, 30% sample. corresponds to test contrasts of [22 25 28 32 35]
            elseif rovingInd==2
                conflictConds=sampleConflictConds{sampleContrastInd};
            end
            hold all
            loadText=['load ',rootFolder,'\PL\psycho_data\',animal,'\allMeanPerf\allMeanPerf_',area,'_',num2str(sampleContrast),'.mat allMeanPerf'];
            eval(loadText)
            if rovingInd==1
                maxV1_1(animalInd)=size(allMeanPerf,1);
            end
            % boxplot([allMeanPerf(1:numN,2) allMeanPerf(end-numN+1:end,2)]);hold on
%             subplot(1,2,animalInd);
            for conflictCondInd=1:length(conflictConds)
                conflictCond=conflictConds(conflictCondInd);
                testContrast=testContrasts(sampleContrastInd,conflictCond);
                colInd=find(testContrast==allConflictConds);
                subplot(2,3,colInd);
                if sampleContrast==30
                    markerFaceType=taskCols(colInd,:);
                elseif sampleContrast==20
                    markerFaceType='none';
                elseif sampleContrast==40
                    markerFaceType=taskColsDark(colInd,:);
                end
                markerText='o';
                markerS=8;
                if strcmp(area,'v1_2_1')
                    x1vals=3+maxV1_1(animalInd);
                    x2vals=2+size(allMeanPerf,1)+maxV1_1(animalInd);
                else
                    if strcmp(animal,'blanco')&&strcmp(area,'v1_1')
                        xlim([0 size(allMeanPerf,1)+maxV1_1(animalInd)+4]);
                    elseif strcmp(animal,'jack')&&strcmp(area,'v1_2_1')
                        yLimVals=get(gca,'ylim');
                        xLimVals=get(gca,'xlim');
                        unitSpace=(yLimVals(2)-yLimVals(1))/30;
                        text('Position',[xLimVals(2)+(xLimVals(2)-xLimVals(1))/25 yLimVals(1)+unitSpace*conflictCond*2],'FontSize',9,'String',[markerText,'  ',num2str(testContrast(conflictCond)),'%'],'Color',colmapText(colInd,:));
                    end
                    x1vals=1;
                    x2vals=size(allMeanPerf,1);
                end
                yvals=allMeanPerf(:,2+conflictCond)*100;
                if rovingInd==1
                    allPerf(1,colInd)={yvals};%for each conflict condition, record the performance
                    allContrasts(1,colInd)=testContrast;%for each conflict condition, record the performance
                end
                if rovingInd==2
                    match=find(allContrasts(1,:)==testContrast);%if test contrast was presented during both non-roving and roving tasks, perform comparison
                    if ~isempty(match)
                        [h,p,ci,stats]=ttest2(allPerf{1,match},yvals);
                        pAll{animalInd}(sampleContrastInd,conflictCondInd)=p;
                    end
                end
                plot(x1vals:x2vals,yvals,'Color',taskColsDark(colInd,:),'LineStyle','none','Marker',markerText,'MarkerFaceColor',markerFaceType,'MarkerEdgeColor',taskColsDark(colInd,:),'MarkerSize',markerS);hold on%'MarkerFaceColor',[1/i 1/i 1/i],
                xplot=(x1vals:x2vals)';
                xTemp=xplot;
                yplot=yvals;
                %linear fit:
                if rovingInd==1
                    startpoint5=0;%1:start at 0.5
                else
                    startpoint5=0;%2:start at the first value
                end
                seedvalslin=[1 1];
                if startpoint5==0
                    f = fittype('a+b*x');
                elseif startpoint5==1
                    f = fittype('50+b*x');
                    seedvalslin=1;
                else
                    f = fittype(sprintf('%d+b*x',allMeanPerf(1,2+conflictCond)*100));
                    seedvalslin=1;
                end
                fit3 = fit(xTemp,yvals,f,'StartPoint',seedvalslin,'Robust','on');
                % plot(fit3,'b:');
                [clinear,goflinear] = fit(xplot,yplot,fit3)
                coefperflinear=coeffvalues(clinear);
                if startpoint5==1
                    coefperflinear=[50 coefperflinear];
                elseif startpoint5==2
                    coefperflinear=[allMeanPerf(1,2+conflictCond)*100 coefperflinear];
                end
                if rovingInd==1
                    xcurvevals=0:0.01:size(allMeanPerf(:,2+conflictCond)*100,1);
                else
                    xcurvevals=3+maxV1_1(animalInd):0.01:2+size(allMeanPerf(:,2+conflictCond)*100,1)+maxV1_1(animalInd);
                end
                ycurvevals=coefperflinear(1)+coefperflinear(2)*xcurvevals;
                plot(xcurvevals,ycurvevals,'Color',taskColsDark(colInd,:));
                if rovingInd==2
                    allSlopes{animalInd}(sampleContrastInd,colInd)=coefperflinear(2);
                end
                if animalInd==1
                    xlim([0 60]);
                else
                    xlim([0 44]);
                end
                ylim([0 100]);
                set(gca,'YTick',[0 50 100],'YTickLabel',[0 0.5 1]);
                title([num2str(testContrast),' %'],'FontSize',18);
                %     xLimVals=get(gca,'xlim');
                %     unitSpace=(yLimVals(2)-yLimVals(1))/30;
            end
        end
    end
end
printFigs=0;
if printFigs==1
    folderPrint=fullfile(rootFolder,'PL','psycho_data','behavioural_figures');
    formats=[{'epsc'} {'png'}];
    for i=1:length(formats)
        format=formats{i};
        for animalInd=1:length(animals)
            for areaInd=1:length(areas)
                figure(animalInd);
                figName=[folderPrint,'\',animals{animalInd},'_perf_conflict_conds'];
                printtext=sprintf('print -d%s %s -r600',format,figName);
                set(gcf, 'PaperPositionMode', 'auto');
                eval(printtext);
            end
        end
    end
end
close all hidden
