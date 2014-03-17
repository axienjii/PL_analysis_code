function bj_conflict_conds
%To compare performance across sample contrasts during stimulus roving, for
%conditions where the identity of the target stimulus was critically
%dependent on the contrast level of the sample stimulus. I.e. when the
%sample contrast was 20%, the conditions that induced this conflict were
%those where the test contrasts were 22, 25 and 28%. When the sample contrast
%was 40%, the conditions that induced this conflict were those where the
%test contrasts were 32, 35 and 38%.   
monochrome=1;
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
AICvals=[];
% fig=figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.2, 0.04, 0.75/(3-length(areas)), 0.5]);
% set(fig,'PaperUnits','centimeters','PaperType','A4','PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 31.5/(3-length(areas)) 15.78]);
%sampleConflictConds=[{6:8} {4:9} {5:7}];%20% sample:22 25 28; 30% sample: 22 25 28 32 35 38; 40% sample: 32 35 38
allConflictConds=[22 25	28 32 35 38];%union of conditions which are to be plotted
% colmapText=colormap(jet(size(allConflictConds,2)));
taskCols=[0.7 0.7 0.7;173/256 216/255 230/256;204/256 153/256 256/256;207/256 253/256 162/256;255/256 222/256 173/256;255/256 153/256 153/256];%lighter versions black; blue;purple;green;orange;red
taskColsDark=[0 0 0;13/256 4/256 137/256;153/256 50/256 204/256;0/256 100/256 0/256;255/256 128/256 0/256;204/256 0 0];%black; blue;purple;green;orange;red
% colmapText=[colmapText(1,:);132/255 22/255 216/255;202/255 65/255 223/255;colmapText(3:7,:);157/255 212/255 61/255;colmapText(10:12,:);178/255 111/255 12/255;colmapText(end,:)];
for animalInd=1:length(animals)
    animal=animals{animalInd};
    fig(animalInd)=figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.2, 0.04, 0.5, 0.5]);
    set(fig,'PaperUnits','centimeters','PaperType','A4','PaperOrientation', 'landscape');
    allPerfRovingInd=0;
    allPerfRoving=[];
    for rovingInd=1:2%do non-roving first, then roving
        %V1_2_1
        area=areas{rovingInd,1};
        [sampleContrasts testContrasts]=area_metadata(area);
        if rovingInd==1
            subplotsAll=1:5;
            conflictConds=5:9;% for non-roving, 30% sample. corresponds to test contrasts of [22 25 28 32 35]
            sampleContrastsConflict=[{1} {1} {1} {1} {1} {1}];
        elseif rovingInd==2
            subplotsAll=1:6;
            conflictConds=5:10;
            sampleContrastsConflict=[{[1 2]} {[1 2]} {[1 2]} {[2 3]} {[2 3]} {[2 3]}];%20% sample:22 25 28; 30% sample: 22 25 28 32 35 38; 40% sample: 32 35 38
        end
        for conflictCondInd=1:length(allConflictConds)
            if find(subplotsAll==conflictCondInd)
                yvalsColumns=[];
                yvalsCombined=[];
                nCombined=[];
                chisqSeparate=[];residualsSeparate=[];xmleSeparate=[];
                mleSeparate=0;AICcombined=[];
                chisqSeparate=0;
                sampleContrastConflict=sampleContrasts(sampleContrastsConflict{conflictCondInd});
                yFittedvalsColumns=[];
                for sampleContrastInd=sampleContrastsConflict{conflictCondInd}(1):sampleContrastsConflict{conflictCondInd}(end)
                    subplot(2,3,conflictCondInd);
                    testContrast=allConflictConds(conflictCondInd);
                    sampleContrast=sampleContrasts(sampleContrastInd);
                    conflictCond=find(testContrasts(sampleContrastInd,:)==testContrast);
                    hold all
                    loadText=['load ',rootFolder,'\PL\psycho_data\',animal,'\allMeanPerf\allMeanPerf_',area,'_',num2str(sampleContrast),'.mat allMeanPerf'];
                    eval(loadText)
                    if rovingInd==1
                        maxV1_1(animalInd)=size(allMeanPerf,1);
                    end
                    % boxplot([allMeanPerf(1:numN,2) allMeanPerf(end-numN+1:end,2)]);hold on
                    %             subplot(1,2,animalInd);
                    colInd=find(testContrast==allConflictConds);
                    if sampleContrast==30
                        markerFaceType=taskCols(colInd,:);
                    elseif sampleContrast==20
                        markerFaceType='none';
                    elseif sampleContrast==40
                        markerFaceType=taskColsDark(colInd,:);
                    end
                    if monochrome==1
                        if sampleContrast==30
                            monochromeTaskCols='k';
                        else
                            monochromeTaskCols='r';
                            monochromeTaskCols=[0.7 0.7 0.7];
                        end
                    end
                    markerText='o';
                    markerS=6;
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
                        allContrasts(1,colInd)=testContrast;%for each conflict condition, record the contrast
                    end
                    if rovingInd==2
                        match=find(allContrasts(1,:)==testContrast);%if test contrast was presented during both non-roving and roving tasks, perform comparison
                        if ~isempty(match)
                            [h,p,ci,stats]=ttest2(allPerf{1,match},yvals);
                            pAll{animalInd}(sampleContrastInd,conflictCondInd)=p;
                        end
                        allPerfRovingInd=allPerfRovingInd+1;
                        if testContrast>sampleContrast
                            allPerfRoving(allPerfRovingInd,:)=yvals;%store perf vals for roving period to plot mean Pcorrect later
                        elseif testContrast<sampleContrast
                            allPerfRoving(allPerfRovingInd,:)=100-yvals;%store perf vals for roving period to plot mean Pcorrect later
                        end
                        yvalsCombined=[yvalsCombined;yvals];%store values across sample contrast conditions
                        yvalsColumns=[yvalsColumns yvals];%store values across sample contrast conditions, but in different columns
                    end
                    %colour/monochrome version:
                    if monochrome==0
                        plot(x1vals:x2vals,yvals,'Color',taskColsDark(colInd,:),'LineStyle','none','Marker',markerText,'MarkerFaceColor',markerFaceType,'MarkerEdgeColor',taskColsDark(colInd,:),'MarkerSize',markerS);hold on;alpha(0.3)%'MarkerFaceColor',[1/i 1/i 1/i],
                    elseif monochrome==1
                        plot(x1vals:x2vals,yvals,'Color',monochromeTaskCols,'LineStyle','none','Marker',markerText,'MarkerFaceColor',monochromeTaskCols,'MarkerEdgeColor',monochromeTaskCols,'MarkerSize',4);hold on;%'MarkerFaceColor',[1/i 1/i 1/i],
                    end
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
                    yFittedvals=coefperflinear(1)+coefperflinear(2)*xplot;
                    yFittedvalsColumns=[yFittedvalsColumns coefperflinear(1)+coefperflinear(2)*xplot];%store fitted values in separate columns
                    if monochrome==0
                        plot(xcurvevals,ycurvevals,'Color',taskColsDark(colInd,:));
                    elseif monochrome==1
                        plot(xcurvevals,ycurvevals,'Color',monochromeTaskCols,'LineWidth',1.5);
                    end
                    %calculate MLE and gof statistics:
                    if rovingInd==2%only do for roving data, not non-roving
                        matPath=[rootFolder,'\PL\psycho_data\',animal,'\allNumTrials_v1_2_',num2str(sampleContrast),'.mat'];
                        if ~exist(matPath,'file')
                            pause
                        else
                            loadText=['load ',rootFolder,'\PL\psycho_data\',animal,'\allNumTrials_v1_2_',num2str(sampleContrast),'.mat perfAll'];
                            eval(loadText)
                            matchSessions=[];
                            for rowInd=1:size(allMeanPerf,1)
                                matchSessions=[matchSessions find(perfAll(:,1)==allMeanPerf(rowInd,1))];
                            end
                            nCorrect=perfAll(matchSessions,1+conflictCond);%perfAll contains number of correct trials
                            nAll=round(nCorrect.*100./yplot);%calculate total number of trials (correct plus incorrect). All values should end in '.0'- rounding is superfluous
                            nCombined=[nCombined;nAll];%store values across sample contrast conditions
                        end
                        [residual1, xmle1, chisq1]=Max_like_estimates_weibull_xing(nAll,yplot./100,yFittedvals./100);%convert percentages into fractions
                        residualsSeparate(sampleContrastInd)=residual1;
%                         xmleSeparate(sampleContrastInd)=xmle1;
%                         tempMLE=mle_calculator(nAll,yplot./100,yFittedvals./100);
%                         mleSeparate=mleSeparate+tempMLE;
%                         chisqSeparate(sampleContrastInd)=chisq1;
                        tempChisq=chisquare_calculator(yplot,yFittedvals);
                        chisqSeparate=chisqSeparate+tempChisq;
                        goodnessOfFit(sampleContrastInd,conflictCondInd)=gof_xing(yFittedvals,yplot);
%                         pFit(sampleContrastInd,conflictCondInd)=chi2cdf(goodnessOfFit(sampleContrastInd,conflictCondInd),length(matchSessions)-1);
                    end
                    set(gca, 'box', 'off')
                    if rovingInd==2
                        allSlopes{animalInd}(sampleContrastInd,colInd)=coefperflinear(2);
                        allIntercepts{animalInd}(sampleContrastInd,colInd)=coefperflinear(1);
                    end
                    if animalInd==1
                        xlim([0 60]);
                    else
                        xlim([0 44]);
                    end
                    ylim([0 100]);
                    set(gca,'YTick',[0 50 100],'YTickLabel',[0 0.5 1]);
                    %     xLimVals=get(gca,'xlim');
                    %     unitSpace=(yLimVals(2)-yLimVals(1))/30;
                end
            end
            title([num2str(testContrast),' %'],'FontSize',18);
            if rovingInd==2
                %calculate residuals when data are combined across sample
                %contrast conditions:
                xvalsCombined=[xplot;xplot];
                %linear fit:
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
                fit3 = fit(xvalsCombined,yvalsCombined,f,'StartPoint',seedvalslin,'Robust','on');
                % plot(fit3,'b:');
                [clinear,goflinear] = fit(xvalsCombined,yvalsCombined,fit3)
                coefperflinear=coeffvalues(clinear);
                if startpoint5==1
                    coefperflinear=[50 coefperflinear];
                elseif startpoint5==2
                    coefperflinear=[allMeanPerf(1,2+conflictCond)*100 coefperflinear];
                end
                yFittedvals=coefperflinear(1)+coefperflinear(2)*xvalsCombined;
                %calculate MLE and gof statistics:
                [residualsCombined, xmle2, chisq2]=Max_like_estimates_weibull_xing(nCombined,yvalsCombined./100,yFittedvals./100);%convert percentages into fractions
                goodnessOfFitCombined(conflictCondInd)=gof_xing(yFittedvals,yvalsCombined);
%                 chisqCombined(conflictCondInd)=chisq2;
%                 mleCombined(conflictCondInd)=mle_calculator(nCombined,yvalsCombined./100,yFittedvals./100);
%                 AICcombined=2*2-2*mleCombined(conflictCondInd);%AIC=2*df-2*ln(L) = chi-squared value + 2*df. note that 'log' function in matlab takes natural log, ie.e 'ln'
%                 AICseparate=2*4-2*mleSeparate;
                chisqCombined(conflictCondInd)=chisquare_calculator(yvalsCombined,yFittedvals);
                AICcombined=chisqCombined(conflictCondInd)+2*2;%AIC= chi-squared value + 2*df. 
                AICseparate=chisqSeparate+2*4;
                AICvals(:,conflictCondInd)=[AICseparate;AICcombined];
%                 diffChi=sum(chisqSeparate)-chisq2;%difference between chi squared errors
%                 p_gof(animalInd,conflictCondInd)=1-chi2cdf(diffChi,2);
                %attempt to swap param values between 2 separate fits:
%                 goodnessOfFitSwap(1,conflictCondInd)=gof_xing(yFittedvalsColumns(:,2),yvalsColumns(:,1));%check gof for 2nd sample contrast data, based on 1st sample fitting
%                 goodnessOfFitSwap(2,conflictCondInd)=gof_xing(yFittedvalsColumns(:,1),yvalsColumns(:,2));%check gof for 1st sample contrast data, based on 2nd sample fitting
%                 pFitSwap(1,conflictCondInd)=chi2cdf(goodnessOfFitSwap(1),length(matchSessions)-1);
%                 pFitSwap(2,conflictCondInd)=chi2cdf(goodnessOfFitSwap(2),length(matchSessions)-1);
                %                 dfs=[4 2];%two free parameters: a and b
                summedResiduals=sum(residualsSeparate);%add together residuals from two separate linear fits
                compareResiduals(:,conflictCondInd)=[summedResiduals;residualsCombined];
%                 p_gof(conflictCondInd,:)=1-chi2cdf(compareResiduals,length(nCombined)-(1+dfs));%df equals to number of sessions minus 1 minus the number of free parameters involved
%                 xmleCombined(conflictCondInd)=xmle2;
%                 sum(xmleSeparate)                
            end
        end
    end
    allPerfRovingMonkeys{animalInd}=allPerfRoving;    
    goodnessOfFitSummed=sum(goodnessOfFit);
    goodnessOfFitDiff(animalInd,:)=goodnessOfFitCombined-goodnessOfFitSummed;%positive values indicate that the fit for combined data is worse than for separate data
    goodnessOfFitMonkeys{animalInd}=goodnessOfFit;
%     goodnessOfFitSwapMonkeys{animalInd}=goodnessOfFitSwap;
%     pFitMonkeys{animalInd}=pFit;
%     pSwapMonkeys{animalInd}=pFitSwap;
    AICMonkeys{animalInd}=AICvals;
%     diffAIC=AICvals(2,:)-AICvals(1,:);
    minAIC=min(AICvals);
    diffAIC2=AICvals(2,:)-minAIC;    
    diffAIC1=AICvals(1,:)-minAIC;    
    for conflictCondInd=1:6
        wAICMonkeys(animalInd,conflictCondInd)=exp(-0.5.*diffAIC1(conflictCondInd))./(exp(-0.5.*diffAIC1(conflictCondInd))+exp(-0.5.*diffAIC2(conflictCondInd)));%note that exp(0)=1
    end
%     diffAIC2=AICvals(1,:)-minAIC;    
%     for conflictCondInd=1:6
%         wAICMonkeys2{animalInd}(conflictCondInd)=exp(-0.5.*diffAIC2(conflictCondInd))./(exp(-0.5.*diffAIC2(conflictCondInd))+exp(0));%note that exp(0)=1
%     end
%     ratioWeights(animalInd,:)=wAICMonkeys2{animalInd}./(wAICMonkeys2{animalInd}+wAICMonkeys{animalInd});
end
figure
for conflictCondInd=1:length(allConflictConds)
%     plot([1 2],AICMonkeys{1,1}(:,conflictCondInd),'ko','LineStyle','-','Color','k');hold on
%     plot([3 4],AICMonkeys{1,2}(:,conflictCondInd),'ko','LineStyle','-','Color','k');
    for animalInd=1:2
        plot(AICMonkeys{1,animalInd}(1,conflictCondInd),AICMonkeys{1,animalInd}(2,conflictCondInd),'ko','LineStyle','-','Color','k');hold on
        AICdiff(animalInd,conflictCondInd)=AICMonkeys{1,animalInd}(2,conflictCondInd)-AICMonkeys{1,animalInd}(1,conflictCondInd);%find difference between AICseparate and AICcombined values
        AICweight(animalInd,conflictCondInd)=exp(-0.5*AICdiff(animalInd,conflictCondInd))/(exp(-0.5*AICdiff(animalInd,conflictCondInd))+1);%calculate weight value according to value in Mehdi's thesis, Ch 4 page 160
    end
end
fig2=figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.2, 0.04, 0.5, 0.5]);
set(fig2,'PaperUnits','centimeters','PaperType','A4','PaperOrientation', 'landscape');
for animalInd=1:length(animals)
    subplot(1,2,animalInd);
    for condInd=1:size(allPerfRovingMonkeys{animalInd},1)
        plot(1:length(allPerfRovingMonkeys{animalInd}(condInd,:)),allPerfRovingMonkeys{animalInd}(condInd,:),'k.');
        hold on        
    end
end
fig3=figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.2, 0.04, 0.17, 0.3]);
set(fig3,'PaperUnits','centimeters','PaperType','A4','PaperOrientation', 'landscape');
for animalInd=1:length(animals)
    for condInd=1:size(allPerfRovingMonkeys{animalInd},1)
        markerFill=[{'none'} {'k'}];
        plot(1:length(allPerfRovingMonkeys{animalInd}(condInd,:)),mean(allPerfRovingMonkeys{animalInd},1),'ko','MarkerFaceColor',markerFill{animalInd});
        hold on
    end
    plotLinear=0;
    startpoint5=0;
    bj_linearexpo_fitting(testContrast,mean(allPerfRovingMonkeys{animalInd},1)',8,startpoint5,'ROC',plotLinear,[],[0 0 0]);
end
xlim([0 35]);
xlabel('session number');
ylabel('Pcorrect for conflict conditions');
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
