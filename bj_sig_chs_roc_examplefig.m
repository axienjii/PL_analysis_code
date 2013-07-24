function bj_sig_chs_roc_examplefig
%Written by Xing 17/05/13
%Set useISI to 1: based on pre-test vs test, not on sample vs test.
%Set useISI to 0: sample vs test.
%Calculate ROC values based on cumulative spike data across channels, not
%just on that from individual channels.
%To compare cumulatively-calculated ROC values obtained using 2 methods:
%1. sglroc3 (blue, old) and 2. the mean trial-wise higher/lower activity (red, new).
%Also plots distributions of sample- and test-evoked activity, and ROC
%curves, for each of the two methods.
%Set exampleFig to 0 to plot ROC curves for new and old methods, set to 1
%to only plot example figures of distributions of stimulus-evoked activity
%and condition-dependent ROC curves.
excludeSessHighSSE=0;
useColMap=1;
analysisType='ROC_zero_one';
sglroc3IndividualChs=1;%set to 0 to read ROC values for individual channels and calculate mean ROC across channels; set to 1 to calculate ROCs based on pooled activity across channels
onExternalHD=0;
if onExternalHD==1
    rootFolder='G:\PL_backup_060413';
else
    rootFolder='F:';
end
calculateTangent=1;
plotSlopeFig=0;
if plotSlopeFig==1
    %example figures for channels with significant changes in slope at 30%:
    animals=[{'blanco'} {'jack'} {'jack'} {'blanco'} {'jack'}];
    animals=[{'blanco'} {'jack'} {'jack'} {'blanco'}];
    allChannels=[{[12 36 51 52 53]} {[1 3 5 6 8 10 24 39 40 41 52 53 54 56]} {[18 19 21 51]} {[2 24 40 42 49]} {[28]}];
    allChannels=[{[12 36 51 52 53]} {[5 6 8 10 24 39 40 41 52 53 54 56]} {[18 21 51]} {[2 24 40 42 49]}];%partial correlation results
    allChannels=[{[12 36 50 51 52]} {[6 8 10 24 41 52 53]} {[18 21 30 51]} {[2 7 24 49]}];%partial correlation results
    areas=[{'v4_1'} {'v4_1'} {'v1_1'} {'v4_1'} {'v1_1'}];
    areas=[{'v4_1'} {'v4_1'} {'v1_1'} {'v4_1'}];
    if excludeSessHighSSE==0
        animals=[{'blanco'} {'jack'} {'jack'} {'blanco'} {'jack'}];
        allChannels=[{[36 51 52]} {[6 8 10 24 41 52 53]} {[18 19 21 51]} {[24 49]} {[21 26]}];%Spearman's correlation results, Bonferroni corrected
        areas=[{'v4_1'} {'v4_1'} {'v1_1'} {'v4_1'} {'v1_1'}];
    end
    figSlope=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.6, 0.8]); %
    set(figSlope, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
    sampleContrast=30;
    allChInd=0;
    for animalInd=1:length(animals)
        animal=animals{animalInd};
        channels=allChannels{animalInd};
        area=areas{animalInd};
        for chInd=1:length(channels)
            allChInd=allChInd+1;
            %         figROCnew=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
            %         set(figROCnew, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
            [sampleContrasts testContrasts]=area_metadata(area);
            testContrast=testContrasts;
            matname=['ROC_Ch',num2str(channels(chInd)),'_',num2str(sampleContrast),'_1024_to_1536'];
            pathname=fullfile(rootFolder,'PL',analysisType,'ROC',animal,area,matname);
            loadText=['load ',pathname,'.mat'];
            eval(loadText)
            %         subplot(ceil(29/5),5,allChInd);
            subplot(ceil(20/5),5,allChInd);
            xvals=testContrast(1):1:testContrast(end);
%             copperCols=colormap(copper(size(ROCmat,1)));
%             copperCols=colormap(cool(size(ROCmat,1)+ceil(size(ROCmat,1)*0.25)));
            copperCols1=[];
            copperCols2=[];
            for colMapInd=1:ceil(size(ROCmat,1)/2)
                copperCols1(colMapInd,:)=[1 0 (colMapInd-1)/ceil((size(ROCmat,1))/2)];
            end
            for colMapInd=1:size(ROCmat,1)-size(ROCmat,1)/2
                copperCols2(colMapInd,:)=[1-colMapInd/floor((size(ROCmat,1))/2) 0 1];
            end
            copperCols=[copperCols1;copperCols2];
            for sessionInd=1:size(ROCmat,1)
                datavals=ROCmat{sessionInd,3};
                if sum(datavals(1:3))<=sum(datavals(end-2:end))
                    X0=[2 30 0.2 0.1];
                elseif sum(datavals(1:3))>sum(datavals(end-2:end))
                    X0=[-2 30 0.2 0.1];
                end
                options = optimset('Display','off','MaxFunEvals',10^4,'MaxIter',10^4,'TolFun',1.0E-6,'TolX',1.0E-6);
                X1=fminsearch(@fit_weibull,X0,options,testContrast,datavals,[],'least_square',[1 1 1 0],[10 100 1 0],[1 1 0 0],[-20 0 0 0]);
                %         X=fminsearch(@fit_weibull,X1,options,testContrast,datavals,[],'mle',[1 1 1 1],[2 100 1 0.2],[1 1 0 0],[-20 0 0 0]);
                X=fminsearch(@fit_weibull,X1,options,testContrast,datavals,[],'mle',[1 1 1 0],[10 100 1 1],[1 1 0 0],[-10 0 0 0]);
                fitted_yvals=1-X(4)-X(3).*exp(-(testContrast./X(2)).^X(1));
                residuals=datavals-fitted_yvals;
                sseCRF=sum(residuals.^2);
                chSSE(chInd,1:2)=[chInd sseCRF];
                if sseCRF>0.1%if the fit seems poor, try a variety of values for the upper and lower limits
                    coefEsts2 = zeros(6,4);
                    InitVar=1;
                    for upperMax=[X1(3) X1(3)+0.1]
                        for lowerMin=[X1(4) X1(4)+0.1]
                            [coefEsts2(InitVar,:)]=fminsearch(@fit_weibull,X1,options,testContrast,datavals,[],'mle',[1 1 1 0],[20 100 upperMax 0],[1 1 0 1],[-20 0 0 lowerMin]);
                            fitted_yvals=1-coefEsts2(InitVar,4)-coefEsts2(InitVar,3).*exp(-(testContrast./coefEsts2(InitVar,2)).^coefEsts2(InitVar,1));
                            residuals=datavals-fitted_yvals;
                            sseCRFtemp(InitVar)=sum(residuals.^2);
                            InitVar=InitVar+1;
                        end
                    end
                    [minSSE I]=min(sseCRFtemp);
                    if minSSE<sseCRF
                        X=coefEsts2(I,:);
                        chSSE(chInd,1:2)=[chInd minSSE];
                    end
                end
                if calculateTangent==0
                    slopeNeuro(1,chInd)=X(1);
                elseif calculateTangent==1
                    slopeNeuro(1,chInd)=X(1)*X(3)*exp(-(sampleContrast/X(2))^X(1))*sampleContrast^(X(1)-1)*(1/X(2))^X(1);
                end
                
                yvals=1-X(4)-X(3).*exp(-(xvals./X(2)).^X(1));
                xvalsFine=testContrast(1):0.01:testContrast(end);
                yvalsFine=1-X(4)-X(3).*exp(-(xvalsFine./X(2)).^X(1));
                if max(yvals)<0.5%out of range
                    c50(1,chInd)=100;
                elseif min(yvals)>0.5
                    c50(1,chInd)=0;
                else
                    diffTemp=yvalsFine-0.5;
                    [tempVal columnInd]=min(abs(diffTemp));
                    c50(1,chInd)=xvalsFine(columnInd);
                end
                hold on
                yvals=1-X(4)-X(3).*exp(-(xvals./X(2)).^X(1));
                if useColMap==1
                    plot(xvals,yvals,'Color',copperCols(sessionInd,:));
                    line(c50(1,chInd),0:0.01:1,'Color',copperCols(sessionInd,:));
                else
                    plot(xvals,yvals,'Color',[1-sessionInd/size(ROCmat,1) 0 sessionInd/size(ROCmat,1)]);
                    line(c50(1,chInd),0:0.01:1,'Color',[1-sessionInd/size(ROCmat,1) 0 sessionInd/size(ROCmat,1)]);
                end
                hold on
                if strcmp(area,'v4_1')
                    xlim([10 60]);
                elseif strcmp(area,'v1_1')
                    xlim([5 90]);
                end
                title(num2str(channels(chInd)));
                title(allChInd);
                if allChInd==1
                    xlabel('contrast (%)');
                    ylabel('AUROC');
                end
            end
        end
    end
    imagename='example_sig_chs_change_slope';
    if excludeSessHighSSE==0
        imagename=[imagename,'_allSess'];
    end
    if useColMap==1
        imagename=[imagename,'_copper'];
        imagename=[imagename,'_cool'];
    end
    pathname=fullfile(rootFolder,'PL',analysisType,imagename);
    printtext=sprintf('print -dpng %s.png',pathname);
    set(gcf,'PaperPositionMode','auto')
    eval(printtext);
end

%example figures for channels with significant changes in PNE:
animals=[{'blanco'} {'jack'} {'blanco'} {'jack'}];
allChannels=[{[24 22 36 50 51 52 60]} {[8 10 24 41 53]} {[26]} {[32 55 18 51]}];%partial correlation results
areas=[{'v4_1'} {'v4_1'} {'v1_1'} {'v1_1'}];
animals=[{'blanco'} {'jack'} {'jack'} {'blanco'} {'jack'}];
allChannels=[{[7 24 36 60]} {[8 10 24 41 53]} {[18]} {[26]} {[9 17 32 55]}];
areas=[{'v4_1'} {'v4_1'} {'v1_1'} {'v1_1'} {'v1_1'}];
if excludeSessHighSSE==0
    animals=[{'blanco'} {'jack'} {'blanco'} {'blanco'} {'jack'}];
    allChannels=[{[7 36 51 60]} {[8 10 41 53]} {24} {26} {[9 17 32 55]}];
    areas=[{'v4_1'} {'v4_1'} {'v4_1'} {'v1_1'} {'v1_1'}];
end
figPNE=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0, 0.6, 0.6]); %
set(figPNE, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
sampleContrast=30;
allChInd=0;
for animalInd=1:length(animals)
    animal=animals{animalInd};
    channels=allChannels{animalInd};
    area=areas{animalInd};
    for chInd=1:length(channels)
        allChInd=allChInd+1;
        %         figROCnew=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
        %         set(figROCnew, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
        [sampleContrasts testContrasts]=area_metadata(area);
        testContrast=testContrasts;
        matname=['ROC_Ch',num2str(channels(chInd)),'_',num2str(sampleContrast),'_1024_to_1536'];
        pathname=fullfile(rootFolder,'PL',analysisType,'ROC',animal,area,matname);
        loadText=['load ',pathname,'.mat'];
        eval(loadText)
        %         subplot(ceil(29/5),5,allChInd);
        subplot(ceil(15/5),5,allChInd);
        xvals=testContrast(1):1:testContrast(end);
%         copperCols=colormap(copper(size(ROCmat,1)));
%         copperCols=[colormap(cool(ceil(size(ROCmat,1)/2)));colormap(spring(floor(size(ROCmat,1)/2)))];
        copperCols1=[];
        copperCols2=[];
        for colMapInd=1:ceil(size(ROCmat,1)/2)
            copperCols1(colMapInd,:)=[1 0 (colMapInd-1)/ceil((size(ROCmat,1))/2)];
        end
        for colMapInd=1:size(ROCmat,1)-size(ROCmat,1)/2
            copperCols2(colMapInd,:)=[1-colMapInd/floor((size(ROCmat,1))/2) 0 1];
        end
        copperCols=[copperCols1;copperCols2];
        for sessionInd=1:size(ROCmat,1)
            datavals=ROCmat{sessionInd,3};
            if sum(datavals(1:3))<=sum(datavals(end-2:end))
                X0=[2 30 0.2 0.1];
            elseif sum(datavals(1:3))>sum(datavals(end-2:end))
                X0=[-2 30 0.2 0.1];
            end
            options = optimset('Display','off','MaxFunEvals',10^4,'MaxIter',10^4,'TolFun',1.0E-6,'TolX',1.0E-6);
            X1=fminsearch(@fit_weibull,X0,options,testContrast,datavals,[],'least_square',[1 1 1 0],[10 100 1 0],[1 1 0 0],[-20 0 0 0]);
            %         X=fminsearch(@fit_weibull,X1,options,testContrast,datavals,[],'mle',[1 1 1 1],[2 100 1 0.2],[1 1 0 0],[-20 0 0 0]);
            X=fminsearch(@fit_weibull,X1,options,testContrast,datavals,[],'mle',[1 1 1 0],[10 100 1 1],[1 1 0 0],[-10 0 0 0]);
            fitted_yvals=1-X(4)-X(3).*exp(-(testContrast./X(2)).^X(1));
            residuals=datavals-fitted_yvals;
            sseCRF=sum(residuals.^2);
            chSSE(chInd,1:2)=[chInd sseCRF];
            if sseCRF>0.1%if the fit seems poor, try a variety of values for the upper and lower limits
                coefEsts2 = zeros(6,4);
                InitVar=1;
                for upperMax=[X1(3) X1(3)+0.1]
                    for lowerMin=[X1(4) X1(4)+0.1]
                        [coefEsts2(InitVar,:)]=fminsearch(@fit_weibull,X1,options,testContrast,datavals,[],'mle',[1 1 1 0],[20 100 upperMax 0],[1 1 0 1],[-20 0 0 lowerMin]);
                        fitted_yvals=1-coefEsts2(InitVar,4)-coefEsts2(InitVar,3).*exp(-(testContrast./coefEsts2(InitVar,2)).^coefEsts2(InitVar,1));
                        residuals=datavals-fitted_yvals;
                        sseCRFtemp(InitVar)=sum(residuals.^2);
                        InitVar=InitVar+1;
                    end
                end
                [minSSE I]=min(sseCRFtemp);
                if minSSE<sseCRF
                    X=coefEsts2(I,:);
                    chSSE(chInd,1:2)=[chInd minSSE];
                end
            end
            if calculateTangent==0
                slopeNeuro(1,chInd)=X(1);
            elseif calculateTangent==1
                slopeNeuro(1,chInd)=X(1)*X(3)*exp(-(sampleContrast/X(2))^X(1))*sampleContrast^(X(1)-1)*(1/X(2))^X(1);
            end
            
            yvals=1-X(4)-X(3).*exp(-(xvals./X(2)).^X(1));
            xvalsFine=testContrast(1):0.01:testContrast(end);
            yvalsFine=1-X(4)-X(3).*exp(-(xvalsFine./X(2)).^X(1));
            if max(yvals)<0.5%out of range
                c50(1,chInd)=100;
            elseif min(yvals)>0.5
                c50(1,chInd)=0;
            else
                diffTemp=yvalsFine-0.5;
                [tempVal columnInd]=min(abs(diffTemp));
                c50(1,chInd)=xvalsFine(columnInd);
            end
            hold on
            if useColMap==1
                plot([c50(1,chInd) c50(1,chInd)],[0 1],'Color',copperCols(sessionInd,:));
            else
                plot([c50(1,chInd) c50(1,chInd)],[0 1],'Color',[1-sessionInd/size(ROCmat,1) 0 sessionInd/size(ROCmat,1)]);
            end
            yvals=1-X(4)-X(3).*exp(-(xvals./X(2)).^X(1));
            if useColMap==1
                plot(xvals,yvals,'Color',copperCols(sessionInd,:));
                line(c50(1,chInd),0:0.01:1,'Color',copperCols(sessionInd,:));
            else
                plot(xvals,yvals,'Color',[1-sessionInd/size(ROCmat,1) 0 sessionInd/size(ROCmat,1)]);
                line(c50(1,chInd),0:0.01:1,'Color',[1-sessionInd/size(ROCmat,1) 0 sessionInd/size(ROCmat,1)]);
            end
            hold on
            if strcmp(area,'v4_1')
                xlim([10 60]);
            elseif strcmp(area,'v1_1')
                xlim([5 90]);
            end
            title(num2str(channels(chInd)));
            title(allChInd);
            if allChInd==1
                xlabel('contrast (%)');
                ylabel('AUROC');
            end
        end
    end
end
imagename='example_sig_chs_change_PNE';
if excludeSessHighSSE==0
    imagename=[imagename,'_allSess'];
end
if useColMap==1
    imagename=[imagename,'_copper'];
    imagename=[imagename,'_cool'];
end
pathname=fullfile(rootFolder,'PL',analysisType,imagename);
printtext=sprintf('print -dpng %s.png',pathname);
set(gcf,'PaperPositionMode','auto')
eval(printtext);
