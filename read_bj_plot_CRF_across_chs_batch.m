function read_bj_plot_CRF_across_chs_batch
%Written on 07/03/13.

test_epochs={0 529 529*2 529*3};durSpon=150;
allChCoefs=cell(1,6);
saveImageAppendText={'_slope_vs_session' '_C50_vs_session' '_neuro_vs_psycho' '_diff_C50_vs_session' '_minRate_vs_session' '_maxRate_vs_session'};
areas=[{'v4_1'} {'v1_1'}];
animals=[{'blanco'} {'jack'}];
plotFig=1;
figureHandles=[{'fighandle1'} {'fighandle2'} {'fighandle3'} {'fighandle4'} {'fighandle5'} {'fighandle6'}];
totalSigChs=0;
for animalInd=1:2
    for areaInd=1:2
        animal=animals{animalInd};
        area=areas{areaInd};
        if strcmp(animal,'jack')
            if strcmp(area,'v4_1')
                testContrasts=[10 15 20 25 27 28 29 31 32 33 35 40 50 60];
                sampleContrasts=30;
            elseif strcmp(area,'v1_1')
                testContrasts=[5 10 15 20 22 25 28 32 35 40 45 50 60 90];
                sampleContrasts=30;
            elseif strcmp(area,'v1_2')
                testContrasts=[5 10 12 15 18 22 25 28 35 45 63 90;5 10 15 22 25 28 32 35 38 45 60 90;5 10 15 25 32 35 38 42 45 50 60 90];
                sampleContrasts=[20 30 40];
            end
        elseif strcmp(animal,'blanco')
            if strcmp(area,'v4_1')
                testContrasts=[10 15 20 25 27 28 29 31 32 33 35 40 50 60];
                sampleContrasts=30;
            elseif strcmp(area,'v1_1')
                testContrasts=[5 10 15 20 22 25 28 32 35 40 45 50 60 90];
                sampleContrasts=30;
            elseif strcmp(area,'v1_2')
                testContrasts=[5 10 12 15 18 22 25 28 35 45 63 90;5 10 15 22 25 28 32 35 38 45 60 90;5 10 15 25 32 35 38 42 45 50 60 90];
                sampleContrasts=[20 30 40];
            end
        end  
        channels = main_channels(animal,area);
        for epoch=1:size(test_epochs,2)
            if epoch==1
                periods=[-durSpon 0];
            else
                periods=[test_epochs{epoch-1} test_epochs{epoch}(1)];
            end
            for subPeriod=1:length(periods)-1
                startEndTime=['_',num2str(periods(subPeriod)),'_to_',num2str(periods(subPeriod+1))];
                for j=1:length(sampleContrasts)
                    sampleContrast=sampleContrasts(j);
                    testContrast=testContrasts(j,:);      
                    if plotFig==1
                        fighandle1=  figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
                        set(fighandle1, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                        fighandle2=  figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
                        set(fighandle2, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                        fighandle3=  figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
                        set(fighandle3, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                        fighandle4=  figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
                        set(fighandle4, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                        fighandle5=  figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
                        set(fighandle5, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                        fighandle6=  figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
                        set(fighandle6, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                    end
                    for i=1:length(channels)
                        chNum=channels(i);
                        coefsFileName=[num2str(chNum),'_',num2str(sampleContrast),startEndTime,'_crf_coefs_',area];%,'_goodSSE_no_outliers_sl3_C503','.mat'];
                        coefsPathName=fullfile('F:','PL','CRF',animal,'CRF_coef_mat',coefsFileName);
                        loadText=['load ',coefsPathName];
                        eval(loadText);
                        if find(coefficients(2,:)<0.05)
                            totalSigChs=totalSigChs+1;
                        end
                        xVals=[{1:length(slopeNeuro)} {1:length(c50)} {matchPsycho} {1:length(diffc50)}  {1:length(minRate)}  {1:length(maxRate)}];
                        yVals=[{slopeNeuro} {c50} {slopeNeuro} {diffc50} {minRate} {maxRate}];
                        for k=1:size(coefficients,2)
                            if plotFig==1
                                figText=sprintf('figure(%s)',figureHandles{k});
                                eval(figText)
                                subplot(ceil(length(channels)/5),5,i);
                                plot(xVals{k},yVals{k},'ok');
                                yLimVals=get(gca,'YLim');
                                xLimVals=get(gca,'XLim');
                                if k==1||k==4
                                    if xLimVals(1)>=0
                                        xlim([0 xLimVals(2)]);
                                    end
                                    if yLimVals(1)>=0
                                        ylim([0 yLimVals(2)]);
                                    end
                                end
                                ptext=sprintf('R: %.3f p: %.4f',coefficients(1,k),coefficients(2,k));
                                if coefficients(2,k)<0.05
                                    colText='g';
                                else
                                    colText='k';
                                end
                                yLimVals=get(gca,'YLim');
                                xLimVals=get(gca,'XLim');
                                title(num2str(chNum));
                                text('Position',[xLimVals(1)+(xLimVals(2)-xLimVals(1))/5 yLimVals(1)+(yLimVals(2)-yLimVals(1))/6],'FontSize',9,'String',ptext,'Color',colText);
                            end
                            allChCoefs{k}=[allChCoefs{k};{animal},{area},chNum,coefficients(1,k),coefficients(2,k)];
                        end
                    end
                end
            end
            saveMatName=['all_',animal,'_',area,'_CRF_coefs'];
            saveMatFolder=fullfile('F:','PL','CRF',animal,'CRF_coef_mat');
            saveMatPath=fullfile(saveMatFolder,saveMatName);
            saveText=['save ',saveMatPath,'.mat allChCoefs'];
            eval(saveText)
            for k=1:size(coefficients,2)
                figText=sprintf('figure(%s)',figureHandles{k});
                eval(figText)
                saveCoefImageName=['all_',animal,'_',area,'_CRF_coefs',startEndTime,saveImageAppendText{k}];
                saveCoefImageFolder=fullfile('F:','PL','CRF',animal,'CRF_coef_images');
                saveCoefImagePath=fullfile(saveCoefImageFolder,saveCoefImageName);
                printtext=sprintf('print -dpng %s.png',saveCoefImagePath);
                set(gcf, 'PaperOrientation', 'landscape');
                set(gcf, 'PaperPositionMode', 'auto');
                eval(printtext);
            end
            close all
        end
    end
end

