function read_bj_plot_CRF_combine_chs_batch
%Written on 07/03/13. Calculates mean CRF values across channels.

test_epochs={0 529 529*2 529*3};durSpon=150;
epochNames={'spontaneous' 'sample' 'ISI' 'test'};
areas=[{'v4_1'} {'v1_1'}];
animals=[{'blanco'} {'jack'}];
plotFig=1;
colTexts='ymcrgbk';
markerTexts='+x';
markerSizes=[5 6];
for animalInd=1:2
    for areaInd=1:2
        animal=animals{animalInd};
        area=areas{areaInd};
        sessions=main_raw_sessions_final(animal,area);
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
        for j=1:length(sampleContrasts)
            sampleContrast=sampleContrasts(j);
            testContrast=testContrasts(j,:);
            figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
            for epoch=1:size(test_epochs,2)
                epochTitle=epochNames{epoch};
                if epoch==1
                    periods=[-durSpon 0];
                else
                    periods=[test_epochs{epoch-1} test_epochs{epoch}(1)];
                end
                for subPeriod=1:length(periods)-1
                    startEndTime=['_',num2str(periods(subPeriod)),'_to_',num2str(periods(subPeriod+1))];
                    periodTitle=[epochTitle,' (',num2str(periods(subPeriod)),' to ',num2str(periods(subPeriod+1)),' ms)'];
                    clear allCRFvals
                    [allCRFvals{1,1:length(sessions)}] = deal(zeros(1));
                    for i=1:length(channels)
                        startEndTime=['_',num2str(periods(subPeriod)),'_to_',num2str(periods(subPeriod+1))];
                        CRFmatName=['CRF_Ch',num2str(channels(i)),'_',num2str(sampleContrast),startEndTime,'.mat'];
                        CRFmatPath=fullfile('F:','PL','CRF',animal,area,CRFmatName);
                        loadText=['load ',CRFmatPath,' CRFmat'];
                        eval(loadText);
                        for rowInd=1:size(CRFmat,1)
                            allCRFvals{rowInd}=allCRFvals{rowInd}+CRFmat{rowInd,3};
                        end
                    end
                    for rowInd=1:size(CRFmat,1)
                        allCRFvals{rowInd}=allCRFvals{rowInd}/length(channels);
                    end
                    saveMatName=['mean_',animal,'_',area,'_CRF_across_chs',startEndTime];
                    saveMatFolder=fullfile('F:','PL','CRF',animal);
                    saveMatPath=fullfile(saveMatFolder,saveMatName);
                    saveText=['save ',saveMatPath,'.mat allCRFvals'];
                    eval(saveText)
                    %proportion of trials monkey reported that test was higher contrast,
                    colmapText=colormap(jet(size(testContrast,2)));
                    colmapText=[colmapText(1,:);132/255 22/255 216/255;202/255 65/255 223/255;colmapText(3:7,:);157/255 212/255 61/255;colmapText(10:12,:);178/255 111/255 12/255;colmapText(end,:)];
                    subplot(2,2,epoch);
                    markerText=markerTexts(2);markerS=8;
                    for i=1:length(testContrast)
                        CRFperCond=[];
                        for sessionInd=1:size(allCRFvals,2)
                            CRFperCond=[CRFperCond allCRFvals{sessionInd}(i)];
                        end
                        plot(1:size(allCRFvals,2),CRFperCond,'Color',colmapText(i,:),'LineStyle','none','Marker',markerText,'MarkerFaceColor',colmapText(i,:),'MarkerEdgeColor',colmapText(i,:),'MarkerSize',markerS);hold on%'MarkerFaceColor',[1/i 1/i 1/i],
                        hold on
                        %                         [chiseb1linear(i,:) chiseb1(i,:) coefperfb1linear(i,:) coefperfb1(i,:) aRSb1linear(i,:) aRSb1(i,:)]=linearexpo_fitting(1:size(allCRFvals,2),CRFperCond',i,1);
                        %     yLimVals=get(gca,'ylim');
                        %     xLimVals=get(gca,'xlim');
                        %     unitSpace=(yLimVals(2)-yLimVals(1))/30;
                    end
                    title(periodTitle);
                    xlabel('session');
                    ylabel('mean firing rate (spikes/s)');
                    saveCoefImageName=[area,'_',animal,'_4epochs_mean_CRF_across_channels_vs_session_',num2str(sampleContrast)];
                    saveCoefImageFolder=fullfile('F:','PL','CRF',animal);
                    saveCoefImagePath=fullfile(saveCoefImageFolder,saveCoefImageName);
                    printtext=sprintf('print -dpng %s.png',saveCoefImagePath);
                    eval(printtext);
                    %                 saveCoefImagePath
                    %                     ylim([0 100]);
                    %                     legend('hide');
                    %                     xlabel('');
                    %                     ylabel('');
                    %                     pc_condcoefFig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
                    %                     subplot(2,2,1);
                    %                     [rperfcoeff(1:2) pperfcpef(1:2)]=plot_coef_expo_perf(coefperfb1,sampleContrast,testContrast,ind);
                    %
                    %                     %to check values of coefficients:
                    %                     figure;plot(1:14,coefperfb1(:,1));
                    %                     figure;plot(1:14,coefperfb1(:,2));
                    %                     figure;plot(1:14,coefperfb1(:,3));
                    %
                    %                     %to create test data for exp curve:
%                     xtest=1:14;
%                     figure
%                     for testind=1:14
%                         ytest(testind)=exp(-testind);
%                     end
%                     plot(xtest,ytest)
%                     figure
%                     for testind=1:14
%                         ytest(testind)=-exp(-testind);
%                     end
%                     plot(xtest,ytest)
                end
            end
        end
    end
    close all
end

