function bj_calc_roc_diff

%Written by Xing 11/04/13
%to calculate difference between ROC values obtained using 2 methods:
%1. sglroc3 and 2. the mean trial-wise higher/lower activity.
analysisType='ROC';
onExternalHD=0;
if onExternalHD==1
    rootFolder='G:\PL_backup_060413';
else
    rootFolder='F:';
end
animals=[{'blanco'} {'jack'}];
areas=[{'v4_1'} {'v4_2'} {'v1_1'} {'v1_2'}];
areas=[{'v4_1'} {'v1_1'} {'v1_2'}];
test_epochs={0 512 512*2 512*3};durSpon=150;
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        [sampleContrasts testContrasts]=area_metadata(area);
        channels=main_channels(animal,area);
        psychoname=['psycho_constants_',area];
        psychoPathname=fullfile('F:','PL','psycho_data',animal,psychoname);
        figSess=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
        set(figSess, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
        figCond=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
        set(figCond, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
        for i=1:length(channels)
            for sampleContrastsInd=1:length(sampleContrasts)
                sampleContrast=sampleContrasts(sampleContrastsInd);
                testContrast=testContrasts(sampleContrastsInd,:);
                colmapText=colormap(jet(size(testContrast,2)));
                colmapText=[colmapText(1,:);132/255 22/255 216/255;202/255 65/255 223/255;colmapText(3:7,:);157/255 212/255 61/255;colmapText(10:12,:);178/255 111/255 12/255;colmapText(end,:)];
                for epoch=1:size(test_epochs,2)
                    if strcmp(analysisType,'CRF')||strcmp(analysisType,'ROC')&&epoch==4||strcmp(analysisType,'NVP')&&epoch==4
                        if epoch==1
                            periods=[-durSpon 0];
                        else
                            periods=[test_epochs{epoch-1} test_epochs{epoch}(1)];
                        end
                        for subPeriod=1:length(periods)-1
                            startEndTime=['_',num2str(periods(subPeriod)),'_to_',num2str(periods(subPeriod+1))];
                            if strcmp(analysisType,'CRF')||strcmp(analysisType,'ROC')
                                matName=[analysisType,'_Ch',num2str(channels(i)),'_',num2str(sampleContrast),startEndTime,'.mat'];
                                matPath=fullfile('F:','PL',analysisType,animal,area,matName);
                                loadText=['load ',matPath,' ',analysisType,'mat'];
                                eval(loadText);
                                dataArrayNew=ROCmat;  
                                if strcmp(area,'v4_1')||strcmp(area,'v1_1')
                                    matName=[analysisType,'_Ch',num2str(channels(i)),'_',num2str(sampleContrast),'_1058_to_1587','.mat'];
                                end
                                matPath=fullfile('F:','PL','ROC_sglroc3',analysisType,animal,area,matName);
                                loadText=['load ',matPath,' ',analysisType,'mat'];
                                eval(loadText);
                                dataArrayOld=ROCmat;
                            elseif strcmp(analysisType,'NVP')
                                matName=['ROC_Ch',num2str(channels(i)),'_',num2str(sampleContrast),startEndTime,'.mat'];
                                matPath=fullfile('F:','PL','ROC',animal,area,matName);
                                loadText=['load ',matPath,' ROCmat'];
                                eval(loadText);
                            end
                            if strcmp(analysisType,'CRF')
                                dataArray=CRFmat;
                            elseif strcmp(analysisType,'NVP')
                                dataArray=ROCmat;
                            end
                            figure(figSess);
                            subplot(ceil(length(channels)/5),5,i);
                            diff=[];
                            for sessInd=1:size(dataArrayNew,1)
                                rowInd=[];
                                for rowIndOld=1:size(dataArrayOld,1)
                                    if dataArrayNew{sessInd,1}==dataArrayOld{rowIndOld,1}
                                        rowInd=rowIndOld;
                                    end
                                end
                                if ~isempty(rowInd)
                                    diff(sessInd,:)=dataArrayNew{sessInd,3}-dataArrayOld{rowInd,3};
                                    for cond=1:length(testContrast)
                                        plot(cond,diff(sessInd,cond),'x','Color',colmapText(cond,:));hold on
                                    end
                                end
                            end
                            ylim([-0.3 0.3]);
                            figure(figCond);
                            subplot(ceil(length(channels)/5),5,i);
                            meanDiff=mean(diff,1);
                            for cond=1:length(testContrast)
                                plot(cond,meanDiff(cond),'x','Color',colmapText(cond,:));hold on
                            end
                            plot([0 length(testContrast)+1],[0 0],'k:');
                            xlim([0 length(testContrast)+1]);
                            ylim([-0.1 0.1]);
                        end
                    end
                end
            end
        end
        figure(figSess);
        imagename=['diff_ROCs_old_new_',area];
        pathname=fullfile(rootFolder,'PL',analysisType,animal,imagename);
        printtext=sprintf('print -dpng %s.png',pathname);
        set(gcf,'PaperPositionMode','auto')
        eval(printtext);
        figure(figCond);
        imagename=['diff_ROCs_old_new_session_means_',area];
        pathname=fullfile(rootFolder,'PL',analysisType,animal,imagename);
        printtext=sprintf('print -dpng %s.png',pathname);
        set(gcf,'PaperPositionMode','auto')
        eval(printtext);
    end
end
