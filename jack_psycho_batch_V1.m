function jack_psycho_batch_V1
%Similar to blanco_SE_roc_batch, contains list of sessions and data
%directories, however, as this analyses psychometric performance, does not
%take spike activity into account and thus the function
%'blanco_2target_psycho_EV_v2' is only run once per session.
%(blanco_SE_roc_batch runs the function blanco_SE_roc 64 times per session- once for each channel).
%If intention is to create data file containing list of Weibull constants
%from each session, must first activate lines of code in
%blanco_2target_psycho_EV_v2 so that values are written to file.

% allMeanPerf=[];save allMeanPerf_20.mat allMeanPerf
% allMeanPerf=[];save allMeanPerf_30.mat allMeanPerf
% allMeanPerf=[];save allMeanPerf_40.mat allMeanPerf

%set accordingly:
% monkeys=[{'blanco'} {'jack'}];
% areas=[{'v4_1'} {'v4_2'} {'v1'}];
roving=1;
readData=0;
% for monkeyCount=1:2
%     for areaCount=1:3
%         monkey=monkeys{monkeyCount};
%         area=areas{areaCount};
        monkey='blanco';area='v1_2';animal=monkey;
        %no need to adjust:
        sampleContrast=30;
        appendText=[];
        splitSess=0;
        if strcmp(monkey,'jack')
            appendText='j_';
            if strcmp(area,'v4_1')
                sessionNums=[1:5];
                sessionNums=[24 26 28:49];
                sessionNums=24:50;
            elseif strcmp(area,'v4_2')
                sessionNums=73:77;
            elseif strcmp(area,'v1')
                sessionNums=51:72;
            elseif strcmp(area,'v1_2')
                %         sessionNums=[78:93];
                %         sessionNums=94:115;
                %         sessionNums=116:119;
                sessionNums=78:119;
                sessionNums=116:120;
            elseif strcmp(area,'v1_3')
                sessionNums=128:129;%126:127; v1_3 is bu/pl control task. 
            elseif strcmp(area,'v1_4')%v1_4 is jack's control task with stim in blanco's v1 location
                sessionNums=183:183;%139 is non-roving, 14 conditions. 140 onwards is roving, 36 conditions. 163 onwards is with flankers.
            end
        elseif strcmp(monkey,'blanco')
            splitSessions=[355 405 435];
            sessionNums=main_raw_sessions_final(animal,area);
        end
        analysisFolderAppend=[];
        for i=1:length(sessionNums)
            if strcmp(monkey,'jack')
                [file_of_int,testContrast,sampleContrasts,roving]=session_metadata(sessionNums(i),monkey);
                if sessionNums(i)==22
                    file_of_int='2161398.1';
                    testContrast=[10 15 20 25 35 40 60 90];
                elseif sessionNums(i)<22
                    file_of_int='2161392.1';
                    testContrast=[10 90];
                elseif sessionNums(i)==23
                    file_of_int='21613912.1';%wrongly named in Cheetah encode- actual name should be 21613612.1
                    testContrast=[10 15 20 25 27 29 31 33 35 40 50 60];
                elseif sessionNums(i)>23&&sessionNums(i)<51
                    file_of_int='21613614.1';
                    testContrast=[10 15 20 25 27 28 29 31 32 33 35 40 50 60];
                elseif sessionNums(i)>50&&sessionNums(i)<73
                    file_of_int='47553914.1';
                    testContrast=[5 10 15 20 22 25 28 32 35 40 45 50 60 90];
                elseif sessionNums(i)>72&&sessionNums(i)<78
                    file_of_int='21613614.1';
                    testContrast=[10 15 20 25 27 28 29 31 32 33 35 40 50 60];
                elseif sessionNums(i)>77&&sessionNums(i)<94
                    if sessionNums(i)==89
                        file_of_int='47523412.2';
                    else
                        file_of_int='47523412.1';
                    end
                    roving=1;
                    analysisFolderAppend='_4';
                elseif sessionNums(i)>93&&sessionNums(i)<116
                    file_of_int='4723412f.1';
                    roving=1;
                    analysisFolderAppend='_5';
                elseif sessionNums(i)>115&&sessionNums(i)<121
                    file_of_int='47523412.1';
                    roving=1;
                    analysisFolderAppend='_6';
                elseif sessionNums(i)>120&&sessionNums(i)<126
                    file_of_int='bu_jack.1';
%                     roving=1;
%                     analysisFolderAppend='_6';
                elseif sessionNums(i)==126
                    file_of_int='47553914.1';%PL task with just PL stimuli in RF, no BU stimuli present
                    analysisFolderAppend='_c';
%                     file_of_int='4753914bp.1';%control BU task with PL stimuli in RF, and BU stimuli presented in upper hemisphere
%                     analysisFolderAppend='_b';
                    roving=0;
                    testContrast=[10 15 20 25 27 28 29 31 32 33 35 40 50 60];
                elseif sessionNums(i)>126&&sessionNums(i)<129
                    file_of_int='4753914c.1';%control PL task with PL stimuli in RF, and BU stimuli presented in upper hemisphere
                    analysisFolderAppend='_c';
%                     file_of_int='4753914bp.1';%control BU task with PL stimuli in RF, and BU stimuli presented in upper hemisphere
%                     analysisFolderAppend='_b';
                    roving=0;
                    testContrast=[10 15 20 25 27 28 29 31 32 33 35 40 50 60];
                elseif sessionNums(i)>139&&sessionNums(i)<163%139 has single session of non-roving 30% sample and 14 conditions
                    file_of_int='225234912.1';%control task with jack's PL stimuli in blanco's V1 RF
                    analysisFolderAppend='_7';
                elseif sessionNums(i)>162&&sessionNums(i)<184%roving samples and 36 conditions
                    file_of_int='2223412f.1';%control task with jack's PL stimuli in blanco's V1 RF, plus flankers
                    analysisFolderAppend='_7';
                end
            elseif strcmp(monkey,'blanco')
                [file_of_int,testContrast,sampleContrasts,roving]=session_metadata(sessionNums(i),monkey);
                if sessionNums(i)>305&&sessionNums(i)<343
                    file_of_int='21613614.1';
                    testContrast=[10 15 20 25 27 28 29 31 32 33 35 40 50 60];
                elseif sessionNums(i)>342&&sessionNums(i)<360
                    file_of_int='2353914.1';
                    testContrast=[5 10 15 20 22 25 28 32 35 40 45 50 60 90];
                    if find(sessionNums(i)==splitSessions)
                        splitSess=1;
                    end
                elseif sessionNums(i)>359&&sessionNums(i)<365
                    file_of_int='21613614.1';
                    testContrast=[10 15 20 25 27 28 29 31 32 33 35 40 50 60];
                elseif sessionNums(i)>387&&sessionNums(i)<423
                    if sessionNums(i)==415
                        %             file_of_int='225234912.2';
                    else
                        file_of_int='225234912.1';
                        %                 file_of_int='2223412f.1';
                    end
                    if find(sessionNums(i)==splitSessions)
                        splitSess=1;
                    end
                    roving=1;
                    analysisFolderAppend='_4';
                elseif sessionNums(i)>421&&sessionNums(i)<453
                    file_of_int='2223412f.1';
                    roving=1;
                    analysisFolderAppend='_5';
                    if find(sessionNums(i)==splitSessions)
                        splitSess=1;
                    end
                elseif sessionNums(i)>452&&sessionNums(i)<460
                    file_of_int='225234912.1';
                    roving=1;
                    analysisFolderAppend='_6';
                end
            end
            if readData==1
                if roving==0
                    if splitSess==1
                        cdText=['cd ','F:\blanco\',area,'_events_files\',num2str(sessionNums(i)),'_1'];
                    end
                    eval(cdText);
                    conditions=1:length(testContrast);
                    [vals]=jack_2target_psycho_EV_v3_mex(file_of_int,sampleContrast,testContrast,sessionNums(i),conditions,area,analysisFolderAppend,monkey);
                    if splitSess==1
                        cdText=['cd ','F:\blanco\',area,'_events_files\',num2str(sessionNums(i)),'_2'];
                        eval(cdText);
                        [vals2]=jack_2target_psycho_EV_v3_mex(file_of_int,sampleContrast,testContrast,sessionNums(i),conditions,area,analysisFolderAppend,monkey);
                        vals=[vals;vals2];
                    end
%                     analyseValsjexample(vals,testContrast,conditions,sampleContrast,monkey,area,analysisFolderAppend,sessionNums(i));
                    analyseVals2(vals,testContrast,conditions,sampleContrast,monkey,area,analysisFolderAppend,sessionNums(i),roving);
                    splitSess=0;
                elseif roving==1
                    testContrasts=[5 10 12 15 18 22 25 28 35 45 63 90;5 10 15 22 25 28 32 35 38 45 60 90;5 10 15 25 32 35 38 42 45 50 60 90];
                    sampleContrasts=[20 30 40];
                    allConditions=[13:1:24;1:1:12;25:1:36];
                    for sampleNum=1:3
                        sampleContrast=sampleContrasts(sampleNum);
                        testContrast=testContrasts(sampleNum,:);
                        conditions=allConditions(sampleNum,:);
                        if strcmp(monkey,'jack')
                            cdText=['cd ','V:\thielelab\Groups\ThieleGroup\monkey_data\Jack\_jackgrid\j_',area,'_events_files\',num2str(sessionNums(i))];
                            eval(cdText);
                        elseif strcmp(monkey,'blanco')
                            if splitSess==1
                                cdText=['cd ','F:\blanco\',area,'_events_files\',num2str(sessionNums(i)),'_1'];
                            else
                                cdText=['cd ','F:\blanco\',area,'_events_files\',num2str(sessionNums(i))];
                            end
                            eval(cdText);
                        end
                        %             [vals]=blanco_2target_psycho_EV_v4_mex(file_of_int,sampleContrast,testContrast,sessionNums(i),folder);
                        [vals]=jack_2target_psycho_EV_v3_mex(file_of_int,sampleContrast,testContrast,sessionNums(i),conditions,area,analysisFolderAppend,monkey,roving);
                        if splitSess==1
                            cdText=['cd ','F:\blanco\',area,'_events_files\',num2str(sessionNums(i)),'_2'];
                            eval(cdText);
                            [vals2]=jack_2target_psycho_EV_v3_mex(file_of_int,sampleContrast,testContrast,sessionNums(i),conditions,area,analysisFolderAppend,monkey);
                            vals=[vals;vals2];
                        end
                        savevals{sampleNum}=analyseVals2(vals,testContrast,conditions,sampleContrast,monkey,area,analysisFolderAppend,sessionNums(i),roving);
                        splitSess=0;
                    end
                    valsFileName=['vals_',num2str(sessionNums(i)),'.mat'];
                    valsFolder=fullfile('F:','PL','vals_perf',animal,valsFileName); %#ok<NASGU>
                    saveText=['save ',valsFolder,' vals'];
                    vals=savevals;
                    eval(saveText);
                end
            end
        end
%     end
% end

taskCols=[230/256 230/256 250/256;255/256 222/256 173/256;154/256 205/256 50/256];%purple;orange;green
taskColsDark=[153/256 50/256 204/256;210/256 105/256 30/256;34/256 139/256 34/256];%purple;orange;green

plotRoving=1;
if plotRoving==1
    if strcmp(animal,'jack')
        areas=[{'v1_2'} {'v1_4'}];
        areaText=[{'V1 RFs location (-0.7,-1.3)'} {'control location (-3.5,-3)'}];
    elseif strcmp(animal,'blanco')
        areas={'v1_2'};
        areaText={'V1 RFs location (-3.5,-3)'};
    end
    fig=figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.2, 0.04, 0.75/(3-length(areas)), 0.9]);
    set(fig,'PaperUnits','centimeters','PaperType','A4','PaperOrientation', 'portrait', 'PaperPosition', [0.63452 0.63452 31.5/(3-length(areas)) 28.41]);
    for areaInd=1:length(areas)
        area=areas{areaInd};
        [flankerType flankerSessions]=getFlankerType(animal,area);
        sampleContrasts=20:10:40;
        count=0;
        for flankerInd=1:length(flankerType)
            taskCol=taskCols(flankerInd,:);
            taskColDark=taskColsDark(flankerInd,:);
            markerFaceTypes=[{'none'} {taskCol} {taskColDark}];
            for i=1:length(sampleContrasts)
                perf=[];
                PSE=[];
                slope=[];
                sampleContrast=sampleContrasts(i);
                if sampleContrast==30
                    markerFaceType=markerFaceTypes{2};
                elseif sampleContrast==20
                    markerFaceType=markerFaceTypes{1};
                elseif sampleContrast==40
                    markerFaceType=markerFaceTypes{3};
                end
                loadText=['load F:\PL\psycho_data\',animal,'\allMeanPerf_',area,'_',num2str(sampleContrast),'.mat allMeanPerf'];
                eval(loadText)
                for sessionInd=1:length(flankerSessions{flankerInd})
                    rowInd=find(allMeanPerf==flankerSessions{flankerInd}(sessionInd));
                    perf=[perf allMeanPerf(rowInd,2)];
                    PSE=[PSE allMeanPerf(rowInd,15)];
                    slope=[slope allMeanPerf(rowInd,16)];
                end
                if length(areas)==1
                    subplot(3,length(areas),1);
                else
                    subplot(3,length(areas),areaInd);
                end
                plot(count+1:count+length(flankerSessions{flankerInd}),perf,'LineStyle','none','Marker','o','MarkerEdgeColor',taskColDark,'MarkerFaceColor',markerFaceType);hold on
                if length(areas)==2
                    ylim([0.5 0.95]);
                else
                    ylim([0.6 0.95]);
                end
                title(areaText{areaInd});
                if areaInd==1
                    ylabel('proportion correct');
                end
                if length(areas)==1
                    subplot(3,length(areas),2);
                else
                    subplot(3,length(areas),areaInd+2);
                end
                plot(count+1:count+length(flankerSessions{flankerInd}),PSE,'LineStyle','none','Marker','o','MarkerEdgeColor',taskColDark,'MarkerFaceColor',markerFaceType);hold on
                if length(areas)==2
                    ylim([20 50]);
                else
                    ylim([15 50]);
                end
                if areaInd==1
                    ylabel('PSE');
                end
                if length(areas)==1
                    subplot(3,length(areas),3);
                else
                    subplot(3,length(areas),areaInd+4);
                end
                plot(count+1:count+length(flankerSessions{flankerInd}),slope,'LineStyle','none','Marker','o','MarkerEdgeColor',taskColDark,'MarkerFaceColor',markerFaceType);hold on
                if length(areas)==2
                    ylim([1 9]);
                else
                    ylim([1 11]);
                end
                if areaInd==1
                    ylabel('slope');
                end
            end
            count=count+length(flankerSessions{flankerInd});
        end
    end
    subplot(3,length(areas),length(areas));
    ylims=get(gca,'YLim');
    xlims=get(gca,'XLim');
    sampleTexts=[{'20% sample'} {'30% sample'} {'40% sample'}];
    flankerTexts=[{'no flankers'} {'flankers'} {'no flankers'}];
    markerFaceTypes=[{'none'} {[180/256 180/256 180/256]} {'k'}];
    taskColsMedium=[147/256 112/256 219/256;255/256 165/256 0/256;154/256 205/256 50/256];%purple;orange;green
    for i=1:3
        if length(areas)==2
            keyText=[sampleTexts{i},'         ',flankerTexts{i}];
        else
            keyText=[sampleTexts{i},'            ',flankerTexts{i}];
        end
        plot(xlims(1)+(xlims(2)-xlims(1))/20,ylims(1)+(ylims(2)-ylims(1))/10*(4-i),'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor',markerFaceTypes{i});
        plot(xlims(1)+(xlims(2)-xlims(1))/3.6,ylims(1)+(ylims(2)-ylims(1))/10*(4-i),'LineStyle','none','Marker','o','MarkerEdgeColor',taskColsMedium(i,:),'MarkerFaceColor',taskColsMedium(i,:));
        text('Position',[xlims(1)+(xlims(2)-xlims(1))/15,ylims(1)+(ylims(2)-ylims(1))/10*(4-i)],'LineStyle','none','FontSize',9,'String',keyText,'Color','k');
    end
    formats=[{'eps'} {'png'}];
    formatsCommand=[{'epsc'} {'png'}];
    for j=1:2
        if length(areas)==2
            printFileName=['F:\PL\psycho_data\',animal,'\','PC_PSE_slope_',areas{1},'_',areas{2},'.',formats{j}];
        elseif length(areas)==1
            printFileName=['F:\PL\psycho_data\',animal,'\','PC_PSE_slope_',areas{1},'.',formats{j}];
        end
        print(sprintf('-d%s',formatsCommand{j}),'-r300',printFileName)
    end
end
close all