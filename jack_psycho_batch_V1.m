function jack_psycho_batch_V1(animal,area,psychoOnly,plotLinear)
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
roving=0;
readData=0;
onExternalHD=0;
if onExternalHD==1
    rootFolder='G:\PL_backup_060413';
else
    rootFolder='F:';
end
% for monkeyCount=1:2
%     for areaCount=1:3
%         animal=monkeys{monkeyCount};
%         area=areas{areaCount};
        %no need to adjust:
        sampleContrast=30;
        splitSess=0;
        if strcmp(animal,'jack')
            appendText='j_';
%             if strcmp(area,'v4_1')
%                 sessionNums=[1:5];
%                 sessionNums=[24 26 28:49];
%                 sessionNums=24:50;
%             elseif strcmp(area,'v4_2')
%                 sessionNums=73:77;
%             elseif strcmp(area,'v1')
%                 sessionNums=51:72;
%             elseif strcmp(area,'v1_2')
%                 %         sessionNums=[78:93];
%                 %         sessionNums=94:115;
%                 %         sessionNums=116:119;
%                 sessionNums=78:119;
%                 sessionNums=116:120;
%             elseif strcmp(area,'v1_3')
%                 sessionNums=128:129;%126:127; v1_3 is bu/pl control task. 
%             elseif strcmp(area,'v1_4')%v1_4 is jack's control task with stim in blanco's v1 location
%                 sessionNums=184:189;%139 is non-roving, 14 conditions. 140 onwards is roving, 36 conditions. 163 to 184 are with flankers. 185 onwards are without flankers again.
%             end
            if psychoOnly==1
                sessionNums=main_raw_sessions_final_psycho(animal,area,[],1);
            else
                sessionNums=main_raw_sessions_final(animal,area,[],1);
            end
        elseif strcmp(animal,'blanco')
            splitSessions=[355 405 435];
            if psychoOnly==1
                sessionNums=main_raw_sessions_final_psycho(animal,area,[],1);
            else
                sessionNums=main_raw_sessions_final(animal,area,[],1);
            end
        end
        analysisFolderAppend=[];
        for i=1:length(sessionNums)
            if strcmp(animal,'jack')
                [sampleContrasts testContrasts]=area_metadata(area);
                if sessionNums(i)<22
                    file_of_int='2161392.1';
                    testContrast=[10 90];
                elseif sessionNums(i)==22
                    file_of_int='2161398.1';
                    testContrast=[10 15 20 25 35 40 60 90];
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
                    roving=1;
                elseif sessionNums(i)>162&&sessionNums(i)<185%roving samples and 36 conditions
                    file_of_int='2223412f.1';%control task with jack's PL stimuli in blanco's V1 RF, plus flankers
                    analysisFolderAppend='_7';
                    roving=1;
                elseif sessionNums(i)>184&&sessionNums(i)<190%roving samples and 36 conditions
                    file_of_int='225234912.1';%control task with jack's PL stimuli in blanco's V1 RF, plus flankers
                    analysisFolderAppend='_7';
                    roving=1;
                elseif sessionNums(i)==192%roving samples and 36 conditions
                    file_of_int='473914.1';%control task with jack's PL stimuli in blanco's V1 RF, plus flankers
                    roving=0;
                    testContrast=testContrasts;
                end
            elseif strcmp(animal,'blanco')
                [file_of_int,testContrast]=session_metadata(sessionNums(i),animal);
                if sessionNums(i)==304
                    file_of_int='21653090.1';
                    testContrast=[5 10 20 25 35 40 60 90];
                elseif sessionNums(i)==305
                    file_of_int='216136.2';
                    testContrast=[10 15 20 25 27 29 31 33 35 40 50 60];
                elseif sessionNums(i)>305&&sessionNums(i)<343
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
                    if strcmp(animal,'jack')
                        cdText=['cd ','V:\thielelab\Groups\ThieleGroup\monkey_data\Jack\_jackgrid\j_',area,'_events_files\',num2str(sessionNums(i))];
                    elseif strcmp(animal,'blanco')
                        cdText=['cd V:\thielelab\Groups\ThieleGroup\monkey_data\blanco\_grid\',area,'_events_files\',num2str(sessionNums(i))];
                    end
                    if splitSess==1
                        cdText=['cd V:\thielelab\Groups\ThieleGroup\monkey_data\blanco\_grid\',area,'_events_files\',num2str(sessionNums(i)),'_1'];
                    end
                    eval(cdText);
                    conditions=1:length(testContrast);
                    [vals]=jack_2target_psycho_EV_v3_mex(file_of_int,sampleContrast,testContrast,sessionNums(i),conditions,area,analysisFolderAppend,animal);
%                     [vals]=bj_2target_psycho_EV_v3_mex(file_of_int,sample
%                     Contrast,testContrast,sessionNums(i),conditions);%can
%                     use this for 304 & 305
                    if splitSess==1
                        cdText=['cd V:\thielelab\Groups\ThieleGroup\monkey_data\blanco\_grid\',area,'_events_files\',num2str(sessionNums(i)),'_2'];
                        eval(cdText);
                        [vals2]=jack_2target_psycho_EV_v3_mex(file_of_int,sampleContrast,testContrast,sessionNums(i),conditions,area,analysisFolderAppend,animal);
                        vals=[vals;vals2];
                    end
%                     analyseValsjexample(vals,testContrast,conditions,sampleContrast,animal,area,analysisFolderAppend,sessionNums(i));
                    analyseVals2(vals,testContrast,conditions,sampleContrast,animal,area,analysisFolderAppend,sessionNums(i),roving,onExternalHD,psychoOnly);
                    splitSess=0;
                elseif roving==1
                    testContrasts=[5 10 12 15 18 22 25 28 35 45 63 90;5 10 15 22 25 28 32 35 38 45 60 90;5 10 15 25 32 35 38 42 45 50 60 90];
                    sampleContrasts=[20 30 40];
                    allConditions=[13:1:24;1:1:12;25:1:36];
                    for sampleNum=1:3
                        sampleContrast=sampleContrasts(sampleNum);
                        testContrast=testContrasts(sampleNum,:);
                        conditions=allConditions(sampleNum,:);
                        if strcmp(animal,'jack')
                            cdText=['cd ','V:\thielelab\Groups\ThieleGroup\monkey_data\Jack\_jackgrid\j_',area,'_events_files\',num2str(sessionNums(i))];
                            eval(cdText);
                        elseif strcmp(animal,'blanco')
                            if splitSess==1
                                cdText=['cd ','V:\thielelab\Groups\ThieleGroup\monkey_data\blanco\_grid\',area,'_events_files\',num2str(sessionNums(i)),'_1'];
                            else
                                cdText=['cd ','V:\thielelab\Groups\ThieleGroup\monkey_data\blanco\_grid\',area,'_events_files\',num2str(sessionNums(i))];
                            end
                            eval(cdText);
                        end
                        %             [vals]=blanco_2target_psycho_EV_v4_mex(file_of_int,sampleContrast,testContrast,sessionNums(i),folder);
                        [vals]=jack_2target_psycho_EV_v3_mex(file_of_int,sampleContrast,testContrast,sessionNums(i),conditions,area,analysisFolderAppend,animal,roving);
                        if splitSess==1
                            cdText=['cd ','V:\thielelab\Groups\ThieleGroup\monkey_data\blanco\_grid\',area,'_events_files\',num2str(sessionNums(i)),'_2'];
                            eval(cdText);
                            [vals2]=jack_2target_psycho_EV_v3_mex(file_of_int,sampleContrast,testContrast,sessionNums(i),conditions,area,analysisFolderAppend,animal);
                            vals=[vals;vals2];
                        end
                        savevals{sampleNum}=analyseVals2(vals,testContrast,conditions,sampleContrast,animal,area,analysisFolderAppend,sessionNums(i),roving,onExternalHD,psychoOnly);
                    end
                end
                splitSess=0;
                valsFileName=['vals_',num2str(sessionNums(i)),'.mat'];
                valsFolder=fullfile(rootFolder,'PL','vals_perf',animal,valsFileName); %#ok<NASGU>
                saveText=['save ',valsFolder,' vals'];
                if roving==1
                    vals=savevals;
                end
                eval(saveText);
            end
        end
%     end
% end

taskCols=[204/256 153/256 256/256;255/256 222/256 173/256;154/256 205/256 50/256];%purple;orange;green
taskColsDark=[153/256 50/256 204/256;210/256 105/256 30/256;34/256 139/256 34/256];%purple;orange;green
taskColsVDark=[84/256 3/256 163/256;166/256 121/256 3/256;79/256 117/256 2/256];%purple;orange;green, for fitted line

plotRoving=0;
preflankerOnly=5;%1: flankers absent, 2: flankers present 3: control task, flankers absent 4: control task, flankers present 0: all roving 5: pre-flanker, flanker & post-flanker V1_2
if plotRoving==1
    if strcmp(animal,'jack')
        if preflankerOnly==1||preflankerOnly==2%show only part of the data
            areas={'v1_2'};
            areaText=[{'V1 RFs location (-0.7,-1.3)'}];
        elseif preflankerOnly==3||preflankerOnly==4%show only part of the control data
            areas={'v1_4'};
            areaText=[{'control location (-3.5,-3)'}];
        elseif preflankerOnly==0
            areas=[{'v1_2'} {'v1_4'}];
            areaText=[{'V1 RFs location (-0.7,-1.3)'} {'control location (-3.5,-3)'}];
        elseif preflankerOnly==5
            areas=[{'v1_2'}];
            areaText=[{'V1 RFs location (-0.7,-1.3)'} {'control location (-3.5,-3)'}];
        end
    elseif strcmp(animal,'blanco')
        areas={'v1_2'};
        areaText={'V1 RFs location (-3.5,-3)'};
    end
    if preflankerOnly==1||preflankerOnly==3
        fig=figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.2, 0.04, 0.5/(3-length(areas)), 0.9]);
    else
        fig=figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.2, 0.04, 0.75/(3-length(areas)), 0.9]);
    end
    set(fig,'PaperUnits','centimeters','PaperType','A4','PaperOrientation', 'portrait', 'PaperPosition', [0.63452 0.63452 31.5/(3-length(areas)) 28.41]);
    for areaInd=1:length(areas)
        area=areas{areaInd};
        [flankerType flankerSessions]=getFlankerType(animal,area);
        sampleContrasts=20:10:40;
        count=0;
        if preflankerOnly==1%show only the pre-flanker data
            flankerType=flankerType(1);
            flankerSessions=flankerSessions(1);
        elseif preflankerOnly==2%show only the pre-flanker and flanker data
            flankerType=flankerType(1:2);
            flankerSessions=flankerSessions(1:2);
        elseif preflankerOnly==3%show only the control pre-flanker data
            flankerType=flankerType(1);
            flankerSessions=flankerSessions(1);
        elseif preflankerOnly==4%show only the control pre-flanker and flanker data
            flankerType=flankerType(1:2);
            flankerSessions=flankerSessions(1:2);
        end
        for flankerInd=1:length(flankerType)
            taskCol=taskCols(flankerInd,:);
            taskColDark=taskColsDark(flankerInd,:);
            markerFaceTypes=[{'none'} {taskCol} {taskColDark}];
            for i=1:length(sampleContrasts)
                sampleContrast=sampleContrasts(i);
                if sampleContrast==30
                    markerFaceType=markerFaceTypes{2};
                    fittedLineCol=taskColsDark(flankerInd,:);
                elseif sampleContrast==20
                    markerFaceType=markerFaceTypes{1};
                    fittedLineCol=taskCols(flankerInd,:);
                elseif sampleContrast==40
                    markerFaceType=markerFaceTypes{3};
                    fittedLineCol=taskColsVDark(flankerInd,:);
                end
                loadText=['load ',rootFolder,'\PL\psycho_data\',animal,'\allMeanPerf\allMeanPerf_',area,'_',num2str(flankerInd),'_',num2str(sampleContrast),'.mat allMeanPerf'];
                eval(loadText)
                perf=allMeanPerf(:,2);
                PSE=allMeanPerf(:,15);
                slope=allMeanPerf(:,16);
                if length(areas)==1
                    subplot(3,length(areas),1);
                else
                    subplot(3,length(areas),areaInd);
                end
                plot(count+1:count+length(flankerSessions{flankerInd}),perf,'LineStyle','none','Marker','o','MarkerEdgeColor',taskColDark,'MarkerFaceColor',markerFaceType);hold on
                bj_linearexpo_fitting([count+1:count+length(flankerSessions{flankerInd})]',perf,0,0,'ROC',plotLinear,count,fittedLineCol);
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
                    subplot(3,length(areas),3);
                else
                    subplot(3,length(areas),areaInd+4);
                end
                plot(count+1:count+length(flankerSessions{flankerInd}),PSE,'LineStyle','none','Marker','o','MarkerEdgeColor',taskColDark,'MarkerFaceColor',markerFaceType);hold on
                bj_linearexpo_fitting([count+1:count+length(flankerSessions{flankerInd})]',PSE,0,0,'ROC',plotLinear,count,fittedLineCol);
                if length(areas)==2
                    ylim([20 50]);
                else
                    ylim([15 50]);
                end
                if areaInd==1
                    ylabel('PSE');
                end
                if length(areas)==1
                    subplot(3,length(areas),2);
                else
                    subplot(3,length(areas),areaInd+2);
                end
                plot(count+1:count+length(flankerSessions{flankerInd}),slope,'LineStyle','none','Marker','o','MarkerEdgeColor',taskColDark,'MarkerFaceColor',markerFaceType);hold on
                bj_linearexpo_fitting([count+1:count+length(flankerSessions{flankerInd})]',slope,0,0,'ROC',plotLinear,count,fittedLineCol);
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
    if strcmp(animal,'blanco')
%         if strcmp(area,'v1_4')
%             subplot(3,1,1);
%             xlabel('session number');
%             subplot(3,1,2);
%             ylim([0 9]);
%             subplot(3,1,3);
%             ylim([20 40]);
%         elseif strcmp(area,'v1_2')
%             subplot(3,1,1);
%             xlabel('session number');
%             ylim([0.5 0.9]);
%             subplot(3,1,2);
%             ylim([0 11]);
%             subplot(3,1,3);
%             ylim([20 45]);
%         end
        if preflankerOnly==1
            subplot(3,1,1);
            ylim([0.7 0.9]);
            xlim([0 35]);
            subplot(3,1,2);
            ylim([0 8]);
            xlim([0 35]);
            subplot(3,1,3);
            ylim([25 40]);
            xlim([0 35]);
        end
        if preflankerOnly==0
            subplot(3,1,1);
            xlim([0 57]);
            subplot(3,1,2);
            xlim([0 57]);
            subplot(3,1,3);
            xlim([0 57]);
        end
    end
    if strcmp(animal,'jack')
        if strcmp(area,'v1_4')&&preflankerOnly~=0
            subplot(3,1,1);
            xlabel('session number');
            subplot(3,1,2);
            ylim([0 9]);
            subplot(3,1,3);
            ylim([20 40]);
        elseif strcmp(area,'v1_2')&&preflankerOnly==2
            subplot(3,1,1);
            xlabel('session number');
            ylim([0.5 0.9]);
            subplot(3,1,2);
            ylim([0 11]);
            subplot(3,1,3);
            ylim([20 45]);
        elseif strcmp(area,'v1_4')&&preflankerOnly==0
            subplot(3,2,1);
            xlabel('session number');
            xlim([0 43]);
            subplot(3,2,2);
            xlim([0 51]);
            subplot(3,2,3);
            ylim([0 9]);
            xlim([0 43]);
            subplot(3,2,4);
            ylim([0 9]);
            xlim([0 51]);
            subplot(3,2,5);
            ylim([20 45]);
            xlim([0 43]);
            subplot(3,2,6);
            ylim([20 45]);
            xlim([0 51]);
        elseif strcmp(area,'v1_4')&&preflankerOnly==4
            subplot(3,1,1);
            xlabel('session number');
            xlim([0 43]);
            subplot(3,1,2);
            xlim([0 43]);
            subplot(3,1,3);
            ylim([20 40]);
            xlim([0 43]);
        elseif strcmp(area,'v1_4')&&preflankerOnly==5
            subplot(3,1,1);
            xlabel('session number');
            xlim([0 43]);
            ylim([0.5 0.9]);
            subplot(3,1,2);
            xlim([0 43]);
            ylim([0 11]);
            subplot(3,1,3);
            ylim([20 45]);
            xlim([0 43]);
        end
        if preflankerOnly==1
            subplot(3,1,1);
            ylim([0.7 0.9]);
            subplot(3,1,2);
            ylim([0 11]);
            subplot(3,1,3);
            ylim([20 40]);
        end
    end
    subplot(3,length(areas),length(areas));
    ylims=get(gca,'YLim');
    xlims=get(gca,'XLim');
    sampleTexts=[{'20% sample'} {'30% sample'} {'40% sample'}];
    flankerTexts=[{'no flankers'} {'flankers'} {'no flankers'}];
    markerFaceTypes=[{'none'} {[180/256 180/256 180/256]} {'k'}];
    markerFaceTypesPreflankerOnly=[{'none'} {[204/256 153/256 256/256]} {[153/256 50/256 204/256]}];%purple
    markerFaceTypesFlankerOnly=[{'none'} {[255/256 222/256 173/256]} {[210/256 105/256 30/256]}];%green
    taskColsMedium=[147/256 112/256 219/256;255/256 165/256 0/256;154/256 205/256 50/256];%purple;orange;green
    for i=1:3
        if preflankerOnly==0
            if length(areas)==2
                keyText=[sampleTexts{i},'         ',flankerTexts{i}];
            else
                keyText=[sampleTexts{i},'            ',flankerTexts{i}];
            end
            plot(xlims(1)+(xlims(2)-xlims(1))/20,ylims(1)+(ylims(2)-ylims(1))/12*(i),'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor',markerFaceTypes{i});
            plot(xlims(1)+(xlims(2)-xlims(1))/3.6,ylims(1)+(ylims(2)-ylims(1))/12*(i),'LineStyle','none','Marker','o','MarkerEdgeColor',taskColsMedium(i,:),'MarkerFaceColor',taskColsMedium(i,:));
            text('Position',[xlims(1)+(xlims(2)-xlims(1))/15,ylims(1)+(ylims(2)-ylims(1))/12*(i)],'LineStyle','none','FontSize',9,'String',keyText,'Color','k');
        elseif preflankerOnly==1||preflankerOnly==2||preflankerOnly==3||preflankerOnly==4
            if length(areas)==2
                keyText=[sampleTexts{i}];
            else
                keyText=[sampleTexts{i}];
            end
            plot(xlims(1)+(xlims(2)-xlims(1))/20,ylims(1)+(ylims(2)-ylims(1))/12*(i),'LineStyle','none','Marker','o','MarkerEdgeColor',taskColsDark(1,:),'MarkerFaceColor',markerFaceTypesPreflankerOnly{i});
            text('Position',[xlims(1)+(xlims(2)-xlims(1))/15,ylims(1)+(ylims(2)-ylims(1))/12*(i)],'LineStyle','none','FontSize',9,'String',keyText,'Color','k');
            if preflankerOnly==2||preflankerOnly==4
                plot(xlims(1)+xlims(2)-(xlims(2)-xlims(1))/4,ylims(1)+(ylims(2)-ylims(1))/12*(i),'LineStyle','none','Marker','o','MarkerEdgeColor',taskColsDark(2,:),'MarkerFaceColor',markerFaceTypesFlankerOnly{i});
                text('Position',[xlims(1)+xlims(2)-(xlims(2)-xlims(1))/5,ylims(1)+(ylims(2)-ylims(1))/12*(i)],'LineStyle','none','FontSize',9,'String',keyText,'Color','k');
            end
        end
    end
    formats=[{'eps'} {'png'}];
    formatsCommand=[{'epsc'} {'png'}];
    for j=1:2
        if length(areas)==2
            printFileName=[rootFolder,'\PL\psycho_data\',animal,'\PC_PSE_slope_',areas{1},'_',areas{2},'.',formats{j}];
        elseif length(areas)==1
            if preflankerOnly==0
                printFileName=[rootFolder,'\PL\psycho_data\',animal,'\PC_PSE_slope_',areas{1},'.',formats{j}];
            elseif preflankerOnly==1
                printFileName=[rootFolder,'\PL\psycho_data\',animal,'\PC_PSE_slope_',areas{1},'preflankers.',formats{j}];
            elseif preflankerOnly==2
                printFileName=[rootFolder,'\PL\psycho_data\',animal,'\PC_PSE_slope_v1_2_1_v1_2_2.',formats{j}];
            elseif preflankerOnly==3
                printFileName=[rootFolder,'\PL\psycho_data\',animal,'\PC_PSE_slope_',areas{1},'preflankers.',formats{j}];
            elseif preflankerOnly==4
                printFileName=[rootFolder,'\PL\psycho_data\',animal,'\PC_PSE_slope_v1_4_1_v1_4_2.',formats{j}];
            end
            print(sprintf('-d%s',formatsCommand{j}),'-r300',printFileName)
        end
    end
end
close all

%analyse performance between pairs of sample contrasts
plotPairSample=1;
preflankerOnly=1;%1: flankers absent, 2: flankers present 
if plotPairSample==1
    animals=[{'blanco'} {'jack'}];
    animalText=[{'Monkey 1'} {'Monkey 2'}];
    if preflankerOnly==1
        fig=figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.04, 0.04, 0.4, 0.35]);
    else
        fig=figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.04, 0.04, 0.4 0.35]);
    end
    set(fig,'PaperUnits','centimeters','PaperType','A4','PaperOrientation', 'portrait');
    for animalInd=1:length(animals)
        animal=animals{animalInd};
        markerCols='kcm';%20 vs 30; 30 vs 40; 20 vs 40
        comparisons=[1 2;2 3;1 3];
        [flankerType flankerSessions]=getFlankerType(animal,area);
        sampleContrasts=20:10:40;
        if preflankerOnly==1%show only the pre-flanker data
            flankerType=flankerType(1);
            flankerSessions=flankerSessions(1);
        elseif preflankerOnly==2%show only the pre-flanker data
            flankerType=flankerType(2);
            flankerSessions=flankerSessions(2);
        end
%         markerFaceTypes=[{'none'} {taskCol} {taskColDark}];
        perf=[];PSE=[];slope=[];
        for i=1:length(sampleContrasts)
            sampleContrast=sampleContrasts(i);
            loadText=['load ',rootFolder,'\PL\psycho_data\',animal,'\allMeanPerf\allMeanPerf_v1_2_',num2str(flankerType),'_',num2str(sampleContrast),'.mat allMeanPerf'];
            eval(loadText)
            perf(i,:)=allMeanPerf(:,2);
            PSE(i,:)=allMeanPerf(:,15);
            slope(i,:)=allMeanPerf(:,16);
        end
        for comparisonInd=1:3
            subplot(1,2,animalInd);
%             if sampleContrast==30
%                 markerFaceType=markerFaceTypes{2};
%             elseif sampleContrast==20
%                 markerFaceType=markerFaceTypes{1};
%             elseif sampleContrast==40
%                 markerFaceType=markerFaceTypes{3};
%             end
            plot(perf(comparisons(comparisonInd,2),:),perf(comparisons(comparisonInd,1),:),'LineStyle','none','Marker','o','MarkerEdgeColor',markerCols(comparisonInd),'MarkerFaceColor',markerCols(comparisonInd));hold on
            hold on
            [rPairSamp(animalInd,comparisonInd),pPairSamp(animalInd,comparisonInd)]=corr(perf(comparisons(comparisonInd,1),:)',perf(comparisons(comparisonInd,2),:)');
            dfs(animalInd,comparisonInd)=length(perf(comparisons(comparisonInd,1),:))-2;
            title(animalText{animalInd});
            if animalInd==1
                ylabel('performance');
                xlabel('performance');
            end
            %perform partial correlation to control for possible
            %effect of session number:
            [r p]=partialcorr([zscore(perf(comparisons(comparisonInd,2),:))' zscore(perf(comparisons(comparisonInd,1),:))'],[1:size(perf,2)]','type','Spearman');
            rPartial(animalInd,comparisonInd)=r(2);
            pPartial(animalInd,comparisonInd)=p(2);
        end
        if preflankerOnly==1
            subplot(1,2,1)
            xlim([0.75 0.9]);
            ylim([0.7 0.9]);
            subplot(1,2,2)
            ylim([0.7 0.9]);
        elseif preflankerOnly==2
            subplot(1,2,1)
            xlim([0.7 0.95]);
            ylim([0.6 0.92]);
            subplot(1,2,2)
            ylim([0.5 0.85]);
        end
        subplot(1,2,1);
        ylims=get(gca,'YLim');
        xlims=get(gca,'XLim');
        sampleTexts=[{'20% vs 30%'} {'30% vs 40%'} {'20% vs 40%'}];
        flankerTexts=[{'no flankers'} {'flankers'} {'no flankers'}];
        markerFaceTypes=[{'none'} {[180/256 180/256 180/256]} {'k'}];
        markerFaceTypesPreflankerOnly=[{'none'} {[204/256 153/256 256/256]} {[153/256 50/256 204/256]}];
        taskColsMedium=[147/256 112/256 219/256;255/256 165/256 0/256;154/256 205/256 50/256];%purple;orange;green
        for i=1:3
            if preflankerOnly==0
                if length(areas)==2
                    keyText=[sampleTexts{i}];
                else
                    keyText=[sampleTexts{i}];
                end
                plot(xlims(1)+(xlims(2)-xlims(1))/20,ylims(1)+(ylims(2)-ylims(1))/10*(4-i),'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor',markerFaceTypes{i});
                plot(xlims(1)+(xlims(2)-xlims(1))/3.6,ylims(1)+(ylims(2)-ylims(1))/10*(4-i),'LineStyle','none','Marker','o','MarkerEdgeColor',taskColsMedium(i,:),'MarkerFaceColor',taskColsMedium(i,:));
                text('Position',[xlims(1)+(xlims(2)-xlims(1))/15,ylims(1)+(ylims(2)-ylims(1))/10*(4-i)],'LineStyle','none','FontSize',9,'String',keyText,'Color','k');
            elseif preflankerOnly==1||preflankerOnly==2
                keyText=[sampleTexts{i}];
                plot(xlims(1)+(xlims(2)-xlims(1))/20,ylims(1)+(ylims(2)-ylims(1))/15*(4-i),'LineStyle','none','Marker','o','MarkerEdgeColor',markerCols(i),'MarkerFaceColor',markerCols(i));
                text('Position',[xlims(1)+(xlims(2)-xlims(1))/15,ylims(1)+(ylims(2)-ylims(1))/15*(4-i)],'LineStyle','none','FontSize',9,'String',keyText,'Color','k');
            end
        end
    end
    formats=[{'eps'} {'png'}];
    formatsCommand=[{'epsc'} {'png'}];
    for j=1:2
        if preflankerOnly==0
            printFileName=[rootFolder,'\PL\psycho_data\roving_sample_PC_',areas{1},'.',formats{j}];
        elseif preflankerOnly==1
            printFileName=[rootFolder,'\PL\psycho_data\roving_sample_PC_v1_2_1.',formats{j}];
        elseif preflankerOnly==2
            printFileName=[rootFolder,'\PL\psycho_data\roving_sample_PC_v1_2_2.',formats{j}];
        end
        set(gcf, 'PaperPositionMode', 'auto');
        print(sprintf('-d%s',formatsCommand{j}),'-r300',printFileName)
    end
end
close all