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
readData=1;
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
                sessionNums=163:164;%139 is non-roving, 14 conditions. 140 onwards is roving, 36 conditions. 163 onwards is with flankers.
            end
        elseif strcmp(monkey,'blanco')
            splitSessions=[355 405 435];
            if strcmp(area,'v4_1')
                sessionNums=[307,308,311,312,313,314,316,317,318,320,321,322,323,324,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342];
            elseif strcmp(area,'v4_2')
                sessionNums=360:364;
            elseif strcmp(area,'v1')
                sessionNums=[343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359];
            elseif strcmp(area,'v1_2')
                sessionNums=[388:1:401 403 404 406:1:422 431:434 436:443 451:452 453:1:459];
            end
        end
        analysisFolderAppend=[];
        for i=1:length(sessionNums)
            if strcmp(monkey,'jack')
                cdText=['cd ','V:\thielelab\Groups\ThieleGroup\monkey_data\Jack\_jackgrid\j_',area,'_events_files\',num2str(sessionNums(i))];
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
                elseif sessionNums(i)>162&&sessionNums(i)<165%roving samples and 36 conditions
                    file_of_int='2223412f.1';%control task with jack's PL stimuli in blanco's V1 RF, plus flankers
                    analysisFolderAppend='_7';
                end
            elseif strcmp(monkey,'blanco')
                cdText=['cd ','F:\blanco\',area,'_events_files\',num2str(sessionNums(i))];
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
                if splitSess==1
                    cdText=['cd ','F:\blanco\',area,'_events_files\',num2str(sessionNums(i)),'_1'];
                end
                eval(cdText);
                if roving==0
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
%                         splitSess=0;
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
