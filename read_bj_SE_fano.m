function read_bj_SE_fano(roving)
%Written by Xing 16/09/10
minusSpontan=0;
markerTexts='+x';

animalTexts=[{'subject B'} {'subject J'}];
animals=[{'blanco'} {'jack'}];
if roving==0
    areaTexts=[{'V4'} {'V1'}];
    areas=[{'v4_1'} {'v1_1'}];
elseif roving==1
    areaTexts=[{'V1_2_1'} {'V1_2_2'} {'V1_2_3'}];
    areas=[{'v1_2_1'} {'v1_2_2'} {'v1_2_3'}];
end
epoch=4;
for animalInd=1:2
    animal=animals{animalInd};
    for areaInd=1:2
        area=areas{areaInd};
        %combine PSTH activity values across channels and sessions:
        if nargin<3 || isempty(channels)
            channels = main_channels(animal,area);
        end
        if nargin<4 || isempty(sessionNums)
            sessionNums = main_raw_sessions_final(animal,area,[],0);
        end
        [sampleContrasts allTestContrasts]=area_metadata(area);
        if minusSpontan==1
            subfolder=['fano',num2str(epoch),'_mspontan_',area];%folder for stimulus-evoked responses minus spontaneous activity levels
        elseif minusSpontan==0
            subfolder=['fano',num2str(epoch),'_wspontan_',area];%folder for stimulus-evoked responses minus spontaneous activity levels
        end
        allFanoSessions=[];
        colmapText=colormap(jet(size(allTestContrasts,2)));
        colmapText=[colmapText(1,:);132/255 22/255 216/255;202/255 65/255 223/255;colmapText(3:7,:);157/255 212/255 61/255;colmapText(10:12,:);178/255 111/255 12/255;colmapText(end,:)];
        analyseSeparateChannels=1;
        if analyseSeparateChannels==1
            sigCoefs=[];
            increase=zeros(1,length(sampleContrasts));
            decrease=zeros(1,length(sampleContrasts));
            for sampleContrastInd=1:length(sampleContrasts)%combine across sample contrasts (if >1) because highest test contrast is always the same (60% in V4, 90% in V1)
                sampleContrast=sampleContrasts(sampleContrastInd);
                if roving==0&&animalInd==1&&areaInd==1
                    fig3=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.4, 0.6]); %
                    set(fig3, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                elseif roving==1&&animalInd==1
                    fig3(areaInd)=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.4, 0.6]); %
                    set(fig3(areaInd), 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                end
                fig1=figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.05,0.05, 0.9, 0.9]);
                for i=1:length(channels)
                    fanoSessions=[];
                    for j=1:length(sessionNums)
                        matFanoName=[num2str(channels(i)),'_',num2str(sessionNums(j)),'_',num2str(sampleContrast),'_',area,'_fano_vals.mat'];
                        matFanoPath=fullfile('F:','PL','fano',animal,subfolder,matFanoName);
                        if exist(matFanoPath,'file')
                            loadText=['load ',matFanoPath,' fano_vals'];
                            eval(loadText);
                            fanoSessions=[fanoSessions;fano_vals];
                            allFanoSessions=[allFanoSessions;channels(i) sessionNums(j) fano_vals];
                        end
                    end
                    subplot(ceil(length(channels)/5),5,i);
                    markerText=markerTexts(2);markerS=8;
                    for cond=1:size(allTestContrasts,2)
                        plot(1:size(fanoSessions,1),fanoSessions(:,cond),'Color',colmapText(cond,:),'LineStyle','none','Marker',markerText,'MarkerFaceColor',colmapText(cond,:),'MarkerEdgeColor',colmapText(cond,:),'MarkerSize',markerS);hold on%'MarkerFaceColor',[1/i 1/i 1/i],
                        hold on
                        a=[(1:size(fanoSessions,1))' fanoSessions(:,cond)];%check for change in slope with time- expect a +ve correlation, if any
                        [coefficients1 p1]=corrcoef(a);
                        allfanoCoefs(i,cond)={[coefficients1(1,2) p1(1,2)]};
                        if p1(1,2)<0.05
                            sigCoefs=[sigCoefs;channels(i) sampleContrast cond coefficients1(1,2) p1(1,2) fano_vals];
                        end
                    end
                    
                    %run ANOVA to check whether there is a main effect of
                    %period (early vs late), while taking test contrast as
                    %a second factor:
                    earlylateIndividualChFano=[];
                    condArr=[];
                    earlylateArr=[];
                    num30Sess=floor(0.3*size(fanoSessions,1));%number of sessions in first and last half
                    for sessCount=1:num30Sess%early sessions
                        earlylateIndividualChFano=[earlylateIndividualChFano fanoSessions(sessCount,:)];
                        condArr=[condArr 1:size(fanoSessions,2)];
                        earlylateArr=[earlylateArr zeros(1,size(fanoSessions,2))+1];%1 for early; 2 for late
                    end
                    for sessCount=size(fanoSessions,1)-num30Sess+1:size(fanoSessions,1)%late sessions
                        earlylateIndividualChFano=[earlylateIndividualChFano fanoSessions(sessCount,:)];
                        condArr=[condArr 1:size(fanoSessions,2)];
                        earlylateArr=[earlylateArr zeros(1,size(fanoSessions,2))+2];%1 for early; 2 for late
                    end
                    [p,t,stats]=anovan(earlylateIndividualChFano,{earlylateArr,condArr},'model','full');%check for differences in fano factor values based on 2 factors: training period (early vs late sessions) and condition number
                    allpIndCh{i}=p;%row 1: main effect of training period; row 2: main effect of condition; row 3: interaction between the two factors
                    pIndCh(i)=p(1);
                    tIndCh{i}=t;
                    statsIndCh{i}=stats;
                    statsColumnIndCh(i)=stats.dfe;
                    if pIndCh(i)<.05
                        c=multcompare(stats,'dimension',1);
                        if c(4)>0
                            decrease(sampleContrastInd)=decrease(sampleContrastInd)+1;
                        elseif c(4)<0%negative value means that early session values were smaller (c(4) equals to mean of group 1 minus mean of group 2
                            increase(sampleContrastInd)=increase(sampleContrastInd)+1;%so if value is negative, group 2 mean was higher than group 1 mean, i.e. later sessions had higher mean than early sessions
                        end
                    end
                    
                    if strcmp(animal,'blanco')&&strcmp(area,'v1_2')
                        ylim([0 3]);
                    end
                    %         xlabel('session');
                    %         ylabel('Fano factor');
                    title(num2str(channels(i)));
                end
                saveMeanImageName=[area,'_',animal,'_epoch',num2str(epoch),'_between_channels_vs_session_',num2str(sampleContrast)];
                saveMeanImageFolder=fullfile('F:','PL','fano',animal);
                saveMeanImagePath=fullfile(saveMeanImageFolder,saveMeanImageName);
                printtext=sprintf('print -dpng -r600 %s.png',saveMeanImagePath);
                set(gcf,'PaperPositionMode','auto')
                eval(printtext);
                matGrandName=['fanoList_',area,'_',num2str(sampleContrast)];
                matGrandPath=fullfile('F:','PL','fano',animal,subfolder,matGrandName);
                saveText=['save ',matGrandPath,' allFanoSessions'];
                eval(saveText);
                matCoefName=['fanoCoefsList_',area,'_',num2str(sampleContrast)];
                matCoefPath=fullfile('F:','PL','fano',animal,subfolder,matCoefName);
                saveText=['save ',matCoefPath,' allfanoCoefs'];
                eval(saveText);
                increase(sampleContrastInd)
                decrease(sampleContrastInd)
                totalSig(sampleContrastInd)=increase(sampleContrastInd)+decrease(sampleContrastInd)
            end
            
            matSigCoefsName=['sigCoefs_',area];
            matSigCoefsPath=fullfile('F:','PL','fano',animal,subfolder,matSigCoefsName);
            saveText=['save ',matSigCoefsPath,' sigCoefs totalSig increase decrease'];
            eval(saveText);
            close all hidden
        end
        
        analyseMeanChannels=0;
        if analyseMeanChannels==1
            %carry out correlation analysis on fano factors averaged across channels:
            meanSigCoefs=[];
            for sampleContrastInd=1:length(sampleContrasts)%combine across sample contrasts (if >1) because highest test contrast is always the same (60% in V4, 90% in V1)
                sampleContrast=sampleContrasts(sampleContrastInd);
                fig2=figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.05,0.05, 0.9, 0.9]);
                meanFanoSessions=[];
                earlyFanoSessions=[];
                lateFanoSessions=[];
                numSessPeriod=floor(length(sessionNums)*0.3);
                for j=1:length(sessionNums)
                    fanoSessions=[];
                    for i=1:length(channels)
                        matFanoName=[num2str(channels(i)),'_',num2str(sessionNums(j)),'_',num2str(sampleContrast),'_',area,'_fano_vals.mat'];
                        matFanoPath=fullfile('F:','PL','fano',animal,subfolder,matFanoName);
                        if exist(matFanoPath,'file')
                            loadText=['load ',matFanoPath,' fano_vals'];
                            eval(loadText);
                            fanoSessions=[fanoSessions;fano_vals];
                            if j<=numSessPeriod
                                earlyFanoSessions=[earlyFanoSessions;fano_vals];%combine across all channels and early sessions for calculations of SD
                            elseif j>length(sessionNums)-numSessPeriod
                                lateFanoSessions=[lateFanoSessions;fano_vals];%combine across all channels and late sessions for calculations of SD
                            end
                        end
                    end
                    meanFanoSessions=[meanFanoSessions;sessionNums(j) mean(fanoSessions,1)];
                end
                markerText=markerTexts(2);markerS=8;
                meanELFanoSessions=[];%compare early and late sessions
                stdELFanoSessions=[];
                condArr=[];
                for cond=1:size(allTestContrasts,2)
                    plot(1:size(meanFanoSessions,1),meanFanoSessions(:,cond+1),'Color',colmapText(cond,:),'LineStyle','none','Marker',markerText,'MarkerFaceColor',colmapText(cond,:),'MarkerEdgeColor',colmapText(cond,:),'MarkerSize',markerS);hold on%'MarkerFaceColor',[1/i 1/i 1/i],
                    hold on
                    a=[(1:size(meanFanoSessions,1))' meanFanoSessions(:,cond+1)];%check for change in slope with time- expect a +ve correlation, if any. Note: 'cond+1' because first column contains session number
                    [coefficients1 p1]=corrcoef(a);
                    allMeanfanoCoefs(cond)={[coefficients1(1,2) p1(1,2)]};
                    if p1(1,2)<0.05
                        meanSigCoefs=[meanSigCoefs;sampleContrast cond coefficients1(1,2) p1(1,2) fano_vals];
                    end
                    [h,p,ci,stats]=ttest(earlyFanoSessions(:,cond),lateFanoSessions(:,cond))
                    p2(cond)=p;
                    stats2{cond}=stats;
                    meanELFanoSessions(cond,1:2)=[mean(earlyFanoSessions(:,cond)) mean(lateFanoSessions(:,cond))];%calculate mean and SD across channels and sessions in early and late phases 
                    stdELFanoSessions(cond,1:2)=[std(earlyFanoSessions(:,cond)) std(lateFanoSessions(:,cond))];   
                    condArr=[condArr ones(size(earlyFanoSessions,1),1)*cond]; 
                end 
                
                %
                
                figure(fig2);
                earlylateFanoSessions=[earlyFanoSessions;lateFanoSessions];
                condArr=[condArr;condArr];%array containing condition numbers for early (upper half) and late (lower half) sessions
                earlylateArr=[earlyFanoSessions*0+1;earlyFanoSessions*0+2];%array of ones, representing first 30% of sessions, followed by twos, representing last 30%
                earlylateFanoSessions=reshape(earlylateFanoSessions,1,size(earlylateFanoSessions,1)*size(earlylateFanoSessions,2));
                condArr=reshape(condArr,1,size(earlylateFanoSessions,1)*size(earlylateFanoSessions,2));
                earlylateArr=reshape(earlylateArr,1,size(earlylateFanoSessions,1)*size(earlylateFanoSessions,2));
                [p,t,stats]=anovan(earlylateFanoSessions,{earlylateArr,condArr},'model','full');%check for differences in fano factor values based on 2 factors: training period (early vs late sessions) and condition number
                p3{animalInd+2*(areaInd-1),sampleContrastInd}=p;%row 1: main effect of training period; row 2: main effect of condition; row 3: interaction between the two factors
                p4(animalInd+2*(areaInd-1),sampleContrastInd)=p(1);%row 1: main effect of training period; row 2: main effect of condition; row 3: interaction between the two factors
                t3{animalInd+2*(areaInd-1),sampleContrastInd}=t;
                stats3{animalInd+2*(areaInd-1),sampleContrastInd}=stats;
                statsColumn(animalInd+2*(areaInd-1),sampleContrastInd)=stats.dfe;
                if areaInd==1
                    title(['Monkey ',num2str(animalInd)]);
                    xlim([0 70]);
                    if animalInd==1
                        xlabel('contrast (%)');
                        ylabel('Fano factor');
                    end
                end
                figure;
                [cRTceANOVA{animalInd+2*(areaInd-1),sampleContrastInd},m,h]=multcompare(stats,'dimension',[1 2])
                Fs(animalInd+2*(areaInd-1),sampleContrastInd)=t{2,6};
                if strcmp(animal,'blanco')&&strcmp(area,'v1_2')
                    ylim([0 3]);
                end
                %         xlabel('session');
                %         ylabel('Fano factor');
                title(num2str(channels(i)));
                saveallMeanImageName=[area,'_',animal,'_epoch',num2str(epoch),'_mean_channels_vs_session_',num2str(sampleContrast)];
                saveallMeanImageFolder=fullfile('F:','PL','fano',animal);
                saveallMeanImagePath=fullfile(saveallMeanImageFolder,saveallMeanImageName);
                printtext=sprintf('print -dpng -r600 %s.png',saveallMeanImagePath);
                set(gcf,'PaperPositionMode','auto')
                eval(printtext);
                matGrandName=['meanfanoList_',area,'_',num2str(sampleContrast)];
                matGrandPath=fullfile('F:','PL','fano',animal,subfolder,matGrandName);
                saveText=['save ',matGrandPath,' meanFanoSessions'];
                eval(saveText);
                matAllCoefName=['fanoList_',area,'_',num2str(sampleContrast)];
                matAllCoefPath=fullfile('F:','PL','fano',animal,subfolder,matAllCoefName);
                saveText=['save ',matAllCoefPath,' allMeanfanoCoefs'];
                eval(saveText);
                if roving==0
                    figure(fig3(sampleContrastInd));
                    subplot(2,2,animalInd+2*(areaInd-1));
                elseif roving==1
                    figure(fig3(areaInd));
                    subplot(3,2,animalInd+2*(sampleContrastInd-1));
                end
                for cond=1:size(allTestContrasts,2)
%                     if p2(cond)<.05
                        MarkerType='o';
%                     else
%                         MarkerType='.';
%                     end
                    plot(allTestContrasts(sampleContrastInd,cond),meanELFanoSessions(cond,1),'LineStyle','none','Marker',MarkerType,'Color','k');hold on
                    plot(allTestContrasts(sampleContrastInd,cond),meanELFanoSessions(cond,2),'LineStyle','none','Marker',MarkerType,'Color','r');hold on
                end
                errorbar(allTestContrasts(sampleContrastInd,:),meanELFanoSessions(:,1),stdELFanoSessions(:,1),'LineStyle','none','Marker','none','Color','k');hold on
                errorbar(allTestContrasts(sampleContrastInd,:),meanELFanoSessions(:,2),stdELFanoSessions(:,2),'LineStyle','none','Marker','none','Color','r');hold on
                if roving==0&&areaInd==1&&sampleContrastInd==1
                    title(['Monkey ',num2str(animalInd)]);
                    xlim([0 70]);
                    if animalInd==1
                        xlabel('contrast (%)');
                        ylabel('Fano factor');
                    end
                elseif roving==1&&sampleContrastInd==1
                    title(['Monkey ',num2str(animalInd)]);
                    if animalInd==1
                        xlabel('contrast (%)');
                        ylabel('Fano factor');
                    end
                end
            end            
            matSigCoefsName=['meanSigCoefs_',area];
            matSigCoefsPath=fullfile('F:','PL','fano',animal,subfolder,matSigCoefsName);
            saveText=['save ',matSigCoefsPath,' meanSigCoefs'];
            eval(saveText);
%             close all
        end        
    end
end
anovaTable=[statsColumn(:,1) Fs(:,1) p4(:,1) statsColumn(:,2) Fs(:,2) p4(:,2) statsColumn(:,3) Fs(:,3) p4(:,3)];
figure(fig3);
subplot(2,2,3);ylim([0.5 2.2])

areas={'v1_2'};
figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.9, 0.1, 0.7, 0.4]);
belowZero=zeros(length(animals),length(areas));
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        [sampleContrasts allTestContrasts]=area_metadata(area);
        colmapText=colormap(jet(size(allTestContrasts,2)));
        colmapText=[colmapText(1,:);132/255 22/255 216/255;202/255 65/255 223/255;colmapText(3:7,:);157/255 212/255 61/255;colmapText(10:12,:);178/255 111/255 12/255;colmapText(end,:)];
        markerText=markerTexts(2);markerS=8;
        if minusSpontan==1
            subfolder=['fano',num2str(epoch),'_mspontan_',area];%folder for stimulus-evoked responses minus spontaneous activity levels
        elseif minusSpontan==0
            subfolder=['fano',num2str(epoch),'_wspontan_',area];%folder for stimulus-evoked responses minus spontaneous activity levels
        end
        matSigCoefsName=['sigCoefs_',area];
        matSigCoefsPath=fullfile('F:','PL','fano',animal,subfolder,matSigCoefsName);
        loadText=['load ',matSigCoefsPath,' sigCoefs'];
        eval(loadText);
        subplotInd=subplot(2,2,animalInd+2*(areaInd-1));
        for rowInd=1:size(sigCoefs,1)
            plot(1,sigCoefs(4,sigCoefs(rowInd,3)),'Color',colmapText(sigCoefs(rowInd,3),:),'LineStyle','none','Marker',markerText,'MarkerFaceColor',colmapText(sigCoefs(rowInd,3),:),'MarkerEdgeColor',colmapText(sigCoefs(rowInd,3),:),'MarkerSize',markerS);hold on%'MarkerFaceColor',[1/i 1/i 1/i],
            if sigCoefs(rowInd,4)<0
                belowZero(animalInd+2*(areaInd-1))=belowZero(animalInd+2*(areaInd-1))+1;
            end
        end
    end
end
belowZero