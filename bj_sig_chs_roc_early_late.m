function bj_sig_chs_roc_early_late
%Written by Xing 14/08/13
%Compare PNE and contrast at half-height during first vs last 30% of
%sessions.
excludeSessHighSSE=0;
roving=0;
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
plotSlopeFig=1;
if plotSlopeFig==1
    animals=[{'blanco'} {'jack'}];
    areas=[{'v4_1'} {'v1_1'}];
    figELPNE=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.6, 0.8]); %
    set(figELPNE, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
    figELCHalf=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.6, 0.8]); %
    set(figELCHalf, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
    sampleContrast=30;
    allChInd=0;
    for animalInd=1:length(animals)
        animal=animals{animalInd};
        for areaInd=1:length(areas)
            area=areas{areaInd};
            channels=main_channels(animal,area);
            sessionNums=main_raw_sessions_final(animal,area,[],0);
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
                for sessionInd=1:length(sessionNums)
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
                    %for comparison across channels during early and late
                    %sessions:
                    allPNENeuro(sessionInd,chInd)=c50(1,chInd);
                    if max(yvals)<0.5%out of range
                        allCHalfNeuro(sessionInd,chInd)=100;
                    elseif min(yvals)>0.5
                        allCHalfNeuro(sessionInd,chInd)=0;
                    else
                        diffTemp=yvalsFine-((max(yvals)-min(yvals))/2+min(yvals));%find contrast at half height
                        [tempVal columnInd]=min(abs(diffTemp));
                        allCHalfNeuro(sessionInd,chInd)=xvalsFine(columnInd);
                    end
                end
            end
            figure(figELPNE);
            if roving==0
                subplot(2,2,animalInd+(areaInd-1)*2);
            elseif roving==1
            end
            earlyPNE=[];
            latePNE=[];
            numSess30=floor(size(ROCmat,1)*0.3);%number of sessions during first & last 30% of sessions
            for numSess30Ind=1:numSess30
                plot(numSess30Ind,allPNENeuro(numSess30Ind,:),'r.');hold on
                earlyPNE=[earlyPNE allPNENeuro(numSess30Ind,:)];
                plot(size(ROCmat,1)-numSess30Ind+1,allPNENeuro(size(ROCmat,1)-numSess30Ind+1,:),'b.');hold on
                latePNE=[latePNE allPNENeuro(size(ROCmat,1)-numSess30Ind+1,:)];
            end
            if areaInd==1
                title(['Monkey ',num2str(animalInd)]);
                if animalInd==1
                    xlabel('early/ late sessions');
                    ylabel('PNE');
                end
            end
            [h,pPNE(areaInd,animalInd),ci,stats]=ttest2(earlyPNE,latePNE);
            statsPNE(areaInd,animalInd)=stats;
            [pPNE30E(areaInd,animalInd),h,stats]=signrank(earlyPNE,30);
            statsPNEE(areaInd,animalInd)=stats;
            [pPNE30L(areaInd,animalInd),h,stats]=signrank(latePNE,30);
            statsPNEL(areaInd,animalInd)=stats;

            %then look at contrast at which half-height occurs
            figure(figELCHalf);
            if roving==0
                subplot(2,2,animalInd+(areaInd-1)*2);
            elseif roving==1
            end
            earlyCHalf=[];
            lateCHalf=[];
            numSess30=floor(size(ROCmat,1)*0.3);%number of sessions during first & last 30% of sessions
            for numSess30Ind=1:numSess30
                plot(numSess30Ind,allCHalfNeuro(numSess30Ind,:),'r.');hold on
                earlyCHalf=[earlyCHalf allCHalfNeuro(numSess30Ind,:)];
                plot(size(ROCmat,1)-numSess30Ind+1,allCHalfNeuro(size(ROCmat,1)-numSess30Ind+1,:),'b.');hold on
                lateCHalf=[lateCHalf allCHalfNeuro(size(ROCmat,1)-numSess30Ind+1,:)];
            end
            if areaInd==1
                title(['Monkey ',num2str(animalInd)]);
                if animalInd==1
                    xlabel('early/ late sessions');
                    ylabel('C_h_a_l_f_-_h_e_i_g_h_t');
                end
            end
            [h,pCHalf(areaInd,animalInd),ci,stats]=ttest2(earlyCHalf,lateCHalf);
            statsCHalf(areaInd,animalInd)=stats;
            [pCHalf30E(areaInd,animalInd),h,stats]=signrank(earlyPNE,30);
            statsCHalfE(areaInd,animalInd)=stats;
            [pCHalf30L(areaInd,animalInd),h,stats]=signrank(latePNE,30);
            statsCHalfL(areaInd,animalInd)=stats;
        end
    end
    figure(figELPNE);
    imagename='early_late_chs_PNE';
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
    
    figure(figELCHalf);
    imagename='early_late_chs_CHalf';
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

