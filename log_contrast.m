function log_contrast
samePlot=1;
originalSessList=1;
errorType='mle';
errorType='least_square';
shiftPSEContrast=1;%check what happens to logistic function when midpoint changes
shiftSlopeContrast=0;%check what happens to logistic function when slope changes
shiftLogContrastLogistic=0;%shift the midpoint for log(contrast), check whether slope of fitted logistic function changes
shiftLogContrastWeibull=0;%shift the midpoint for log(contrast), check whether slope of fitted weibull function changes
animals=[{'blanco'} {'jack'}];
areas=[{'v4_1'} {'v1_1'}];
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        allSessR=[];
        folderName=['F:\PL\sample_test_activity\',animal,'_',area];
        channels = main_channels(animal,area);
        sessionNums = main_raw_sessions_final(animal,area,[],0);
        [sampleContrasts testContrasts]=area_metadata(area);
        if originalSessList==1
            if strcmp(animal,'blanco')&&strcmp(area,'v4_1')
                sessionNums=[307 308 311 313 314 318 320 321 329 330 331:1:341];% blanco V4
            end
        end
        for sampleInd=1:length(sampleContrasts)
            sampleContrast=sampleContrasts(sampleInd);
            testContrast=testContrasts(sampleInd,:);
            for chInd=1:length(channels)
                loadText=['F:\PL\ROC\',animal,'\',area,'\ROC_Ch',num2str(channels(chInd)),'_',num2str(sampleContrast),'_1024_to_1536.mat'];
                load(loadText)
                sessList=[];%arrange ROCmat rows according to ascending session number
                for sessionInd=1:size(ROCmat,1)
                    sessList=[sessList ROCmat{sessionInd,1}];                
                end
                [sessListSorted ind]=sort(sessList);
                ROCmat2=cell(size(ROCmat,1),3);
                for sessionInd=1:size(ROCmat,1)
                    ROCmat2(sessionInd,:)=[ROCmat(ind(sessionInd),:)];                
                end
                ROCmat=ROCmat2;
                saveText=['F:\PL\ROC\',animal,'\',area,'\ROC_Ch',num2str(channels(chInd)),'_',num2str(sampleContrast),'_1024_to_1536.mat ROCmat'];
                save(saveText)
                allX=[];
                PSE=[];
                allX2=[];
                PSE2=[];
                figure1=figure('Color',[1,1,1],'Units','Normalized','Position',[0, 0, 1, 1]);
                if samePlot==0
                    figure2=figure('Color',[1,1,1],'Units','Normalized','Position',[0, 0, 1, 1]);
                end
                for sessionInd=1:size(ROCmat,1)
                    dummy=ROCmat{sessionInd,3};
                    figure(figure1)
                    if samePlot==1
                        subplot(1,2,1);
                    else
                        subplot(ceil(size(ROCmat,1)/6),6,sessionInd);
                    end
                    plot(testContrast,dummy,'Color',[sessionInd/size(ROCmat,1) 0 sessionInd/size(ROCmat,1)],'Marker','o','MarkerFaceColor',[sessionInd/size(ROCmat,1) 0 sessionInd/size(ROCmat,1)],'LineStyle','none');hold on
                    X0=[0.2 sampleContrast];
                    X=fminsearch(@fit_logistic,X0,[],testContrast,dummy,[],errorType,[1 1],[20 100],[0 1],[0 -100]);
                    allX(sessionInd,1)=X(1);%multiply by 100 as prop_corr given as fraction, not percentage
                    xvals=0:1:testContrast(end)+10;
                    yvals=1./(1+10.^(-X(1).*(xvals-X(2))));
                    PSE(sessionInd)=X(2);
                    %    X=fminsearch(@fit_weibull,X0,[],testContrasts,dummy,[],errorType,[0 0 0 0],[],[0 0 0 0],[]);
                    %    allX(sessionInd,1)=100*X(1)*X(3)*exp(-(sampleContrast/X(2))^X(1))*sampleContrast^(X(1)-1)*(1/X(2))^X(1);%multiply by 100 as prop_corr given as fraction, not percentage
                    %    xvals=0:1:testContrasts(end)+10;
                    %    yvals=1-X(4)-X(3).*exp(-(xvals./X(2)).^X(1));
                    %    PSE(sessionInd)=X(2).*(-log((0.5-X(4))/X(3))).^(1/X(1))
                    plot([PSE(sessionInd) PSE(sessionInd)],[0 1],'Color',[sessionInd/size(ROCmat,1) 0 sessionInd/size(ROCmat,1)]);
                    plot(xvals,yvals,'Color',[sessionInd/size(ROCmat,1) 0 sessionInd/size(ROCmat,1)]);
                    if samePlot==1
                        subplot(1,2,2);
                    else
                        figure(figure2)
                        subplot(ceil(size(ROCmat,1)/6),6,sessionInd);
                    end
                    plot(log(testContrast),dummy,'Color',[sessionInd/size(ROCmat,1) 0 sessionInd/size(ROCmat,1)],'Marker','o','MarkerFaceColor',[sessionInd/size(ROCmat,1) 0 sessionInd/size(ROCmat,1)],'LineStyle','none');hold on
                    X0=[0.5 log(sampleContrast)];
                    X2=fminsearch(@fit_logistic,X0,[],log(testContrast),dummy,[],errorType,[0 1],[0 5],[1 1],[-0.2 -2]);
                    allX2(sessionInd,1)=X2(1);%multiply by 100 as prop_corr given as fraction, not percentage
                    xvals=2:0.001:log(testContrast(end))+0.1;
                    yvals=1./(1+10.^(-X2(1).*(xvals-X2(2))));
                    PSE2(sessionInd)=X2(2);
                    plot([PSE2(sessionInd) PSE2(sessionInd)],[0 1],'Color',[sessionInd/size(ROCmat,1) 0 sessionInd/size(ROCmat,1)]);
                    plot(xvals,yvals,'Color',[sessionInd/size(ROCmat,1) 0 sessionInd/size(ROCmat,1)]);
                    
                    %check slope when midpoint shifts
                    X0=[0.2 sampleContrast];
                    X3=fminsearch(@fit_logistic,X0,[],testContrast-10,dummy,[],errorType,[1 1],[20 100],[0 1],[0 -100]);
                    allX3(sessionInd,1)=X3(1);%multiply by 100 as prop_corr given as fraction, not percentage
                    PSE3(sessionInd)=X3(2);
                    X4=fminsearch(@fit_logistic,X0,[],testContrast+30,dummy,[],errorType,[1 1],[20 100],[0 1],[0 -100]);
                    allX4(sessionInd,1)=X4(1);%multiply by 100 as prop_corr given as fraction, not percentage
                    PSE4(sessionInd)=X4(2);
                    if shiftPSEContrast==1
                        %check what happens to slope for log(contrast) when
                        %midpoint shifts for ordinary contrast:
                        figCompare=figure('Color',[1,1,1],'Units','Normalized','Position',[0, 0, 1, 1]);
                        %plotted against contrast
                        subplot(2,2,1);
                        %weibull fitting:
                        title('weibull');
                        X0=[2 sampleContrast 0.6 0.1];
                        if mean(dummy(1:3))>mean(dummy(end-2:end))
                            X0(1)=-X0(1);
                        end
                        Xw=fminsearch(@fit_weibull,X0,[],testContrast,dummy,[],errorType,[0 0 0 0],[20 100],[0 0 0 0],[0 -1]);
                        allXw(sessionInd,1)=Xw(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSEw(sessionInd)=Xw(2);
%                         strl=['slope: ',num2str(round(100*Xw(1))/100)];
                        strl=[num2str(round(100*Xw(1))/100)];
                        text(Xw(2),dummy(end)+0.05,strl,'FontSize',6);hold on
                        xvals=2:0.001:testContrast(end)+0.1;
                        yvals=1-Xw(4)-Xw(3).*exp(-(xvals./Xw(2)).^Xw(1));
                        plot([PSEw(sessionInd) PSEw(sessionInd)],[0 1],'Color','k');
                        plot(xvals,yvals,'Color','k');
                        plot(testContrast,dummy,'Color','k','Marker','o','MarkerFaceColor','k','LineStyle','none');hold on
                        X5=fminsearch(@fit_weibull,X0,[],testContrast-10,dummy,[],errorType,[0 0 0 0],[20 100],[0 0 0 0],[0 -100]);
                        allX5(sessionInd,1)=X5(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSE5(sessionInd)=X5(2);
%                         str5=['slope: ',num2str(round(100*X5(1))/100)];
                        str5=[num2str(round(100*X5(1))/100)];
                        text(X5(2),dummy(end)+0.05,str5,'FontSize',6);hold on
                        xvals=2:0.001:testContrast(end)+0.1;
                        yvals=1-X5(4)-X5(3).*exp(-(xvals./X5(2)).^X5(1));
                        plot([PSE5(sessionInd) PSE5(sessionInd)],[0 1],'Color','b');
                        plot(xvals,yvals,'Color','b');
                        plot(testContrast-10,dummy,'Color','b','Marker','o','MarkerFaceColor','b','LineStyle','none');hold on
                        X6=fminsearch(@fit_weibull,X0,[],testContrast+30,dummy,[],errorType,[0 1 0 0],[20 100],[1 0 0 0],[0 0]);
                        allX6(sessionInd,1)=X6(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSE6(sessionInd)=X6(2);
%                         str6=['slope: ',num2str(round(100*X6(1))/100)];
                        str6=[num2str(round(100*X6(1))/100)];
                        text(X6(2),dummy(end)+0.05,str6,'FontSize',6);hold on
                        xvals=2:0.001:testContrast(end)+30;
                        yvals=1-X6(4)-X6(3).*exp(-(xvals./X6(2)).^X6(1));
                        plot([PSE6(sessionInd) PSE6(sessionInd)],[0 1],'Color','r');
                        plot(xvals,yvals,'Color','r');
                        plot(testContrast+30,dummy,'Color','r','Marker','o','MarkerFaceColor','r','LineStyle','none');hold on
                        %logistic fitting:
                        subplot(2,2,2);
                        title('logistic');
                        X0=[0.2 sampleContrast];
                        if mean(dummy(1:3))>mean(dummy(end-2:end))
                            X0(1)=-X0(1);
                        end
                        Xl=fminsearch(@fit_logistic,X0,[],testContrast,dummy,[],errorType,[0 0],[20 100],[0 1],[0 -100]);
                        allXl(sessionInd,1)=Xl(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSEl(sessionInd)=Xl(2);
%                         strl=['slope: ',num2str(round(100*Xl(1))/100)];
                        strl=[num2str(round(100*Xl(1))/100)];
                        text(Xl(2),dummy(end)+0.05,strl,'FontSize',6);hold on
                        xvals=2:0.001:testContrast(end)+0.1;
                        yvals=1./(1+10.^(-Xl(1).*(xvals-Xl(2))));
                        plot([PSEl(sessionInd) PSEl(sessionInd)],[0 1],'Color','k');
                        plot(xvals,yvals,'Color','k');
                        plot(testContrast,dummy,'Color','k','Marker','o','MarkerFaceColor','k','LineStyle','none');hold on
                        X5=fminsearch(@fit_logistic,X0,[],testContrast-10,dummy,[],errorType,[0 0],[20 100],[0 0],[0 -100]);
                        allX5(sessionInd,1)=X5(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSE5(sessionInd)=X5(2);
%                         str5=['slope: ',num2str(round(100*X5(1))/100)];
                        str5=[num2str(round(100*X5(1))/100)];
                        text(X5(2),dummy(end)+0.05,str5,'FontSize',6);hold on
                        xvals=2:0.001:testContrast(end)+0.1;
                        yvals=1./(1+10.^(-X5(1).*(xvals-X5(2))));
                        plot([PSE5(sessionInd) PSE5(sessionInd)],[0 1],'Color','b');
                        plot(xvals,yvals,'Color','b');
                        plot(testContrast-10,dummy,'Color','b','Marker','o','MarkerFaceColor','b','LineStyle','none');hold on
                        X6=fminsearch(@fit_logistic,X0,[],testContrast+30,dummy,[],errorType,[0 0],[20 100],[0 0],[0 -100]);
                        allX6(sessionInd,1)=X6(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSE6(sessionInd)=X6(2);
%                         str6=['slope: ',num2str(round(100*X6(1))/100)];
                        str6=[num2str(round(100*X6(1))/100)];
                        text(X6(2),dummy(end)+0.05,str6,'FontSize',6);hold on
                        xvals=2:0.001:testContrast(end)+30;
                        yvals=1./(1+10.^(-X6(1).*(xvals-X6(2))));
                        plot([PSE6(sessionInd) PSE6(sessionInd)],[0 1],'Color','r');
                        plot(xvals,yvals,'Color','r');
                        plot(testContrast+30,dummy,'Color','r','Marker','o','MarkerFaceColor','r','LineStyle','none');hold on
                        %plotted against log(contrast):
                        subplot(2,2,3);
                        %weibull fitting:
                        title('weibull, log(contrast)');
                        X0=[2 log(sampleContrast) 0.6 0.1];
                        if mean(dummy(1:3))>mean(dummy(end-2:end))
                            X0(1)=-X0(1);
                        end
                        Xw=fminsearch(@fit_weibull,X0,[],log(testContrast),dummy,[],errorType,[0 0 0 0],[20 100],[0 0 0 0],[0 -1]);
                        allXw(sessionInd,1)=Xw(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSEw(sessionInd)=Xw(2);
%                         strl=['slope: ',num2str(round(100*Xw(1))/100)];
                        strl=[num2str(round(100*Xw(1))/100)];
                        text(Xw(2),dummy(end)+0.05,strl,'FontSize',6);hold on
                        xvals=2:0.001:log(testContrast(end))+0.5;
                        yvals=1-Xw(4)-Xw(3).*exp(-(xvals./Xw(2)).^Xw(1));
                        plot([PSEw(sessionInd) PSEw(sessionInd)],[0 1],'Color','k');
                        plot(xvals,yvals,'Color','k');
                        plot(log(testContrast),dummy,'Color','k','Marker','o','MarkerFaceColor','k','LineStyle','none');hold on
                        fittedY=1-Xw(4)-Xw(3).*exp(-(log(testContrast)./Xw(2)).^Xw(1));
                        weibullFit(1)=gof_xing(fittedY,dummy);  
%                         X0=[4 log(sampleContrast-10) 0.9 0.2];
%                         if mean(dummy(1:3))>mean(dummy(end-2:end))
%                             X0(1)=-X0(1);
%                         end
                        X5=fminsearch(@fit_weibull,X0,[],log(testContrast-10),dummy,[],errorType,[0 0 0 0],[20 100],[0 0 0 0],[0 -100]);
                        allX5(sessionInd,1)=X5(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSE5(sessionInd)=X5(2);
%                         str5=['slope: ',num2str(round(100*X5(1))/100)];
                        str5=[num2str(round(100*X5(1))/100)];
                        text(X5(2),dummy(end)+0.05,str5,'FontSize',6);hold on
                        xvals=2:0.001:log(testContrast(end))+0.5;
                        yvals=1-X5(4)-X5(3).*exp(-(xvals./X5(2)).^X5(1));
                        plot([PSE5(sessionInd) PSE5(sessionInd)],[0 1],'Color','b');
                        plot(xvals,yvals,'Color','b');
                        plot(log(testContrast-10),dummy,'Color','b','Marker','o','MarkerFaceColor','b','LineStyle','none');hold on
                        fittedY=1-X5(4)-X5(3).*exp(-(log(testContrast-10)./X5(2)).^X5(1));
                        weibullFit(2)=gof_xing(fittedY,dummy);
                        X0=[10 log(sampleContrast) 0.6 0.1];
                        if mean(dummy(1:3))>mean(dummy(end-2:end))
                            X0(1)=-X0(1);
                        end                        
                        X6=fminsearch(@fit_weibull,X0,[],log(testContrast+30),dummy,[],errorType,[0 0 0 0],[20 100],[0 1 0 0],[0 -100]);
                        allX6(sessionInd,1)=X6(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSE6(sessionInd)=X6(2);
%                         str6=['slope: ',num2str(round(100*X6(1))/100)];
                        str6=[num2str(round(100*X6(1))/100)];
                        text(X6(2),dummy(end)+0.05,str6,'FontSize',6);hold on
                        xvals=2:0.001:log(testContrast(end))+0.5;
                        yvals=1-X6(4)-X6(3).*exp(-(xvals./X6(2)).^X6(1));
                        plot([PSE6(sessionInd) PSE6(sessionInd)],[0 1],'Color','r');
                        plot(xvals,yvals,'Color','r');
                        plot(log(testContrast+30),dummy,'Color','r','Marker','o','MarkerFaceColor','r','LineStyle','none');hold on
                        fittedY=1-X6(4)-X6(3).*exp(-(log(testContrast+30)./X6(2)).^X6(1));
                        weibullFit(3)=gof_xing(fittedY,dummy);
                        %logistic fitting:
                        subplot(2,2,4);
                        title('logistic, log(contrast)');
                        X0=[0.2 log(sampleContrast)];
                        if mean(dummy(1:3))>mean(dummy(end-2:end))
                            X0(1)=-X0(1);
                        end
                        Xl=fminsearch(@fit_logistic,X0,[],log(testContrast),dummy,[],errorType,[0 0],[20 100],[0 1],[0 -100]);
                        allXl(sessionInd,1)=Xl(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSEl(sessionInd)=Xl(2);
%                         strl=['slope: ',num2str(round(100*Xl(1))/100)];
                        strl=[num2str(round(100*Xl(1))/100)];
                        text(Xl(2),dummy(end)+0.05,strl,'FontSize',6);hold on
                        xvals=2:0.001:log(testContrast(end))+0.5;
                        yvals=1./(1+10.^(-Xl(1).*(xvals-Xl(2))));
                        plot([PSEl(sessionInd) PSEl(sessionInd)],[0 1],'Color','k');
                        plot(xvals,yvals,'Color','k');
                        plot(log(testContrast),dummy,'Color','k','Marker','o','MarkerFaceColor','k','LineStyle','none');hold on
                        fittedY=1./(1+10.^(-Xl(1).*(log(testContrast)-Xl(2))));
                        logisticFit(1)=gof_xing(fittedY,dummy);  
                        X0=[0.8 log(sampleContrast-10)];
                        if mean(dummy(1:3))>mean(dummy(end-2:end))
                            X0(1)=-X0(1);
                        end       
                        X5=fminsearch(@fit_logistic,X0,[],log(testContrast-10),dummy,[],errorType,[0 0],[20 100],[0 0],[0 -100]);
                        allX5(sessionInd,1)=X5(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSE5(sessionInd)=X5(2);
%                         str5=['slope: ',num2str(round(100*X5(1))/100)];
                        str5=[num2str(round(100*X5(1))/100)];
                        text(X5(2),dummy(end)+0.05,str5,'FontSize',6);hold on
                        xvals=2:0.001:log(testContrast(end))+0.5;
                        yvals=1./(1+10.^(-X5(1).*(xvals-X5(2))));
                        plot([PSE5(sessionInd) PSE5(sessionInd)],[0 1],'Color','b');
                        plot(xvals,yvals,'Color','b');
                        plot(log(testContrast-10),dummy,'Color','b','Marker','o','MarkerFaceColor','b','LineStyle','none');hold on
                        fittedY=1./(1+10.^(-X5(1).*(log(testContrast-10)-X5(2))));
                        logisticFit(2)=gof_xing(fittedY,dummy);    
                        X6=fminsearch(@fit_logistic,X0,[],log(testContrast+30),dummy,[],errorType,[0 0],[20 100],[0 1],[0 -100]);
                        allX6(sessionInd,1)=X6(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSE6(sessionInd)=X6(2);
%                         str6=['slope: ',num2str(round(100*X6(1))/100)];
                        str6=[num2str(round(100*X6(1))/100)];
                        text(X6(2),dummy(end)+0.05,str6,'FontSize',6);hold on
                        xvals=2:0.001:log(testContrast(end))+0.5;
                        yvals=1./(1+10.^(-X6(1).*(xvals-X6(2))));
                        plot([PSE6(sessionInd) PSE6(sessionInd)],[0 1],'Color','r');
                        plot(xvals,yvals,'Color','r');
                        plot(log(testContrast+30),dummy,'Color','r','Marker','o','MarkerFaceColor','r','LineStyle','none');hold on
                        fittedY=1./(1+10.^(-X6(1).*(log(testContrast+30)-X6(2))));
                        logisticFit(3)=gof_xing(fittedY,dummy);                 
                        %save figures of weibull vs logistic, contrast vs
                        %log(contrast) plots:
                        formats=[{'eps'} {'png'}];
                        for j=2:2%1
                                imageFolderName=fullfile('F:','PL','logistic_versus_weibull',animal,area);
                            if ~exist(imageFolderName,'dir')
                                mkdir(imageFolderName);
                            end
                            imageName=[num2str(channels(chInd)),'_',num2str(sessListSorted(sessionInd)),'_weibull_vs_logistic_C_vs_logC'];
                            imagePathName=fullfile(imageFolderName,imageName);
                            %                 export_fig(imageName,formats{j});
                            printtext=sprintf('print -d%s %s',formats{j},imagePathName);
                            eval(printtext)
                        end
                        close(figCompare)
                    end
                    if shiftSlopeContrast==1
                        %check what happens to shape of logistic function for log(contrast) when
                        %slope shifts for ordinary contrast:
                        figure('Color',[1,1,1],'Units','Normalized','Position',[0, 0, 1, 1]);
                        %weibull:
                        subplot(2,2,1);
                        title('weibull');
                        X0=[2 sampleContrast 0.2 0.1];
                        if mean(dummy(1:3))>mean(dummy(end-2:end))
                            X0(1)=-X0(1);
                        end
                        Xw=fminsearch(@fit_weibull,X0,[],testContrast,dummy,[],errorType,[0 0 0 0],[20 100],[0 1 0 0],[0 -1]);
                        allXw(sessionInd,1)=Xw(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSEw(sessionInd)=Xw(2);
%                         strl=['slope: ',num2str(round(100*Xw(1))/100)];
%                         text(Xw(2),dummy(end)+0.05,strl);hold on
                        xvals=2:0.001:testContrast(end)+30;
                        yvals=1-Xw(4)-Xw(3).*exp(-(xvals./Xw(2)).^Xw(1));
%                         plot([PSEw(sessionInd) PSEw(sessionInd)],[0 1],'Color','k');
%                         plot(xvals,yvals,'Color','k');
%                         plot(testContrast,dummy,'Color','k','Marker','o','MarkerFaceColor','none','LineStyle','none');hold on
                        steepTestContrast=testContrast;
                        steepTestContrast(4:11)=[26 28 29 29.5 30.5 31 32 34];
                        X5=fminsearch(@fit_weibull,X0,[],steepTestContrast,dummy,[],errorType,[0 0 0 0],[20 100],[0 0 0 0],[0 -100]);
                        allX5(sessionInd,1)=X5(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSE5(sessionInd)=X5(2);
                        str5=['slope: ',num2str(round(100*X5(1))/100)];
                        text(X5(2),dummy(end)+0.05,str5);hold on
                        xvals=2:0.001:testContrast(end)+30;
                        yvals=1-X5(4)-X5(3).*exp(-(xvals./X5(2)).^X5(1));
                        plot([PSE5(sessionInd) PSE5(sessionInd)],[0 1],'Color','b');
                        plot(xvals,yvals,'Color','b');
                        plot(steepTestContrast,dummy,'Color','b','Marker','x','MarkerFaceColor','none','LineStyle','none');hold on
                        shallowTestContrast=[testContrast(1:7)-9 testContrast(8:14)+10];
%                         shallowTestContrast=[23 25 26 27 33 34 35 37];
                        X6=fminsearch(@fit_weibull,X0,[],shallowTestContrast,dummy,[],errorType,[0 0 0 0],[20 100],[0 1 0 0],[0 -100]);
                        allX6(sessionInd,1)=X6(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSE6(sessionInd)=X6(2);
                        str6=['slope: ',num2str(round(100*X6(1))/100)];
                        text(X6(2),dummy(end)+0.05,str6);hold on
                        xvals=2:0.001:testContrast(end)+30;
                        yvals=1-X6(4)-X6(3).*exp(-(xvals./X6(2)).^X6(1));
                        plot([PSE6(sessionInd) PSE6(sessionInd)],[0 1],'Color','r');
                        plot(xvals,yvals,'Color','r');
                        plot(shallowTestContrast,dummy,'Color','r','Marker','o','MarkerFaceColor','none','LineStyle','none');hold on
                        %logistic
                        subplot(2,2,2);
                        title('logistic');
                        X0=[0.2 sampleContrast];
                        if mean(dummy(1:3))>mean(dummy(end-2:end))
                            X0(1)=-X0(1);
                        end
                        steepTestContrast=testContrast;
                        steepTestContrast(4:11)=[26 28 29 29.5 30.5 31 32 34];
                        X5=fminsearch(@fit_logistic,X0,[],steepTestContrast,dummy,[],'least_square',[0 0],[20 100],[0 0],[0 -100]);
                        allX5(sessionInd,1)=X5(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSE5(sessionInd)=X5(2);
                        str5=['s: ',num2str(round(100*X5(1))/100)];
                        text(X5(2),dummy(end)+0.05,str5);hold on
                        xvals=2:0.001:testContrast(end)+0.1;
                        yvals=1./(1+10.^(-X5(1).*(xvals-X5(2))));
                        plot([PSE5(sessionInd) PSE5(sessionInd)],[0 1],'Color','b');hold on
                        plot(xvals,yvals,'Color','b');
                        plot(steepTestContrast,dummy,'Color','b','Marker','x','MarkerFaceColor','none','LineStyle','none');hold on
                        shallowTestContrast=testContrast;
%                         shallowTestContrast(4:11)=[23 25 26 27 33 34 35 37];
                        shallowTestContrast=[testContrast(1:7)-9 testContrast(8:14)+10];
                        X6=fminsearch(@fit_logistic,X0,[],shallowTestContrast,dummy,[],'least_square',[0 0],[20 100],[0 1],[0 -100]);
                        allX6(sessionInd,1)=X6(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSE6(sessionInd)=X6(2);
                        str6=['s: ',num2str(round(100*X6(1))/100)];
                        text(X6(2),dummy(end)+0.05,str6);hold on
                        xvals=2:0.001:testContrast(end)+0.1;
                        yvals=1./(1+10.^(-X6(1).*(xvals-X6(2))));
                        plot([PSE6(sessionInd) PSE6(sessionInd)],[0 1],'Color','r');
                        plot(xvals,yvals,'Color','r');
                        plot(shallowTestContrast,dummy,'Color','r','Marker','o','MarkerFaceColor','none','LineStyle','none');hold on
                        %log(contrast):
                        %weibull:
                        subplot(2,2,3);
                        title('weibull, log(contrast)');
                        X0=[4 log(sampleContrast) 0.2 0.1];
                        if mean(dummy(1:3))>mean(dummy(end-2:end))
                            X0(1)=-X0(1);
                        end
                        Xw=fminsearch(@fit_weibull,X0,[],log(testContrast),dummy,[],errorType,[0 1 0 0],[20 5],[0 1 0 0],[0 -1]);
                        allXw(sessionInd,1)=Xw(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSEw(sessionInd)=Xw(2);
%                         strl=['slope: ',num2str(round(100*Xw(1))/100)];
%                         text(Xw(2),dummy(end)+0.05,strl);hold on
                        xvals=2:0.001:log(testContrast(end))+0.5;
                        yvals=1-Xw(4)-Xw(3).*exp(-(xvals./Xw(2)).^Xw(1));
%                         plot([PSEw(sessionInd) PSEw(sessionInd)],[0 1],'Color','k');
%                         plot(xvals,yvals,'Color','k');
%                         plot(log(testContrast),dummy,'Color','k','Marker','o','MarkerFaceColor','none','LineStyle','none');hold on
                        steepTestContrast=testContrast;
                        steepTestContrast(4:11)=[26 28 29 29.5 30.5 31 32 34];
                        X5=fminsearch(@fit_weibull,X0,[],log(steepTestContrast),dummy,[],errorType,[0 1 0 0],[20 5],[0 0 0 0],[0 -100]);
                        allX5(sessionInd,1)=X5(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSE5(sessionInd)=X5(2);
                        str5=['s: ',num2str(round(100*X5(1))/100)];
                        text(X5(2),dummy(end)+0.05,str5);hold on
                        xvals=2:0.001:log(testContrast(end))+0.5;
                        yvals=1-X5(4)-X5(3).*exp(-(xvals./X5(2)).^X5(1));
                        plot([PSE5(sessionInd) PSE5(sessionInd)],[0 1],'Color','b');
                        plot(xvals,yvals,'Color','b');
                        plot(log(steepTestContrast),dummy,'Color','b','Marker','x','MarkerFaceColor','none','LineStyle','none');hold on
                        shallowTestContrast=testContrast;
%                         shallowTestContrast(4:11)=[23 25 26 27 33 34 35 37];
                        shallowTestContrast=[testContrast(1:7)-9 testContrast(8:14)+10];
                        X0=[3 log(sampleContrast) 0.2 0.1];
                        if mean(dummy(1:3))>mean(dummy(end-2:end))
                            X0(1)=-X0(1);
                        end
                        X6=fminsearch(@fit_weibull,X0,[],log(shallowTestContrast),dummy,[],errorType,[0 1 0 0],[20 5],[0 0 0 0],[0 0]);
                        allX6(sessionInd,1)=X6(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSE6(sessionInd)=X6(2);
                        str6=['s: ',num2str(round(100*X6(1))/100)];
                        text(X6(2),dummy(end)+0.05,str6);hold on
                        xvals=2:0.001:log(testContrast(end))+0.5;
                        yvals=1-X6(4)-X6(3).*exp(-(xvals./X6(2)).^X6(1));
                        plot([PSE6(sessionInd) PSE6(sessionInd)],[0 1],'Color','r');
                        plot(xvals,yvals,'Color','r');
                        plot(log(shallowTestContrast),dummy,'Color','r','Marker','o','MarkerFaceColor','none','LineStyle','none');hold on
                        %fit logistic
                        subplot(2,2,4);
                        title('logistic, log(contrast)');
                        X0=[0.2 log(sampleContrast)];
                        if mean(dummy(1:3))>mean(dummy(end-2:end))
                            X0(1)=-X0(1);
                        end
                        steepTestContrast=testContrast;
                        steepTestContrast(4:11)=[26 28 29 29.5 30.5 31 32 34];
                        X5=fminsearch(@fit_logistic,X0,[],log(steepTestContrast),dummy,[],'least_square',[0 0],[20 100],[0 0],[0 -100]);
                        allX5(sessionInd,1)=X5(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSE5(sessionInd)=X5(2);
                        str5=['s: ',num2str(round(100*X5(1))/100)];
                        text(X5(2),dummy(end)+0.05,str5);hold on
                        xvals=2:0.001:log(testContrast(end))+0.1;
                        yvals=1./(1+10.^(-X5(1).*(xvals-X5(2))));
                        plot([PSE5(sessionInd) PSE5(sessionInd)],[0 1],'Color','b');hold on
                        plot(xvals,yvals,'Color','b');
                        plot(log(steepTestContrast),dummy,'Color','b','Marker','x','MarkerFaceColor','none','LineStyle','none');hold on
                        shallowTestContrast=testContrast;
%                         shallowTestContrast(4:11)=[23 25 26 27 33 34 35 37];
                        shallowTestContrast=[testContrast(1:7)-9 testContrast(8:14)+10];
                        X6=fminsearch(@fit_logistic,X0,[],log(shallowTestContrast),dummy,[],'least_square',[0 0],[20 100],[0 1],[0 -100]);
                        allX6(sessionInd,1)=X6(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSE6(sessionInd)=X6(2);
                        str6=['s: ',num2str(round(100*X6(1))/100)];
                        text(X6(2),dummy(end)+0.05,str6);hold on
                        xvals=2:0.001:log(testContrast(end))+0.1;
                        yvals=1./(1+10.^(-X6(1).*(xvals-X6(2))));
                        plot([PSE6(sessionInd) PSE6(sessionInd)],[0 1],'Color','r');
                        plot(xvals,yvals,'Color','r');
                        plot(log(shallowTestContrast),dummy,'Color','r','Marker','o','MarkerFaceColor','none','LineStyle','none');hold on
                    end
                    if shiftLogContrastLogistic==1
                        %check what happens to slope for log(contrast) when
                        %midpoint shifts for log(contrast), when fitting logistic:
                        figure('Color',[1,1,1],'Units','Normalized','Position',[0, 0, 1, 1]);
                        subplot(1,2,2);
                        X0=[0.2 log(sampleContrast)];
                        X5=fminsearch(@fit_logistic,X0,[],log(testContrast)-1,dummy,[],errorType,[0 0],[20 100],[0 0],[0 -100]);
                        allX5(sessionInd,1)=X5(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSE5(sessionInd)=X5(2);
                        str1=['slope: ',num2str(round(100*X5(1))/100)];
                        text(X5(2),dummy(end)+0.05,str1);hold on
                        xvals=2:0.001:log(testContrast(end))+0.1;
                        xvals=xvals-1;
                        yvals=1./(1+10.^(-X5(1).*(xvals-X5(2))));
                        plot([PSE5(sessionInd) PSE5(sessionInd)],[0 1],'Color','b');
                        plot(xvals,yvals,'Color','b');hold on
                        plot(log(testContrast)-1,dummy,'Color','b','Marker','o','MarkerFaceColor','b','LineStyle','none');hold on
                        X6=fminsearch(@fit_logistic,X0,[],log(testContrast)+1,dummy,[],errorType,[0 0],[20 100],[0 1],[0 -100]);
                        allX6(sessionInd,1)=X6(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSE6(sessionInd)=X6(2);
                        str2=['slope: ',num2str(round(100*X6(1))/100)];
                        text(X6(2),dummy(end)+0.05,str2);
                        xvals=2:0.001:log(testContrast(end))+0.1;
                        xvals=xvals+1;
                        yvals=1./(1+10.^(-X6(1).*(xvals-X6(2))));
                        plot([PSE6(sessionInd) PSE6(sessionInd)],[0 1],'Color','r');
                        plot(xvals,yvals,'Color','r');
                        plot(log(testContrast)+1,dummy,'Color','r','Marker','o','MarkerFaceColor','r','LineStyle','none');hold on
                        Xl=fminsearch(@fit_logistic,X0,[],log(testContrast),dummy,[],errorType,[0 0],[20 100],[0 1],[0 -100]);
                        allXl(sessionInd,1)=Xl(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSEl(sessionInd)=Xl(2);
                        strl=['slope: ',num2str(round(100*Xl(1))/100)];
                        text(Xl(2),dummy(end)+0.05,strl);
                        xvals=2:0.001:log(testContrast(end))+0.1;
                        yvals=1./(1+10.^(-Xl(1).*(xvals-Xl(2))));
                        plot([PSEl(sessionInd) PSEl(sessionInd)],[0 1],'Color','k');
                        plot(xvals,yvals,'Color','k');
                        plot(log(testContrast),dummy,'Color','k','Marker','o','MarkerFaceColor','k','LineStyle','none');hold on
                    end
                    if shiftLogContrastWeibull==1
                        %check what happens to slope for log(contrast) when
                        %midpoint shifts for log(contrast), when fitting weibull:
                        subplot(1,2,1);
                        X0=[2 log(sampleContrast) 0.2 0.1];
                        X5=fminsearch(@fit_weibull,X0,[],log(testContrast)-1,dummy,[],errorType,[0 0 0 0],[20 100],[0 0 0 0],[0 -100]);
                        allX5(sessionInd,1)=X5(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSE5(sessionInd)=X5(2);
                        str1=['slope: ',num2str(round(100*X5(1))/100)];
                        text(X5(2),dummy(end)+0.05,str1);hold on
                        xvals=2:0.001:log(testContrast(end))+0.1;
                        xvals=xvals-1;
                        yvals=1-X5(4)-X5(3).*exp(-(xvals./X5(2)).^X5(1));
                        plot([PSE5(sessionInd) PSE5(sessionInd)],[0 1],'Color','b');hold on
                        plot(xvals,yvals,'Color','b');hold on
                        plot(log(testContrast)-1,dummy,'Color','b','Marker','o','MarkerFaceColor','b','LineStyle','none');hold on
                        X0=[4 log(sampleContrast) 0.2 0.1];
                        X6=fminsearch(@fit_weibull,X0,[],log(testContrast)+1,dummy,[],errorType,[0 1 0 0],[10 7],[0 0 0 0],[-10 -7]);
                        allX6(sessionInd,1)=X6(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSE6(sessionInd)=X6(2);
                        str2=['slope: ',num2str(round(100*X6(1))/100)];
                        text(X6(2),dummy(end)+0.05,str2);
                        xvals=2:0.001:log(testContrast(end))+0.1;
                        xvals=xvals+1;
                        yvals=1-X6(4)-X6(3).*exp(-(xvals./X6(2)).^X6(1));
                        plot([PSE6(sessionInd) PSE6(sessionInd)],[0 1],'Color','r');
                        plot(xvals,yvals,'Color','r');
                        plot(log(testContrast)+1,dummy,'Color','r','Marker','o','MarkerFaceColor','r','LineStyle','none');hold on
                        X0=[4 log(sampleContrast) 0.2 0.1];
                        Xw=fminsearch(@fit_weibull,X0,[],log(testContrast),dummy,[],errorType,[0 0 0 0],[20 100],[0 0 0 0],[0 -100]);
                        allXw(sessionInd,1)=Xw(1);%multiply by 100 as prop_corr given as fraction, not percentage
                        PSEw(sessionInd)=Xw(2);
                        strw=['slope: ',num2str(round(100*Xw(1))/100)];
                        text(Xw(2),dummy(end)+0.05,strw);
                        xvals=2:0.001:log(testContrast(end))+0.1;
                        yvals=1-Xw(4)-Xw(3).*exp(-(xvals./Xw(2)).^Xw(1));
                        plot([PSEw(sessionInd) PSEw(sessionInd)],[0 1],'Color','k');
                        plot(xvals,yvals,'Color','k');hold on
                        plot(log(testContrast),dummy,'Color','k','Marker','o','MarkerFaceColor','k','LineStyle','none');hold on
                    end
                end
                if samePlot==1
                    subplot(1,2,1);
                    title('contrast');
                    subplot(1,2,2);
                    title('log(contrast)');
                end
                
                figure3=figure('Color',[1,1,1],'Units','Normalized','Position',[0.3, 0.3, 0.3, 0.5]);
                subplot(2,2,1);
                plot(1:size(ROCmat,1),allX,'Marker','o','MarkerFaceColor','k','LineStyle','none');
                title('slope, contrast');
                subplot(2,2,2);
                plot(1:size(ROCmat,1),PSE,'Marker','o','MarkerFaceColor','k','LineStyle','none');hold on
                plot([0 size(ROCmat,1)],[sampleContrast sampleContrast],'Color','r');
                title('PSE, contrast');
                subplot(2,2,3);
                plot(1:size(ROCmat,1),allX2,'Marker','o','MarkerFaceColor','k','LineStyle','none');
                title('slope, log(contrast)');
                subplot(2,2,4);
                plot(1:size(ROCmat,1),PSE2,'Marker','o','MarkerFaceColor','k','LineStyle','none');hold on
                plot([0 size(ROCmat,1)],[log(sampleContrast) log(sampleContrast)],'Color','r');
                title('PSE, log(contrast)');                
                formats=[{'eps'} {'png'}];
                for j=2:2%1
                    imageFolderName=fullfile('F:','PL','logistic_versus_weibull',animal,area);
                    if ~exist(imageFolderName,'dir')
                        mkdir(imageFolderName);
                    end
                    imageName=[num2str(channels(chInd)),'_sl_PSE_weibull_vs_logistic_C_vs_logC'];
                    imagePathName=fullfile(imageFolderName,imageName);
                    %                 export_fig(imageName,formats{j});
                    printtext=sprintf('print -d%s %s',formats{j},imagePathName);
                    eval(printtext)
                end
                
                % figure
                % subplot(1,2,1);
                % plot(1:size(ROCmat,1),allX,'Marker','o','MarkerFaceColor','k','LineStyle','none');
                % subplot(1,2,2);
                % plot(1:size(ROCmat,1),PSE,'Marker','o','MarkerFaceColor','k','LineStyle','none');hold on
                % plot([0 size(ROCmat,1)],[sampleContrast sampleContrast],'Color','r');
                % figure
                % subplot(1,2,1);
                % plot(1:size(ROCmat,1),allX2,'Marker','o','MarkerFaceColor','k','LineStyle','none');
                % subplot(1,2,2);
                % plot(1:size(ROCmat,1),PSE2,'Marker','o','MarkerFaceColor','k','LineStyle','none');hold on
                % plot([0 size(ROCmat,1)],[log(sampleContrast) log(sampleContrast)],'Color','r');
            end
            close all
        end
    end
end