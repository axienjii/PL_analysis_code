function [slopeNeuro,c50,diffc50,minRate,maxRate,chSSE]=nr_fitting(datavals,sampleContrast,testContrast,i,slopeNeuro,chSSE,c50,minRate,maxRate,diffc50,plotDiffC50_30,calculateTangent,startEndTime,animal,area)
%Naka-Rushton counterpart of weibull_fitting function
        xvals=testContrast(1):1:testContrast(end);
        hold on
        fittingFunction='n_r';
        options = optimset('Display','off','MaxFunEvals',10^6,'MaxIter',10^6,'TolFun',1.0E-6,'TolX',1.0E-6);
        if strcmp(animal,'blanco')
            X0=[max(datavals) testContrast(end)-1 0.1 min(datavals)];
            if strcmp(area,'v1_1')
                if i==10
                    X0=[max(datavals) testContrast(end)-1 0.1 min(datavals)+5];
                end
            end
        elseif strcmp(animal,'jack')
            if strcmp(area,'v4_1')
                X0=[max(datavals) 50 1 min(datavals)];%jack v1 X0=[max(datavals) 20 0.1 min(datavals)];
            else
                X0=[max(datavals) 20 0.1 min(datavals)];%jack v1 X0=[max(datavals) 20 0.1 min(datavals)];
            end
        end
        if mean(datavals(1:3))>mean(datavals(end-2:end))%||chNum==13.2||chNum==24||chNum==42
            X0=[max(datavals) 1 -0.1 min(datavals)];%negative slope
        end
        [X]=fminsearch(fittingFunction,X0,options,testContrast,datavals);
        
        fitted_yvals=X(1)*(testContrast.^X(3)./(testContrast.^X(3)+X(2).^X(3)))+X(4);
        residuals=datavals-fitted_yvals;
        sseCRF=sum(residuals.^2);
        chSSE(i,1:2)=[i sseCRF];
        meanFirstThree=mean(datavals(1:3));
        meanMiddleThree=mean(datavals(length(datavals)/2-1):length(datavals)/2+1);
        meanLastThree=mean(datavals(end-2:end));
        %if sse is high and range exceeds 10 spikes/s, or sse is high and
        %period is test, or sse is high and activity seems to be
        %contrast-dependent
        poorSSEcounter=0;
        while sseCRF>5&&poorSSEcounter<30
            if max(datavals)-min(datavals)>10||strcmp(startEndTime,'_1024_to_1536')||meanFirstThree<meanMiddleThree&&meanMiddleThree<meanLastThree||meanFirstThree>meanMiddleThree&&meanMiddleThree>meanLastThree%if the fit seems poor, try a variety of values for the upper and lower limits
                coefEsts2 = zeros(3,4);
                InitVar=1;
                if strcmp(animal,'jack')
                    upperMaxs=[X(1)+1 X(1)+2 X(1)+3];%jack: upperMax=[X(1)+1 X(1)+2 X(1)+3]
                    lowerMins=[X(4)-1];%jack: lowerMin=[X(4)-1]
                    c50Seeds=X(2);%blanco: c50Seed=[5 50]
                end
                if strcmp(animal,'blanco')
                    upperMaxs=[X(1)+1 X(1)+2 X(1)+3];%blanco: upperMax=[X(1)+1 X(1)+2 X(1)+3]
                    lowerMins=[X(4)-1];%blanco: lowerMin=[X(4)-1]
                    c50Seeds=[X(2)];%blanco: c50Seed=[5 50]
                    if strcmp(animal,'blanco')
                        X0=[max(datavals) testContrast(end)-1 0.1 min(datavals)];
                        if strcmp(area,'v1_1')
                            if i==12||i==14
                                c50Seeds=[5 50];%blanco: c50Seed=[5 50] 
                            end
                        end
                    end
                end
                for upperMax=upperMaxs%blanco: upperMax=[X(1)+1 X(1)+2 X(1)+3]
                    for lowerMin=lowerMins%blanco: lowerMin=[X(4)-1]
                        for c50Seed=c50Seeds%blanco: c50Seed=[5 50]
                            X(2)=c50Seed;
                            [coefEsts2(InitVar,:)]=fminsearch(fittingFunction,X,options,testContrast,datavals,[1 1 0 0],[upperMax 100 5 0],[0 1 0 1],[0 0 0 lowerMin]);
                            fitted_yvals=coefEsts2(InitVar,1)*(testContrast.^coefEsts2(InitVar,3)./(testContrast.^coefEsts2(InitVar,3)+coefEsts2(InitVar,2).^coefEsts2(InitVar,3)))+coefEsts2(InitVar,4);
                            residuals=datavals-fitted_yvals;
                            sseCRFtemp(InitVar)=sum(residuals.^2);
                            InitVar=InitVar+1;
                        end
                    end
                end
                [minSSE I]=min(sseCRFtemp);
                if minSSE<sseCRF
                    X=coefEsts2(I,:);
                    chSSE(i,1:2)=[i minSSE];
                    sseCRF=minSSE;
                end
            end
            poorSSEcounter=poorSSEcounter+1;
        end
        if calculateTangent==0
            slopeNeuro(1,i)=X(3);
        elseif calculateTangent==1
            slopeNeuro(1,i)=X(2)^X(3)*sampleContrast^(X(3)-1)/(30^X(3)+X(2)^X(3))^2;
        end
        maxFitted=max(fitted_yvals);
        minFitted=min(fitted_yvals);
        halfFitted=mean([minFitted maxFitted]);
        xvalsFine=testContrast(1):0.01:testContrast(end);
        yvalsFine=X(1)*(xvalsFine.^X(3)./(xvalsFine.^X(3)+X(2).^X(3)))+X(4);
        diffTemp=yvalsFine-halfFitted;%find contrast at which response is (half Rmax)+baseline
        [tempVal columnInd]=min(abs(diffTemp));
        c50(1,i)=xvalsFine(columnInd);
        minRate(1,i)=X(4);
        maxRate(1,i)=X(1);
        %     [X,fval]=fminsearch('weib_sim_min_max',X0,options,testContrast,datavals);
        %     slopeNeuro(1,i)=X(2);
        %     c50(1,i)=X(1).*(-log(0.5)).^(1/X(2));
        %     yvals=max(datavals)-(max(datavals)-min(datavals)).*exp(-((xvals./X(1
        %     )).^X(2)));
        if plotDiffC50_30==1
            diffc50(1,i)=abs(c50(1,i)-30);
        else
            diffc50=[];
        end
        yvals=X(1)*(xvals.^X(3)./(xvals.^X(3)+X(2).^X(3)))+X(4);
        plot(xvals,yvals,'Color','k');