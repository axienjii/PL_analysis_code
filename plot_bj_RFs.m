function plot_bj_RFs
%Written by Xing 05/02/2010
%%%%%%%%%%%%%%%%%%%%%%%%%
%to load new RF coordinates and RF sizes, load rfxy.mat, copy and paste
%values from Excel into Matlab workspace, then enter the command
%save rfxy.mat x xSize y ySize SUAx SUAxSize SUAy SUAySize
%%%%%%%%%%%%%%%%%%%%%%%%%

plotUnknown=0;
plotSUA=0;
plotAllRFs=0;
plotFlankers=1;
%can load the following file and enter or adjust values as needed:
animals=[{'blanco'} {'jack'}];
numsSmall=[22 24];numsLarge=[50 44];numsUnknown=55;numsSUA=12;
RFs=[-5 -16;-5 -16];
RF2s=[-3.5 -3;-0.7 -1.3];
sizes=[16 14];
size2s=[3 0.75];
if plotAllRFs==1
    for animalInd=1:length(animals)
        animal=animals{animalInd};
        loadText=['load F:\PL\RFs\',animal,'\rfxy.mat'];
        eval(loadText)
        RF=RFs(animalInd,:);
        RF2=RF2s(animalInd,:);
        size=sizes(animalInd);
        size2=size2s(animalInd);
        numSmall=numsSmall(animalInd);
        numLarge=numsLarge(animalInd);
        xRadius=xSize/2;xSUARadius=SUAxSize/2;
        yRadius=ySize/2;ySUARadius=SUAySize/2;
        figure('Position',[200 10 800 1200]);
        
        stimVertices=[RF(1)+size/2 RF(2)+size/2;RF(1)+size/2 RF(2)-size/2;RF(1)-size/2 RF(2)-size/2;RF(1)-size/2 RF(2)+size/2];
        %patch(stimVertices(:,1),stimVertices(:,2),[0.9 0.9 0.9]);hold on
        %ellipse(size/2,size/2,0,RF(1),RF(2),[0.9 0.9 0.9]);hold on
        rectangle('Position',[RF(1)-size/2,RF(2)-size/2,size,size],'Curvature',[1,1],'FaceColor',[0.9 0.9 0.9],'LineStyle','--');hold on
        rectangle('Position',[RF2(1)-size2/2,RF2(2)-size2/2,size2,size2],'Curvature',[1,1],'FaceColor',[0.9 0.9 0.9],'LineStyle','--');hold on
        
        for i=1:numSmall
            if y(i)<0
                ellipse(xRadius(i),yRadius(i),0,x(i),y(i),'b')
                plot(x(i),y(i),'LineStyle','none','Marker','+','MarkerEdgeColor','b','MarkerFaceColor','g','MarkerSize',6);hold on
            end
        end
        for i=numSmall+1:numLarge
            ellipse(xRadius(i),yRadius(i),0,x(i),y(i),'r');%,'LineStyle',':'
            plot(x(i),y(i),'LineStyle','none','Marker','+','MarkerEdgeColor','r','MarkerFaceColor','g','MarkerSize',6);hold on
        end
        if plotUnknown==1
            for i=numLarge+1:numUnknown
                plot(x(i),y(i),'LineStyle','none','Marker','+','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6);hold on
            end
            for i=numLarge+1:numUnknown-2
                ellipse(xRadius(i),yRadius(i),0,x(i),y(i),'g');
            end
        end
        if plotSUA==1
            for i=1:numSUA
                ellipse(xSUARadius(i),ySUARadius(i),0,SUAx(i),SUAy(i),'m');%,'LineStyle',':'
                plot(SUAx(i),SUAy(i),'LineStyle','none','Marker','+','MarkerEdgeColor','m','MarkerFaceColor','g','MarkerSize',6);hold on
            end
        end
        plot(0,0,'Marker','o','MarkerEdgeColor','k','MarkerSize',5,'MarkerFaceColor','k');
        xGrid=-19:1:9;yGrid=-34:1:4;xVals=-34.9:0.1:4.9;yVals=-19.9:0.1:9.9;
        for i=1:length(xGrid)
            plot(xGrid(i),xVals,'color','k');
        end
        for i=1:length(yGrid)
            plot(yVals,yGrid(i),'color','k');
        end
        grid
        if strcmp(animal,'blanco')
            xlim([-20 5])
            ylim([-30 1])
            xlim([-6 -1.5])
            ylim([-5.5 -1])
        elseif strcmp(animal,'jack')
            xlim([-20 7])
            ylim([-27 1])
            xlim([-1.2 -0.2])
            ylim([-1.8 -0.8])
        end
    end
end

animalLineStyles={'-' '--'};
animalCols={'c' 'm'};
animalCols={[250/255 160/255 160/255] [0.9 0.9 0.9]};
figure('Position',[200 10 800 1200]);
for animalInd=1:length(animals)
    animal=animals{animalInd};
    loadText=['load F:\PL\RFs\',animal,'\rfxy.mat'];
    eval(loadText)
    RF=RFs(animalInd,:);
    RF2=RF2s(animalInd,:);
    size=sizes(animalInd);
    size2=size2s(animalInd);
    rectangle('Position',[RF(1)-size/2,RF(2)-size/2,size,size],'Curvature',[1,1],'LineStyle',animalLineStyles{animalInd},'FaceColor',animalCols{animalInd});hold on
    rectangle('Position',[RF2(1)-size2/2,RF2(2)-size2/2,size2,size2],'Curvature',[1,1],'LineStyle',animalLineStyles{animalInd},'FaceColor',animalCols{animalInd});hold on
    if plotFlankers==1
        rectangle('Position',[RF2(1)-size2/2,RF2(2)-size2*1.5,size2,size2],'Curvature',[1,1],'LineStyle',':');hold on
        rectangle('Position',[RF2(1)-size2/2,RF2(2)+size2*0.5,size2,size2],'Curvature',[1,1],'LineStyle',':');hold on
    end
    plot(0,0,'Marker','o','MarkerEdgeColor','k','MarkerSize',5,'MarkerFaceColor','k');
    xGrid=-19:1:9;yGrid=-34:1:4;xVals=-34.9:0.1:4.9;yVals=-19.9:0.1:9.9;
    for i=1:length(xGrid)
        plot(xGrid(i),xVals,'color','k');
    end
    for i=1:length(yGrid)
        plot(yVals,yGrid(i),'color','k');
    end
    grid
    ylim([-25 2]);
    xlim([-15 5]);
end
pause
