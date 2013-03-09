function plot_psycho_across_sessions(psychoname,sampleContrast,testContrast,excludeSessions,calculateTangent)
fid=fopen(psychoname,'r');%read file and calculate number of sessions
[A,count1]=fscanf (fid,'%s ', inf);
numsessions=count1/(length(testContrast)+5)
fclose(fid);

count=1;
session2=zeros(1,numsessions);
slopePsycho=zeros(1,numsessions);
xvals=testContrast(1):1:testContrast(end);

% X0=[30 2];
X0=[2 30 0.2 0.1];
fid=fopen(psychoname,'r');
while (count<=numsessions)
    [ID,N]=fscanf (fid,'%d', 1);%and session ID
    session2(1,count)=ID;
    [values,count2]=fscanf(fid,'%f', [1,length(testContrast)+4]);
    VALUES(count,1:length(testContrast)+4)=values;
    count=count+1;
end;

for excludeCount=1:length(excludeSessions)
    ind=(session2==excludeSessions(excludeCount));
    VALUES=VALUES(~ind,:);
    session2=session2(~ind);
end

[sessionSorted2 index]=sort(session2);
VALUESsorted=[];
numsessions=length(session2);
for i=1:numsessions
    VALUESsorted(i,:)=VALUES(index(i),:);
end
VALUES=VALUESsorted;

figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
subplotInd=0;
for i=1:numsessions
    subplotInd=subplotInd+1;
    perfvals=VALUES(i,1:length(testContrast));%psychometric performance for each session
    X=fminsearch(@fit_weibull,X0,[],testContrast,perfvals,[],'least_square',[0 0 0 0],[],[0 0 0 0],[]);
    %         [X,fval]=fminsearch('weib_sim_min_max',X0,options,testContrast,perfvals);
    %         slopePsycho(1,i)=X(2);
    %         PSE=X(1).*(-log(0.5)).^(1/X(2));
    %         yvals=max(perfvals)-(max(perfvals)-min(perfvals))*exp(-((xval
    %         s/X(1)).^X(2)));
    if calculateTangent==0
        slopePsycho(1,i)=X(1);
    elseif calculateTangent==1
        slopePsycho(1,i)=X(1)*X(3)*exp(-(30/X(2))^X(1) )*30^(X(1)-1)*(1/X(2))^X(1);
    end
    PSE=X(2).*(-log((0.5-X(4))/X(3))).^(1/X(1));
    subplot(ceil((numsessions)/5),5,subplotInd);
    plot(testContrast,perfvals,'ok');
    hold on
    yvals=1-X(4)-X(3).*exp(-(xvals./X(2)).^X(1));
    plot(xvals,yvals,'r');
    line(PSE,0:0.01:1,'Color','r');
    line(sampleContrast,0:0.01:1);
    ylim([0,max(perfvals)]);
    xlim([0 max(testContrast)]);
    subplottitle=num2str(sessionSorted2(1,i));
    title(subplottitle);
    if i==1
        ptext=sprintf('%s',psychoname);
        orient landscape
        yLimVals=get(gca,'YLim');
        text('Position',[-10 yLimVals(2)+0.2*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
    end
end
printtext=sprintf('print -dpng %s',psychoname);
eval(printtext);