function bj_SE_crf(ch,session,test_epochs,minusSpon,matarray,animal,area,sampleContrast)
%Writes CRF values to mat file:
%session number, time epochs, CRF in spikes/s for each condition durnig
%spontaneous period, during sample, during ISI, and during test.
%
%Set the parameter 'minusSpon' to 1 to subtract spontaneous firing rates:
%spontaneous levels are calculated in units of spikes/s during -150 to 0 ms
%before sample onset. Average firing rates (also in spikes/s) are
%calculated for epochs 4 and 5, and spontaneous rates are subtracted.

durSpon=150;%length of period prior to sample onset from which spontaneous rates are calculated. Can take on a value of up to 512 ms.
minTrials=10;%set value of minumum number of trials for inclusion of session

%             if length(test_epochs)>2
%                 samp_epochs=[0 test_epochs(2)-1024 512];
%             else
%                 samp_epochs=[0 512];
%             end
numTrials=zeros(1,14);
for cond=1:14
    numTrials(cond)=min([length(matarray{cond,1}) length(matarray{cond,2}) length(matarray{cond,3}) length(matarray{cond,4}) length(matarray{cond,5})]);
end
if min(numTrials)>=minTrials
    for cond=1:14
        for epoch=1:size(test_epochs,2)
            if epoch==1
                periods=[-durSpon 0];
            else
                periods=[test_epochs{epoch-1} test_epochs{epoch}(1)];
            end
            for subPeriod=1:length(periods)-1
                    startTime=periods(subPeriod);
                    endTime=periods(subPeriod+1);
                for n=1:numTrials(cond)
                    temp=find(matarray{cond,epoch}{n}<=endTime);
                    temp=find(matarray{cond,epoch}{n}(temp)>startTime);
                    actList{epoch}(n)=length(temp)*1000/(endTime-startTime);
                end
                ave_act{epoch}(subPeriod,cond)=mean(actList{epoch}(:));
            end
        end
    end
    if round(ch)~=ch
        CRFmatName=['CRF_Ch',num2str(round(ch)),'_',num2str(10*(ch-round(ch))),'_',num2str(sampleContrast),'.mat'];
    else
        CRFmatName=['CRF_Ch',num2str(ch),'_',num2str(sampleContrast),'.mat'];
    end
    CRFmatFolder=fullfile('F:','PL','CRF',animal,area);
    if ~exist(CRFmatFolder,'dir')
        mkdir(CRFmatFolder);
    end
    CRFmatPath=fullfile('F:','PL','CRF',animal,area,CRFmatName);
    CRFmatTemp=[{session} {test_epochs} ave_act];
    if ~exist(CRFmatPath,'file')
        CRFmat=CRFmatTemp;
    elseif exist(CRFmatPath,'file')
        loadText=['load ',CRFmatPath,' CRFmat'];
        eval(loadText);
        CRFmat=[CRFmat;CRFmatTemp];
    end
    saveText=['save ',CRFmatPath,' CRFmat'];
    eval(saveText);        
end

