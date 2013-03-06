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

test_epochs=test_epochs{:}';
samp_epochs=test_epochs-529*2;
%             if length(test_epochs)>2
%                 samp_epochs=[0 test_epochs(2)-1024 512];
%             else
%                 samp_epochs=[0 512];
%             end
numTrials=zeros(1,14);
for cond=1:14
    numTrials(cond)=min([length(matarray{cond,1}) length(matarray{cond,4}) length(matarray{cond,5})]);
end
if min(numTrials)>=minTrials
    ave_test_act=zeros(length(test_epochs)-1,14);
    for cond=1:14
        epoch=1;
        spon_act=zeros(length(samp_epochs)-2,length(matarray{cond,epoch}));%list of average activity values from each trial
        for n=1:numTrials(cond)
            temp=find(matarray{cond,epoch}{n}<=0);
            temp=find(matarray{cond,epoch}{n}(temp)>-1*durSpon);
            spon_act(1,n)=length(temp)*1000/durSpon;
        end
        
        %write activity levels for sample presentation period (epoch 2) and post-sample period (epoch 3) to 1 array
        %write activity levels for test presentation period (epoch 4) and post-test period (epoch 5) to 1 array
        epoch=2;
        samp_act=zeros(length(samp_epochs)-2,length(matarray{cond,epoch}));%list of average activity values from each trial
        for i=1:find(samp_epochs==529)-1%samp_epochs spans 0 to 512 to 912, possibly with subdivisions in between- analyse just 0 to 512 for epoch 2, when calculating samp_act
            for n=1:length(matarray{cond,epoch})
                temp=find(matarray{cond,epoch}{n}<=samp_epochs(i+1));
                temp=find(matarray{cond,epoch}{n}(temp)>samp_epochs(i));
                samp_act(i,n)=length(temp)*1000/((samp_epochs(i+1)-samp_epochs(i)));
            end
        end
        epoch=3;
        isi_act=zeros(length(samp_epochs)-2,length(matarray{cond,epoch}));%list of average activity values from each trial
        for i=find(samp_epochs==529)-1:length(samp_epochs)-2%samp_epochs spans 0 to 512 to 912, possibly with subdivisions in between- analyse just 512 to 912 for epoch 3, when calculating post_samp_act
            for n=1:length(matarray{cond,epoch})
                temp=find(matarray{cond,epoch}{n}<=samp_epochs(i+2));
                temp=find(matarray{cond,epoch}{n}(temp)>samp_epochs(i+1));
                isi_act(i,n)=length(temp)*1000/((samp_epochs(i+2)-samp_epochs(i+1)));
            end
        end        
        epoch=4;
        test_act=zeros(length(test_epochs)-2,length(matarray{cond,epoch}));%list of average activity values from each trial
        for i=1:find(test_epochs==529*3)-1%test_epochs spans 1024 to 1536 to 1936, possibly with subdivisions in between- analyse just 1024 to 1536 for epoch 4, when calculating test_act
            for n=1:numTrials(cond)
                temp=find(matarray{cond,epoch}{n}<=test_epochs(i+1));
                temp=find(matarray{cond,epoch}{n}(temp)>test_epochs(i));
                test_act(i,n)=length(temp)*1000/((test_epochs(i+1)-test_epochs(i)));
            end
        end
%         epoch=5;
%         for i=find(test_epochs==1536):length(test_epochs)-1%test_epochs spans 1024 to 1536 to 1936, possibly with subdivisions in between- analyse just 1536 to 1936 for epoch 5, when calculating post_test_act
%             for n=1:numTrials(cond)
%                 temp=find(matarray{cond,epoch}{n}<=test_epochs(i+1));
%                 temp=find(matarray{cond,epoch}{n}(temp)>test_epochs(i));
%                 test_act(i,n)=length(temp)*1000/((test_epochs(i+1)-test_epochs(i)));
%             end
%         end
        for i=1:length(test_epochs)-2
            ave_spon_act(i,cond)=mean(spon_act(i,:));
            ave_samp_act(i,cond)=mean(samp_act(i,:));
            ave_isi_act(i,cond)=mean(isi_act(i,:));
            ave_test_act(i,cond)=mean(test_act(i,:));
            if minusSpon==1;
                ave_test_act(i,cond)=ave_test_act(i,cond)-mean(spon_act(1,:));
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
    for i=1:length(test_epochs)-2
        CRFmatTemp=[{session} {test_epochs} {ave_spon_act(i,:)} {ave_samp_act(i,:)} {ave_isi_act(i,:)} {ave_test_act(i,:)}];
    end
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

