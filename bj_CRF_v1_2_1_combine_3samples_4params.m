function bj_CRF_v1_2_1_combine_3samples_4params
animals=[{'blanco'} {'jack'}];
for animalInd=1:length(animals)
    animal=animals{animalInd};
    allCondSlope=[];
    allCondPNE=[];
    allCondMin=[];
    allCondMax=[];
    sessionArr=[];
    for sampleCond=20:10:40
        loadText=['load F:\PL\CRF\',animal,'\crf_meanchannels\cumulative_CRFs_old_new_v1_2_2_',num2str(sampleCond),'_cutoff10.mat'];
        eval(loadText);
        allCondSlope=[allCondSlope slopeNeuroNew];
        allCondPNE=[allCondPNE PNENew];
        allCondMin=[allCondMin minRateNew];
        allCondMax=[allCondMax maxRateNew];
        sessionArr=[sessionArr 1:length(slopeNeuroNew)];
    end
    spearmanTable=cell(length(animals),3);
    statsObserved=cell(length(animals),4);
    allParams=[allCondSlope;allCondPNE;allCondMin;allCondMax];
    for paramInd=1:4
        [b,stats]=robustfit(sessionArr,allParams(paramInd,:));
        bObserved{animalInd,paramInd}=b(2);
        statsObserved{animalInd,paramInd}=[statsObserved{animalInd,paramInd};b(2) stats.dfe stats.t(2) stats.p(2)];
        allP{animalInd,paramInd}=stats.p(2);
        figure;
        for ind=1:length(allParams(paramInd,:))
            plot(sessionArr(ind),allParams(paramInd,ind),'ko');hold on
        end
        plot([0 length(slopeNeuroNew)],[b(1) b(1)+length(slopeNeuroNew)*b(2)]);
%         [rho p]=corr(chCParray',sessionArr','type','Spearman');
%         if roving==1
%             spearmanTable{animalInd,sampleContrastsInd}=[spearmanTable{animalInd,paramInd};length(chCParray)-2 rho p];
%         elseif roving==0
%             spearmanTable{animalInd,areaInd}=[spearmanTable{animalInd,paramInd};length(chCParray)-2 rho p];
%         end
    end
end