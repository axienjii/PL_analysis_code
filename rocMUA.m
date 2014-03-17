function AUROC=rocMUA(testAct,sampleAct)
%Written by Xing 17/10/13.
%Modified from Mehdi's code

steps = min([testAct sampleAct])-1:1/50:max([testAct sampleAct])+1;
SAMPLE_ROC=zeros(1,length(steps));
TEST_ROC=zeros(1,length(steps));
for i = 1 : length(steps)
    SAMPLE_ROC(i) = (sum(sampleAct>=steps(i)))/ length(sampleAct);
    TEST_ROC(i) = (sum(testAct>=steps(i)))/ length(testAct);
end
AUROC = 0;
for i = 2 : length(TEST_ROC)
    y = (TEST_ROC(i) + TEST_ROC(i-1))/2;
    x = SAMPLE_ROC(i-1) - SAMPLE_ROC(i);
    if y<0||x<0
        pauseHere=1;
    end
    AUROC = AUROC + x*y;
end