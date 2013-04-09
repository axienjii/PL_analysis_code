roc_vec=[];

newroc_vec=[];
corr_vec=[];
x1=30+randn(1,1000)*sqrt(30);
x1_1=sort(x1,'ascend');
x2_1=33+randn(1,1000)*sqrt(30);
[x4,IX]=sort(x2_1);
close all
for hi=-1000:1:1000
    if hi>0
        [x4,IX]=sort(x2_1,'ascend');
    
    x2=x4;
    prl=x4(hi:1000);
    prl2=randperm(1001-hi);
    x2(hi:1000)=x2(hi-1+prl2);
    elseif hi<0
        [x4,IX]=sort(x2_1,'descend');
        hi2=abs(hi);
        x2=x4;
        prl=x4(hi2:1000);
        prl2=randperm(1001-hi2);
        x2(hi2:1000)=x2(hi2-1+prl2);
    end
    
    if hi ~=0
    [roc1,vec1,vec2]= sglroc3 (x1_1,x2);
    roc_vec=[roc_vec roc1];
    
    
    counter=0;
    for jk=1:1000
        if x1_1(jk)<x2(jk)
            counter=counter+1;
        end
    end
    roc2=counter/1000;
    newroc_vec=[newroc_vec roc2];
    RHO = corr(x1_1',x2')
    corr_vec=[corr_vec, RHO];
    end
end
plot(roc_vec, corr_vec,'or');
hold on
plot(newroc_vec, corr_vec,'sb');
xlabel('ROC value')
ylabel('single trial correlation')


