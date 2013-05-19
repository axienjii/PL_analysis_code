roc_vec=[];

newroc_vec=[];
corr_vec=[];

x1=poissrnd(10,1,1000);
[x1_sorted,IX]=sort(x1,'ascend');
x2=poissrnd(12,1,1000);
[x2_sorted,IX]=sort(x2,'descend');
mean(x2)
var(x2)

%x2_1=x1;
close all
for hi=0:1:1000
    x3=x2_sorted;
    
    xpos=randperm(1000,ceil(hi/2)*2);
    
    x3(xpos(1:length(xpos)/2))=x2_sorted(xpos((length(xpos)/2)+1:end));
    x3(xpos((length(xpos)/2)+1:end))=x2_sorted(xpos((1:length(xpos)/2)));
    
    
    [roc1,vec1,vec2]= sglroc3 (x1_sorted,x3);
    roc_vec=[roc_vec roc1];
    counter=0;
    for jk=1:1000
        if x1_sorted(jk)<=x3(jk)
            counter=counter+1;
        end
    end
    roc2=counter/1000;
    newroc_vec=[newroc_vec roc2];
    RHO = corr(x1_1',x3');
    corr_vec=[corr_vec, RHO];
end
[x2_sorted,IX]=sort(x2,'ascend');
for hi=0:1:1000
    x3=x2_sorted;
    
    xpos=randperm(1000,ceil(hi/2)*2);
    
    x3(xpos(1:length(xpos)/2))=x2_sorted(xpos((length(xpos)/2)+1:end));
    x3(xpos((length(xpos)/2)+1:end))=x2_sorted(xpos((1:length(xpos)/2)));
    
    
    [roc1,vec1,vec2]= sglroc3 (x1_sorted,x3);
    roc_vec=[roc_vec roc1];
    counter=0;
    for jk=1:1000
        if x1_sorted(jk)<=x3(jk)
            counter=counter+1;
        end
    end
    roc2=counter/1000;
    newroc_vec=[newroc_vec roc2];
    RHO = corr(x1_1',x3');
    corr_vec=[corr_vec, RHO];
end
plot(roc_vec, corr_vec,'or');
hold on
plot(newroc_vec, corr_vec,'sb');
xlabel('ROC value')
ylabel('single trial correlation')


