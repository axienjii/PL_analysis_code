roc_vec=[];

newroc_vec=[];
corr_vec=[];



%x2_1=x1;
close all
for hi=-1:0.001:1
    m1=30;
    s1=30;
    m2=50;
    s2=50;
    p=hi;
    SIGMA = [s1^2,p*s1*s2;p*s1*s2,s2^2];
    r = mvnrnd([m1 m2],SIGMA,1000);
    minval=min(r);
    minval=min(minval);
    r=r-minval;
    x1=r(:,1);
    x2=r(:,2);
    
    
    
    [roc1,vec1,vec2]= sglroc3 (x1',x2');
    roc_vec=[roc_vec roc1];
    counter=0;
    for jk=1:length(x1)
        if x1(jk)<=x2(jk)
            counter=counter+1;
        end
    end
    roc2=counter/length(x1);
    newroc_vec=[newroc_vec roc2];
    RHO = corr(x1,x2);
    corr_vec=[corr_vec, RHO];
end
test=subplot(1,1,1);
% plot(x1,x2,'or');
% roc2
% counter
% jk
% length(x1)
% hold on
% minval
% plot([-100 150],[-100 150],':k');
% %set(test, 'XLim',[-100 150],'YLim',[-100 150]);
% axis square
% pause
clf
plot(roc_vec, corr_vec,'or');
hold on
plot(newroc_vec, corr_vec,'sb');
xlabel('ROC value')
ylabel('single trial correlation')


