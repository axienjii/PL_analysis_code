roc_vec=[];

newroc_vec=[];
corr_vec=[];
load F:\PL\ch8_343_example_sample_test_act.mat
x3=colormap ('hot');
for j=1:14;
    x3=colormap ('copper');
    x1=epoch2{j};
    x2=epoch4{j};
    test=subplot(1,1,1);
    
    
    rho=corr(x1',x2');
    plot(j,rho,'s','MarkerFaceColor',x3(j*4,:));
    hold on
    axis square
    x3=colormap ('autumn');
    test=subplot(1,1,1);
    [roc1,vec1,vec2]= sglroc3 (x2,x1);
    plot(j,roc1,'o','MarkerFaceColor',x3(j*4,:),'MarkerSize',15);
    hold on 
    axis square
    prob_higher=0;
    for kl=1:length(x1)
        if (x1(kl)<x2(kl))
          prob_higher=prob_higher+1;
        end
    end
    prob_higher=prob_higher/length(x1);
    x3=colormap ('winter');
     plot(j,prob_higher,'o','MarkerFaceColor',x3(j*4,:),'MarkerSize',15);
    
end
