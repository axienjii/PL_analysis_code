roc_vec=[];

newroc_vec=[];
corr_vec=[];
%cd E:\data\blanco_v1_1\blanco_v1_1
cd E:\data\blanco_v4_1\blanco_v4_1
%cd E:\data\jack_v1_1\jack_v1_1
%cd E:\data\jack_v4_1\jack_v4_1
bins=[0:5:100];
cont=[10, 15, 20, 25, 27, 28, 29, 31, 32, 33, 35, 40, 50, 60]; % V4
%cont=[5, 10, 15, 20, 22, 25, 28, 32, 35, 40, 45, 50, 60,  90]; %V1

%days=[343:1:359]; %blanco V1
days=[307 308 311 313 314 318 320 321 329 330 331:1:341];% blanco V4
%days=[51:1:72];% Jack V1
%days=[24 25 27 28:1:49];% Jack V4


chan=dir;
%chan_nums=[8 9 11 15 17 19 21 23 25 26 27 28 29 31 44 45 46 48 61 62 63 64];%V1 blanco

chan_nums=[1 2 3 4 7 12 13 14 18 20 22 24 33 34 36 37 38 40 42 49 50 51 52 53 54 55 57 59 60];%V4 blanco
%chan_nums=[7 9 11 12 13 14 15 16 17 18 19 20 21 22 23 25 26 27 28 29 30 31 32 51 55];% V1 Jack
%chan_nums=[1 2 3 4 5 6 8 10 24 35 37 39 40 41 49 50 52 53 54 56];% V4 Jack
prob_pars=zeros(1,2);
roc_pars=zeros(1,2);
thresholds_roc=zeros(1,4);
thresholds_prob=zeros(1,4);
normalize=1;
for hij=1:length(days)
    x1=cell(14,1);
    x2=cell(14,1);
    counter=0;
    hij
    for mn=1:1:length(chan_nums)
        counter=counter+1;
        searched_channel=sprintf('ch%d_%d_example_sample_test_act.mat',chan_nums(mn), days(hij));
        evaltxt=sprintf('data1=load (''%s'');',searched_channel);
        eval(evaltxt);
        test=subplot(1+ceil(length(chan_nums)/6),6,counter);
        rho_vec=zeros(1,14);
        roc_vec=zeros(1,14);
        prob_vec=zeros(1,14);
        maxval=1;
        if normalize==1
            for j=1:14;
                max1=max(data1.epoch2{j});
                max2=max(data1.epoch4{j});
                max3=max([max1 max2]);
                if maxval<max3
                    maxval=max3;
                end
            end
        end
        for j=1:14;
            if normalize
                ep1=((data1.epoch2{j})/maxval)*100;
                ep2=((data1.epoch4{j})/maxval)*100;
            else
                ep1=(data1.epoch2{j});
                ep2=(data1.epoch4{j});
            end
            if isempty(x1{j})
                x1{j}=[ep1];
                x2{j}=[ep2];
            else
                if length(ep1)==length(x1{j})
                    x1{j}=[cell2mat(x1(j))+ ep1];
                    x2{j}=[cell2mat(x2(j))+ ep2];
                else
                    fid=fopen('trial_mismatch','a+');
                    fprintf(fid,'channel %d day:%d contrast: %d\n', chan_nums(mn),days(hij),cont(j));
                    fclose(fid);
                end
            end
            x3=cell2mat(x1(j))/(counter);
            x4=cell2mat(x2(j))/(counter);
            rho=corr(x3',x4');
            rho_vec(j)=rho;
            rho2=corr(ep1',ep2');
            if rho2>0.7 
                fid=fopen('suspicious_correlation','a+');
                fprintf(fid,'channel %d day:%d contrast: %d  rho:%5.4f\n', chan_nums(mn),days(hij),cont(j), rho2);
                fclose(fid);
            end
            plot(cont(j),(1+rho)/2,'s','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',2);
            hold on
            axis square
            [roc1,vec1,vec2]= sglroc3 (x3,x4);
            plot(cont(j),roc1,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',2);
            hold on
            axis square
            prob_higher=0;
            for kl=1:length(x3)
                if (x3(kl)<x4(kl))
                    prob_higher=prob_higher+1;
                end
            end
            prob_higher=prob_higher/length(x3);
            prob_vec(1,j)=prob_higher;
            roc_vec(1,j)=roc1;
            plot(cont(j),prob_higher,'o','MarkerFaceColor','b','MarkerSize',2);
            set(test,'YLim',[0 1]);
            prob_higher=0;
            x3=ep1;
            x4=ep2;
            for kl=1:length(x3)
                if (x3(kl)<x4(kl))
                    prob_higher=prob_higher+1;
                end
            end
            prob_higher=prob_higher/length(x3);
            [roc1,vec1,vec2]= sglroc3 (x3,x4);
            plot(cont(j),roc1,'d','MarkerFaceColor','m','MarkerEdgeColor','m','MarkerSize',2);
            set(test,'YLim',[0 1]);
        end
        meanrho=mean(rho_vec);
        
        titletext=sprintf('mean R: %4.3f',meanrho);
        title(titletext);
        prl=1;
        X0=[30 1];
        if mn==length(chan_nums)
            X0=[30 1];
            [pars,fval(1)]=fminsearch(@weibull_zero_one,[X0],[],cont,roc_vec);
            roc_pars(hij,:)=pars;
            xvals=[0:0.1:100];
            yvals=1-exp(-((xvals/pars(1)).^pars(2)));
            yvals2=abs(yvals-0.5);
            PNE_roc=find(yvals2==min(yvals2));
            PNE_roc=xvals(PNE_roc);
            roc_pars(hij,1)=PNE_roc;
            
            plot(xvals,yvals,'r');
            X0=[30 1];
            [pars,fval(1)]=fminsearch(@weibull_zero_one,[X0],[],cont,prob_vec);
            yvals=1-exp(-((xvals/pars(1)).^pars(2)));
            plot(xvals,yvals,'b');
            prob_pars(hij,:)=pars;
            yvals2=abs(yvals-0.5);
            PNE_prob=find(yvals2==min(yvals2));
            PNE_prob=xvals(PNE_prob);
            prob_pars(hij,1)=PNE_prob;
            %%%%%%% now calculate thresholds %%%%%%%%%
            lower_rocs=1-roc_vec([7 6 5 4 3 2 1]);
            lower_prob=1-prob_vec([7 6 5 4 3 2 1]);
            lower_cont=abs(cont([7 6 5 4 3 2 1])-30);
            
            X0=[1 2];
            [pars1,fval(1)]=fminsearch(@weibull,[X0],[],lower_cont,lower_rocs);
            thresholds_roc(hij,1:2)=pars1;
            [pars2,fval(1)]=fminsearch(@weibull,[X0],[],lower_cont,lower_prob);
            thresholds_prob(hij,1:2)=pars2;
%             clf
%             test=subplot(2,1,1);
%             plot(lower_cont,lower_prob,'ob');
%             hold on
%             xvals=[0.01:0.01:70];
%             yvals=1-.5*exp(-((xvals./pars2(1)).^pars2(2)));
%             plot(xvals,yvals,'b');
%             plot(lower_cont,lower_rocs,'or');
%             yvals=1-.5*exp(-((xvals./pars1(1)).^pars1(2)));
%             plot(xvals,yvals,'r');
            
            
            higher_rocs=roc_vec([8:14]);
            higher_prob=prob_vec([8:14]);
            higher_cont=abs(cont([8:14])-30);
            [pars1,fval(1)]=fminsearch(@weibull,[X0],[],higher_cont,higher_rocs);
            thresholds_roc(hij,3:4)=pars1;
            [pars2,fval(1)]=fminsearch(@weibull,[X0],[],higher_cont,higher_prob);
            thresholds_prob(hij,3:4)=pars2;
            
%             test=subplot(2,1,2);
%             plot(higher_cont,higher_prob,'ob');
%             hold on
%             xvals=[0.01:0.01:70];
%             yvals=1-.5*exp(-((xvals./pars2(1)).^pars2(2)));
%             plot(xvals,yvals,'b');
%             plot(higher_cont,higher_rocs,'or');
%             yvals=1-.5*exp(-((xvals./pars1(1)).^pars1(2)));
%             plot(xvals,yvals,'r');
%             pause
        end
    end
    
%     %%%%%%%%%%%%%%%%%% now calculate noise correlation between channels
%     %%%%%%%%%%%%%%%%%%% (pairwise)
%     x3=[];
%     x4=[];
%     rho_matr=zeros(length(chan_nums),length(chan_nums));
%     for mn=1:1:length(chan_nums)-1
%         x3=[];
%         x4=[];
%         for no=mn+1:length(chan_nums)
%             %searched_channel=sprintf('ch%d_343_example_sample_test_act.mat',chan_nums(mn));
%             evaltxt=sprintf('data1=load (''%s'')',searched_channel);
%             eval(evaltxt);
%             %searched_channel=sprintf('ch%d_343_example_sample_test_act.mat',chan_nums(no));
%             evaltxt=sprintf('data2=load (''%s'')',searched_channel);
%             eval(evaltxt);
%             for j=1:14;
%                 ep1=(data1.epoch2{j});
%                 ep2=(data2.epoch2{j});
%                 if isempty(x3)
%                     x3=[ep1];
%                     x4=[ep2];
%                 else
%                     x3=[x3 ep1];
%                     x4=[x4 ep2];
%                 end
%             end
%             rho=corr(x3',x4');
%             rho_matr(mn,no)=rho;
%         end
%     end
%     rho_matr
%     pause
%      close all
end
test=subplot(1+ceil(length(chan_nums)/6),6,counter+1);
plot([1:1:length(days)],roc_pars(:,2),'or');
hold on
plot([1:1:length(days)],prob_pars(:,2),'ob');
title('slope')
test=subplot(1+ceil(length(chan_nums)/6),6,counter+2);
plot([1:1:length(days)],roc_pars(:,1),'or');
hold on
plot([1:1:length(days)],prob_pars(:,1),'ob');
title('PNE')
test=subplot(1+ceil(length(chan_nums)/6),6,counter+3);
plot([1:1:length(days)],thresholds_roc(:,1),'or');
hold on
plot([1:1:length(days)],thresholds_prob(:,1),'ob');
title('lower threshold')
test=subplot(1+ceil(length(chan_nums)/6),6,counter+4);
plot([1:1:length(days)],thresholds_roc(:,3),'or');
hold on
plot([1:1:length(days)],thresholds_prob(:,3),'ob');
title('upper threshold')
test=subplot(1+ceil(length(chan_nums)/6),6,counter+5);
plot([1:1:length(days)],thresholds_roc(:,2),'or');
hold on
plot([1:1:length(days)],thresholds_prob(:,2),'ob');
title('lower exponent')
test=subplot(1+ceil(length(chan_nums)/6),6,counter+6);
plot([1:1:length(days)],thresholds_roc(:,1),'or');
hold on
plot([1:1:length(days)],thresholds_prob(:,1),'ob');
title('upper exponent')

