function ROC_and_correlation_xing_data4
%Modified from ROC_and_correlation_xing_data3
animals=[{'blanco'} {'jack'}];
animals=[{'jack'}];
areas=[{'v4_1'} {'v1_1'}];
roc_vec=[];
newroc_vec=[];
corr_vec=[];
allTableStats=[];
allParTableStats=[];
mpallTableStats=[]
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        [sampleContrasts testContrasts]=area_metadata(area);
        channels=main_channels(animal,area);
        sessionNums=main_raw_sessions_final(animal,area,[],0);
        for sampleContrastsInd=1:length(sampleContrasts)
            sampleContrast=sampleContrasts(sampleContrastsInd);
            testContrast=testContrasts(sampleContrastsInd,:);
            %cd F:\PL\sample_test_activity\blanco_v1_1
            % cd F:\PL\sample_test_activity\blanco_v4_1
            % cd F:\PL\sample_test_activity\correct_only\jack_v1_1
            % cd F:\PL\sample_test_activity\jack_v1_1
            %cd F:\PL\sample_test_activity\jack_v4_1
            bins=[0:5:100];
            figure
            
            prob_pars=zeros(1,2);
            roc_pars=zeros(1,2);
            thresholds_roc=zeros(1,4);
            thresholds_prob=zeros(1,4);
            normalize=1;
            equalAct=0;
            for hij=1:length(sessionNums)
                x1=cell(14,1);
                x2=cell(14,1);
                counter=0;
                hij
                for mn=1:1:length(channels)
                    counter=counter+1;
                    fileName=['ch',num2str(channels(mn)),'_',num2str(sessionNums(hij)),'_example_sample_test_act'];
                    loadText=['load F:\PL\sample_test_activity\',animal,'_',area,'\',fileName];
                    eval(loadText)
                    test=subplot(1+ceil(length(channels)/6),6,counter);
                    rho_vec=zeros(1,14);
                    roc_vec=zeros(1,14);
                    prob_vec=zeros(1,14);
                    if normalize==1
                        maxAll=[];%find the highest firing rate across all conditions and trials, across both the sample and test presentation periods
                        for j=1:14;
                            maxAll=[maxAll epoch2{j} epoch4{j}];
                        end
                        maxval=max(maxAll);
                    end
                    for j=1:14;
                        if normalize
                            ep1=((epoch2{j})/maxval)*100;
                            ep2=((epoch4{j})/maxval)*100;
                        else
                            ep1=(epoch2{j});
                            ep2=(epoch4{j});
                        end
                        if isempty(x1{j})
                            x1{j}=[ep1];
                            x2{j}=[ep2];
                        else
                            if length(ep1)==length(x1{j})
                                x1{j}=[cell2mat(x1(j))+ ep1];%sum activity across channels for each trial
                                x2{j}=[cell2mat(x2(j))+ ep2];
                            else
                                fid=fopen('trial_mismatch','a+');
                                fprintf(fid,'channel %d day:%d contrast: %d\n', channels(mn),sessionNums(hij),testContrast(j));
                                fclose(fid);
                            end
                        end
                        x3=cell2mat(x1(j))/(counter);%epoch 2
                        x4=cell2mat(x2(j))/(counter);%epoch 4
                        rho=corr(x3',x4');
                        rho_vec(j)=rho;
                        rho2=corr(ep1',ep2');
                        if rho2>0.7
                            fid=fopen('suspicious_correlation','a+');
                            fprintf(fid,'channel %d day:%d contrast: %d  rho:%5.4f\n', channels(mn),sessionNums(hij),testContrast(j), rho2);
                            fclose(fid);
                        end
                        plot(testContrast(j),(1+rho)/2,'s','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',2);
                        hold on
                        axis square
                        [roc1,vec1,vec2]= sglroc3 (x4,x3);
                        plot(testContrast(j),roc1,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',2);
                        hold on
                        axis square
                        prob_higher=0;
                        prob_lower=0;
                        for kl=1:length(x3)
                            if (x3(kl)<x4(kl))
                                prob_higher=prob_higher+1;
                            elseif (x3(kl)>x4(kl))
                                prob_lower=prob_lower+1;
                            end
                        end
                        prob_higher_fraction=prob_higher/(prob_higher+prob_lower);
                        prob_higher_old=prob_higher/length(x3);
                        if prob_higher_old~=prob_higher_fraction
                            hij
                            mn
                            j
                            equalAct=equalAct+1;
                        end
                        %             prob_vec(1,j)=prob_higher;
                        prob_vec(1,j)=prob_higher_fraction;
                        roc_vec(1,j)=roc1;
                        %             plot(testContrast(j),prob_higher,'o','MarkerFaceColor','b','MarkerSize',2);
                        plot(testContrast(j),prob_higher_fraction,'o','MarkerFaceColor','b','MarkerSize',2);
                        set(test,'YLim',[0 1]);
                        prob_higher=0;
                        prob_lower=0;
                        x3=ep1;
                        x4=ep2;
                        for kl=1:length(x3)
                            if (x3(kl)<x4(kl))
                                prob_higher=prob_higher+1;
                            elseif (x3(kl)>x4(kl))
                                prob_lower=prob_lower+1;
                            end
                        end
                        prob_higher_fraction=prob_higher/(prob_higher+prob_lower);
                        prob_higher_old=prob_higher/length(x3);
                        [roc1,vec1,vec2]= sglroc3 (x4,x3);
                        plot(testContrast(j),roc1,'d','MarkerFaceColor','m','MarkerEdgeColor','m','MarkerSize',2);
                        set(test,'YLim',[0 1]);
                    end
                    meanrho=mean(rho_vec);
                    
                    titletext=sprintf('mean R: %4.3f',meanrho);
                    title(titletext);
                    prl=1;
                    X0=[30 1];
                    if mn==length(channels)
                        X0=[30 1];
                        [pars,fval(1)]=fminsearch(@weibull_zero_one,[X0],[],testContrast,roc_vec);
                        roc_pars(hij,:)=pars;
                        xvals=[0:0.1:100];
                        yvals=1-exp(-((xvals/pars(1)).^pars(2)));
                        yvals2=abs(yvals-0.5);
                        PNE_roc=find(yvals2==min(yvals2));
                        PNE_roc=xvals(PNE_roc);
                        roc_pars(hij,1)=PNE_roc;
                        
                        plot(xvals,yvals,'r');
                        X0=[30 1];
                        [pars,fval(1)]=fminsearch(@weibull_zero_one,[X0],[],testContrast,prob_vec);
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
                        lower_cont=abs(testContrast([7 6 5 4 3 2 1])-30);
                        
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
                        higher_cont=abs(testContrast([8:14])-30);
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
                all_sglroc3_AUROCs(hij,:)=roc_vec;  
                all_AUROCs(hij,:)=prob_vec;             
                %     %%%%%%%%%%%%%%%%%% now calculate noise correlation between channels
                %     %%%%%%%%%%%%%%%%%%% (pairwise)
                %     x3=[];
                %     x4=[];
                %     rho_matr=zeros(length(channels),length(channels));
                %     for mn=1:1:length(channels)-1
                %         x3=[];
                %         x4=[];
                %         for no=mn+1:length(channels)
                %             %searched_channel=sprintf('ch%d_343_example_sample_test_act.mat',channels(mn));
                %             evaltxt=sprintf('data1=load (''%s'')',searched_channel);
                %             eval(evaltxt);
                %             %searched_channel=sprintf('ch%d_343_example_sample_test_act.mat',channels(no));
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
            test=subplot(1+ceil(length(channels)/6),6,counter+1);
            plot([1:1:length(sessionNums)],roc_pars(:,2),'or');
            hold on
            plot([1:1:length(sessionNums)],prob_pars(:,2),'ob');
            title('slope')
            test=subplot(1+ceil(length(channels)/6),6,counter+2);
            plot([1:1:length(sessionNums)],roc_pars(:,1),'or');
            hold on
            plot([1:1:length(sessionNums)],prob_pars(:,1),'ob');
            title('PNE')
            test=subplot(1+ceil(length(channels)/6),6,counter+3);
            plot([1:1:length(sessionNums)],thresholds_roc(:,1),'or');
            hold on
            plot([1:1:length(sessionNums)],thresholds_prob(:,1),'ob');
            title('lower threshold')
            test=subplot(1+ceil(length(channels)/6),6,counter+4);
            plot([1:1:length(sessionNums)],thresholds_roc(:,3),'or');
            hold on
            plot([1:1:length(sessionNums)],thresholds_prob(:,3),'ob');
            title('upper threshold')
            test=subplot(1+ceil(length(channels)/6),6,counter+5);
            plot([1:1:length(sessionNums)],thresholds_roc(:,2),'or');
            hold on
            plot([1:1:length(sessionNums)],thresholds_prob(:,2),'ob');
            title('lower exponent')
            test=subplot(1+ceil(length(channels)/6),6,counter+6);
            plot([1:1:length(sessionNums)],thresholds_roc(:,1),'or');
            hold on
            plot([1:1:length(sessionNums)],thresholds_prob(:,1),'ob');
            title('upper exponent')
            saveText=['save F:\PL\ROC_zero_one\',animal,'\new_vs_old_sglroc3acrosschannels\alex_cumulative_ROCs_',area,'_',num2str(sampleContrast),'.mat prob_pars roc_pars thresholds_roc thresholds_prob all_sglroc3_AUROCs all_AUROCs'];
            eval(saveText)
        end
    end
end
