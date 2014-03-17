function Probmat_PL_noise_correlations_xing
%Modified from Alex's function, to check that results can be replicated.

figure
cols=[0 0 1;0 0 1; 1 0 0; 1 0 0];
for kl=1:4
    correlation_matrix2=[];
    if kl==3
        select_data='blanco_v1_1'
        plt=1;
    elseif kl==1
        select_data='blanco_v4_1'
        plt=2;
    elseif kl==4
        select_data='jack_v1_1'
        plt=1;
    elseif kl==2
        select_data='jack_v4_1'
        plt=2
    end
    if strcmp(select_data,'blanco_v1_1')
        cd F:\PL\sample_test_activity\blanco_v1_1
        cont=[5, 10, 15, 20, 22, 25, 28, 32, 35, 40, 45, 50, 60,  90];
        days=[343:1:359];
        days=[343:1:347 355:1:359] %%% here only the first 5 and last five days are taken
        %days=[343 344  358 359] %%% here only the first 2 and last 2 days are taken
        %days=[343:1:350 352:1:359] %%% here  first half and second half days are taken

        chan_nums=[8 9 10 11 15 17 19 21 23 25 26 27 28 29 31 44 45 46 48 61 62 63 64];
    end
    if strcmp(select_data,'blanco_v4_1')
        cd F:\PL\sample_test_activity\blanco_v4_1
        cont=[10, 15, 20, 25, 27, 28, 29, 31, 32, 33, 35, 40, 50, 60];
        days=[307 308 311 313 314 318 320 321 329 330 331:1:341];% blanco V4
        days=[307 308 311 313 314 337:1:341];% blanco V4 %%% here only the first 5 and last five days are taken
        %days=[307 308  340 341];% blanco V4 %%% here only the first 2 and last 2 days are taken
        %days=[307 308 311 313 314 318 320 321 329 330 332:1:341];% blanco V4 %%% %%% here  first half and second half days are taken
        chan_nums=[1 2 3 4 7 12 13 14 18 20 22 24 33 34 36 37 38 40 42 49 50 51 52 53 54 55 57 59 60];%V4 blanco
    end
    if strcmp(select_data,'jack_v1_1')
        cd F:\PL\sample_test_activity\jack_v1_1
        cont=[5, 10, 15, 20, 22, 25, 28, 32, 35, 40, 45, 50, 60,  90];
        days=[51:1:72];% Jack V1 
        days=[51:1:55 68:1:72];% Jack V1 here only the first 5 and last five days are taken
        %days=[51 52  71 72];% Jack V1 here only the first 2 and last 2 days are taken
        %days=[51:1:72];% Jack V1 here  first half and second half days are taken
        chan_nums=[7 9 11 12 13 14 15 16 17 18 19 20 21 22 23 25 26 27 28 29 30 31 32 51 55];% V1 Jack
    end
    if strcmp(select_data,'jack_v4_1')
        cd F:\PL\sample_test_activity\jack_v4_1
        cont=[10, 15, 20, 25, 27, 28, 29, 31, 32, 33, 35, 40, 50, 60];
        days=[24 25 27 28:1:49];% Jack V4
        days=[24 25 27 28 30   44 45 46 48 49];% Jack V4 %%% here only the first 5 and last five days are taken, but with exclusion of days with very high noise correlation 
        days=[24 25 27 28 29   45 46 47 48 49];% Jack V4 %%% here only the first 5 and last five days are taken
        %days=[24 25 48 49];% Jack V4 %%% here only the first 2 and last 2 days are taken, but with exclusion of days with very high noise correlation
        %days=[24 25 27 28:1:31 33:1:49];% Jack V4 here  first half and second half days are taken
        chan_nums=[1 2 3 4 5 6 8 10 24 35 37 39 40 41 49 50 52 53 54 56];% V4 Jack
        
    end
    
    
    
    correlations=nan(length(chan_nums),length(days));
    coefficients=nan(length(chan_nums),2);
    
    x=[-3:0.01:3];
    correlation_matrix=nan(length(chan_nums),(length(chan_nums))*(length(days)));
    count_sig=1;
    for hij=1:length(days)
        sessionRs=[];
        for  mn=1:1:length(chan_nums)
            %counter=counter+1;
            searched_channel1=sprintf('ch%d_%d_example_sample_test_act.mat',chan_nums(mn), days(hij));
            evaltxt=sprintf('data1=load (''%s'');',searched_channel1);
            eval(evaltxt);
            for  no=mn+1:1:length(chan_nums)
                searched_channel2=sprintf('ch%d_%d_example_sample_test_act.mat',chan_nums(no), days(hij));
                evaltxt=sprintf('data2=load (''%s'');',searched_channel2);
                eval(evaltxt);
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%% Noise correlation will be done based on z-scores form
                %%%%%%% sample test difference activities
                %%%%%%% %%%%%%%%%%%%%%%%%%
                
                ep01=[];%these are to calculate correation within days
                ep02=[];%these are to calculate correation within days
                for j=1:14;
                    prl=data1.epoch2{j}-data1.epoch4{j};
                    prl=data1.epoch2{j};
                    prl=(prl-mean(prl))/std(prl);
                    ep01=[ep01 prl];
                    prl=data2.epoch2{j}-data2.epoch4{j};
                    prl=data2.epoch2{j};
                    prl=(prl-mean(prl))/std(prl);
                    ep02=[ep02 prl];
                end
                [r,p]=corrcoef(ep01, ep02);
                %correlation_matrix(no, mn+(hij-1)*(length(chan_nums)+1))=r(2);
                correlation_matrix(mn, no+(hij-1)*(length(chan_nums)+1))=r(2);
                sessionRs=[sessionRs r(2)];
            end
        end
        for no=1:1:length(chan_nums)
            correlation_matrix(no, mn+1+(hij-1)*(length(chan_nums)+1))=1;
        end
%         test=subplot(1,1,1)
%         prl=correlation_matrix(:,1:length(chan_nums))
%         surf(prl','EdgeColor','none');
%         set(test,  'XLim',[1 length(chan_nums)], 'YLim',[1 length(chan_nums)]);
%         GridVisible = 0
%         view(0,90);
%         colorbar
%         x=1;
        correlation_matrix2(hij,:)=sessionRs;
    end
%     %%%%%%% this would plot day resolved nosie correlations   %%%%%%%
%     test=subplot(4,1,kl);
%     CLim1      = get(test,'CLim');
%     surf(correlation_matrix,'EdgeColor','none');
%     GridVisible = 0
%     view(0,90);
%     set(test, 'Ylim',[1 length(chan_nums)],'CLim',[minC maxC]);
%     colorbar
%     prl=find(correlation_matrix<0.9);
%     maxC=max(correlation_matrix(prl));
%     minC=min(correlation_matrix(prl));
%     set(test, 'Xlim',[1 size(correlation_matrix,2)]);
%     
%     prl=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% this plots noise correlation distributions for first vs. second half
%%%% of days
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test=subplot(2,2,kl);
dist1=(correlation_matrix(:, 1:floor(size(correlation_matrix,2)/2)));
prl=find(dist1<0.9);
dist1=dist1(prl);
prl=find(dist1~=0);
dist1=dist1(prl);
dist2=(correlation_matrix(:, floor(size(correlation_matrix,2)/2):end));
prl=find(dist2<0.9);
dist2=dist2(prl);
prl=find(dist2~=0);
dist2=dist2(prl);
[h,p] = ttest2(dist1(:),dist2(:))
titletext=sprintf('%s p: %6.4f', select_data, p);
hold on

text(-0.01,80,titletext, 'color',cols(kl,:));

bins=[-0.95:0.01:0.95];
[X,N]=hist(dist1(:),bins);
[ newbins,newhisto ] = histo_to_bar( bins,X);
plot(newbins,newhisto,'r');
plot([mean(dist1) mean(dist1)],[0 100],'.-r');
hold on
[X,N]=hist(dist2(:),bins);
[ newbins,newhisto ] = histo_to_bar( bins,X);
plot(newbins,newhisto,'b');
set(test,'XLim',[-0.3 0.5]);%, 'Ylim',[0 80]);
plot([mean(dist2) mean(dist2)],[0 100],'.-b');
% 
% test=subplot(2,2,kl);
% early=correlation_matrix2(1:floor(length(days)/2),:);
% hist(early(:),bins);
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','r','facealpha',0.5)
% hold on
% late=correlation_matrix2(1+floor(length(days)/2):end,:);
% hist(late(:),bins);
% h = findobj(gca,'Type','patch');
% set(h,'facealpha',0.5)
% set(test,'XLim',[-0.3 0.5]);%, 'Ylim',[0 80]);

% test=subplot(2,2,kl);
% hist(dist1(:),bins);
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','r','facealpha',0.5)
% hold on
% hist(dist2(:),bins);
% h = findobj(gca,'Type','patch');
% set(h,'facealpha',0.5)
% set(test,'XLim',[-0.3 0.5]);%, 'Ylim',[0 80]);

end
% 
% 
%     