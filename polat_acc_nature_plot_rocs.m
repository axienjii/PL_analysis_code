function[ROCs, thresholds, slopes, p_flank_no_flank_GoF, p_drug, p_no_flank_flank]=polat_acc_nature_plot_rocs(file_of_int0,sgl_trial_akt, n_trials, contrasts, xvals, col, currdir, Cell_o_int, write_ROCs, elim_vec,drug_vec, combine_same_drug)
global histo sgl_trial_akt n_trials ERR num_cond


fighandle3=  figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.2, 0.3, 0.5, 0.5]);
set(fighandle3, 'NumberTitle', 'off', 'Name', 'AT Data panel');
set(fighandle3, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 26.65 20.305]);
max_cont=max(contrasts);
%%%%%%%%%%%%% in case all ROC curves from the different recordings should
%%%%%%%%%%%%% be plotted
p_drug=ones(1,num_cond/2);
if combine_same_drug==0
    num_comp=size(file_of_int0,2);
    p_flank_no_flank_GoF=ones(num_cond/2,num_comp);
    p_no_flank_flank=ones(num_cond/2,num_comp);
    ROCs=zeros(size(file_of_int0,2),4*num_cond);
    thresholds=zeros(num_cond/2,size(file_of_int0,2));
    slopes=zeros(num_cond/2,size(file_of_int0,2));

    new_ntrials=zeros(size(file_of_int0,2),4*num_cond);
    for h=1:size(file_of_int0,2)
        for j=1:4*num_cond
            prl_trials=[1:n_trials(j,h)];
            prl_trials=ismember(prl_trials,elim_vec{h});
            prl_trials=find(prl_trials==0);
            new_ntrials(h,j)=length(prl_trials);

            if j<9
                arr1=sgl_trial_akt(h,j,prl_trials);
                arr1=arr1(:)';
                arr2=sgl_trial_akt(h,1,prl_trials);
                arr2=arr2(:)';
                [roc,vec1,vec2]= sglroc3 (arr1,arr2);
                ROCs(h,j)=roc;
            elseif j<17 & j>8
                arr1=sgl_trial_akt(h,j,prl_trials);
                arr1=arr1(:)';
                arr2=sgl_trial_akt(h,9,prl_trials);
                arr2=arr2(:)';
                [roc,vec1,vec2]= sglroc3 (arr1,arr2);
                ROCs(h,j)=roc;
            elseif j<25 & j>16
                arr1=sgl_trial_akt(h,j,prl_trials);
                arr1=arr1(:)';
                arr2=sgl_trial_akt(h,17,prl_trials);
                arr2=arr2(:)';
                [roc,vec1,vec2]= sglroc3 (arr1,arr2);
                ROCs(h,j)=roc;
            elseif j<33 & j>24
                arr1=sgl_trial_akt(h,j,prl_trials);
                arr1=arr1(:)';
                arr2=sgl_trial_akt(h,25,prl_trials);
                arr2=arr2(:)';
                [roc,vec1,vec2]= sglroc3 (arr1,arr2);
                ROCs(h,j)=roc;
            end
        end
    end
    %%%%%%%%%%%%% in case ROC curves from the identical drug level should be combined

elseif combine_same_drug==1

    num_drugs=[];
    for j=1:length(drug_vec)
        if isempty(find(drug_vec(j)==num_drugs))
            num_drugs=[num_drugs drug_vec(j)];
        end
    end
    num_comp=size(num_drugs,2);
    p_flank_no_flank_GoF=ones(num_cond/2,num_comp);
    p_no_flank_flank=ones(num_cond/2,num_comp);
    ROCs=zeros(num_comp,4*num_cond);
    thresholds=zeros(num_cond/2,num_comp);
    slopes=zeros(num_cond/2,num_comp);
    new_ntrials=zeros(num_comp,4*num_cond);

    akt2=cell(size(num_drugs,2),4*num_cond);
    for g=1:length(num_drugs)
        for h=1:size(file_of_int0,2)
            if drug_vec(h)==num_drugs(g)
                for j=1:4*num_cond
                    prl_trials=[1:n_trials(j,h)];
                    prl_trials=ismember(prl_trials,elim_vec{h});
                    prl_trials=find(prl_trials==0);
                    prl=sgl_trial_akt(h,j,prl_trials);
                    prl=prl(:)';
                    akt2{drug_vec(h)+1,j}=[akt2{drug_vec(h)+1,j} prl];
                    new_ntrials(drug_vec(h)+1,j)=length(akt2{drug_vec(h)+1,j});
                end
            end
        end
    end
    new_ntrials
    for h=1:num_comp
        for j=1:4*num_cond
            if j<9
                arr1=akt2{h,j};
                arr1=arr1(:)';
                arr2=akt2{h,1};
                arr2=arr2(:)';
                 size(arr1)
                 size(arr2)
                [roc,vec1,vec2]= sglroc3 (arr1,arr2);
                ROCs(h,j)=roc;
            elseif j<17 & j>8
                arr1=akt2{h,j};
                arr1=arr1(:)';
                arr2=akt2{h,9};
                arr2=arr2(:)';
                [roc,vec1,vec2]= sglroc3 (arr1,arr2);
                ROCs(h,j)=roc;
            elseif j<25 & j>16
                arr1=akt2{h,j};
                arr1=arr1(:)';
                arr2=akt2{h,17};
                arr2=arr2(:)';
                [roc,vec1,vec2]= sglroc3 (arr1,arr2);
                ROCs(h,j)=roc;
            elseif j<33 & j>24
                arr1=akt2{h,j};
                arr1=arr1(:)';
                arr2=akt2{h,25};
                arr2=arr2(:)';
                [roc,vec1,vec2]= sglroc3 (arr1,arr2);
                ROCs(h,j)=roc;
            end
        end
    end
end

p_flank_no_flank_GoF=ones(num_cond/2,num_comp);
p_flank_no_flank=ones(num_cond/2,num_comp);

%%%%%% ROC curves from the No Flanker conditions
test=subplot(num_cond/2,2,1);
hold on;
maxval=0;
combined_data=[];
combined_ntrials=[];
combined_contrasts=[];
for h=1:num_comp
    ma=ROCs(h,1:8);
    se=ones(1,8);
    n_s=new_ntrials(h,1:8);
    combined_data=[combined_data ma];
    combined_ntrials=[combined_ntrials n_s];
    combined_contrasts=[combined_contrasts contrasts];
    curr_err=plot(contrasts, ma,'o','color',col(h,:),'Markersize',2);
    X0=[25 2];
    options=optimset('MaxFunEvals',10^100);
    options=optimset('MaxIter',10^100);
    X= fminsearch('weibull',X0,options,contrasts,ma);
    err=ERR;
    ma2=1-0.5*exp(-(contrasts./X(1)).^X(2))
    %%% testing real data vs. the fit, i.e. determining goodness of fit
    [residual1, xmle, chisq]=Max_like_estimates_weibull(n_s,ma,ma2);
    p_flank_no_flank_GoF(1,h)=1-chi2cdf(residual1,length(contrasts)-3);
    if X(1)>100
        X(1)=100;
    end
    thresholds(1,h)=X(1);
    slopes(1,h)=X(2);
    R=1-0.5*exp(-(xvals./X(1)).^X(2));
    plot(xvals,R,':', 'color',col(h,:));
    outtext=sprintf('thr: %.2f sl: %.2f p: %.5f',thresholds(1,h), slopes(1,h), p_flank_no_flank_GoF(1,h));
    text(40,0.5+0.1*(h-1),outtext,'Fontsize', [6],'color',col(h,:));
end
set(test,'XLim',[0 max_cont],'YLim',[0.4 1]);
titletext=sprintf('no flankers');
title(titletext)
prefix1=char(file_of_int0(1));
prefix2=find(prefix1=='.');
prefix1=prefix1(1:prefix2-1);
file_names=sprintf('files: %s, ext: ',prefix1);
for k=1:num_comp
    prefix3=char(file_of_int0(k))
    prefix3(prefix2+1:end)
    file_names=sprintf('%s %s, ',file_names, prefix3(prefix2+1:end));
end
titletext=sprintf('%s',currdir);
text(2,1.01,titletext,'FontSize',[6]);
titletext=sprintf('%s  cell: %d',file_names, Cell_o_int);
text(2,0.98,titletext,'FontSize',[6]);
title('no flankers','Fontsize',[8]);

ylabel('ROC value','FontSize',8)
xlabel('contrast [%]','FontSize',8)
prl=gca;
set(prl,'FontSize',[8]);
%%%%%%%%%%%%%% if data from same drug conditions are combined, then test
%%%%%%%%%%%%%% for drug effect on fits
if combine_same_drug==1
    X0=[25 2];
    options=optimset('MaxFunEvals',10^100);
    options=optimset('MaxIter',10^100);
    X3= fminsearch('weibull',X0,options,combined_contrasts,combined_data);
    err=ERR;
    ma_xpect2=1-0.5*exp(-(combined_contrasts./X3(1)).^X3(2))
    [residual2, xmle, chisq]=Max_like_estimates_weibull(combined_ntrials,combined_data,ma_xpect2);
    abs(residual2-residual1)
    p_drug(1,1)=1-chi2cdf(abs(residual2-residual1),1);
    outtext=sprintf('p drug: %.5f', p_drug(1,1));
    text(40,0.5+0.1*((h+1)-1),outtext,'Fontsize', [6],'color',col(2,:));
end




%%%%%% ROC curves from the Flanker with 1*times wavelength separation conditions

test=subplot(num_cond/2,2,3);
hold on;
maxval=0
combined_data=[];
combined_ntrials=[];
combined_contrasts=[];

for h=1:num_comp
    ma=ROCs(h,9:16);
    se=ones(1,16);
    n_s=new_ntrials(h,9:16);
    combined_data=[combined_data ma];
    combined_ntrials=[combined_ntrials n_s];
    combined_contrasts=[combined_contrasts contrasts];

    curr_err=plot(contrasts, ma,'o','color',col(h,:),'Markersize',2);
    X0=[25 2];
    options=optimset('MaxFunEvals',10^100);
    options=optimset('MaxIter',10^100);
    X= fminsearch('weibull',X0,options,contrasts,ma);
    err=ERR;
    ma2=1-0.5*exp(-(contrasts./X(1)).^X(2));
    %%% testing real data vs. the fit, i.e. determining goodness of fit
    [residual1, xmle, chisq]=Max_like_estimates_weibull(n_s,ma,ma2);
    p_flank_no_flank_GoF(2,h)=1-chi2cdf(residual1,length(contrasts)-3);
    if X(1)>100
        X(1)=100;
    end

    thresholds(2,h)=X(1);
    slopes(2,h)=X(2);
    R=1-0.5*exp(-(xvals./X(1)).^X(2))
    plot(xvals,R,':', 'color',col(h,:));
    outtext=sprintf('thres: %.2f slope: %.2f p: %.5f',thresholds(2,h), slopes(2,h), p_flank_no_flank_GoF(2,h));
    text(40,0.5+0.1*(h-1),outtext,'Fontsize', [6],'color',col(h,:));

end
set(test,'XLim',[0 max_cont],'YLim',[0.4 1]);
titletext=sprintf('flankers, wavelength: 1');
title(titletext)
ylabel('ROC value','FontSize',8)
xlabel('contrast [%]','FontSize',8)
prl=gca;
set(prl,'FontSize',[8]);
%%%%%%%%%%%%%% if data from same drug conditions are combined, then test
%%%%%%%%%%%%%% for drug effect on fits

if combine_same_drug==1
    X0=[25 2];
    options=optimset('MaxFunEvals',10^100);
    options=optimset('MaxIter',10^100);
    X3= fminsearch('weibull',X0,options,combined_contrasts,combined_data);
    err=ERR;
    ma_xpect2=1-0.5*exp(-(combined_contrasts./X3(1)).^X3(2))
    [residual2, xmle, chisq]=Max_like_estimates_weibull(combined_ntrials,combined_data,ma_xpect2);
    abs(residual2-residual1)
    p_drug(1,2)=1-chi2cdf(abs(residual2-residual1),1);
    outtext=sprintf('p drug: %.5f', p_drug(1,2));
    text(40,0.5+0.1*((h+1)-1),outtext,'Fontsize', [6],'color',col(2,:));
end
if num_cond==8
    %%%%%% ROC curves from the Flanker with 2*times wavelength separation conditions
    test=subplot(num_cond/2,2,5);
    hold on;
    maxval=0
    combined_data=[];
    combined_ntrials=[];
    combined_contrasts=[];

    for h=1:num_comp
        ma=ROCs(h,25:32);
        n_s=new_ntrials(h,25:32);
        combined_data=[combined_data ma];
        combined_ntrials=[combined_ntrials n_s];
        combined_contrasts=[combined_contrasts contrasts];

        curr_err=plot(contrasts, ma,'o','color',col(h,:),'Markersize',2);
        X0=[25 2];
        options=optimset('MaxFunEvals',10^100);
        options=optimset('MaxIter',10^100);
        X= fminsearch('weibull',X0,options,contrasts,ma);
        err=ERR;
        ma2=1-0.5*exp(-(contrasts./X(1)).^X(2));
        %%% testing real data vs. the fit, i.e. determining goodness of fit
        [residual1, xmle, chisq]=Max_like_estimates_weibull(n_s,ma,ma2);
        p_flank_no_flank_GoF(3,h)=1-chi2cdf(residual1,length(contrasts)-3);
        if X(1)>100
            X(1)=100;
        end

        thresholds(3,h)=X(1);
        slopes(3,h)=X(2);
        R=1-0.5*exp(-(xvals./X(1)).^X(2))
        plot(xvals,R,':', 'color',col(h,:));
        outtext=sprintf('thres: %.2f slope: %.2f p: %.5f',thresholds(3,h), slopes(3,h), p_flank_no_flank_GoF(3,h));
        text(40,0.5+0.1*(h-1),outtext,'Fontsize', [6],'color',col(h,:));

    end
    set(test,'XLim',[0 max_cont],'YLim',[0.4 1]);
    titletext=sprintf('flankers, wavelength: 2');
    title(titletext)
    ylabel('ROC value','FontSize',8)
    xlabel('contrast [%]','FontSize',8)
    prl=gca;
    set(prl,'FontSize',[8]);
    %%%%%%%%%%%%%% if data from same drug conditions are combined, then test
    %%%%%%%%%%%%%% for drug effect on fits
    if combine_same_drug==1
        X0=[25 2];
        options=optimset('MaxFunEvals',10^100);
        options=optimset('MaxIter',10^100);
        X3= fminsearch('weibull',X0,options,combined_contrasts,combined_data);
        err=ERR;
        ma_xpect2=1-0.5*exp(-(combined_contrasts./X3(1)).^X3(2))
        [residual2, xmle, chisq]=Max_like_estimates_weibull(combined_ntrials,combined_data,ma_xpect2);
        abs(residual2-residual1)

        p_drug(1,3)=1-chi2cdf(abs(residual2-residual1),1);
        outtext=sprintf('p drug: %.5f', p_drug(1,3));
        text(40,0.5+0.1*((h+1)-1),outtext,'Fontsize', [6],'color',col(2,:));
    end


    %%%%%% ROC curves from the Flanker with 3*times wavelength separation conditions
    test=subplot(num_cond/2,2,7);
    hold on;
    maxval=0
    combined_data=[];
    combined_ntrials=[];
    combined_contrasts=[];

    for h=1:num_comp
        ma=ROCs(h,17:24);
        n_s=new_ntrials(h,17:24);
        combined_data=[combined_data ma];
        combined_ntrials=[combined_ntrials n_s];
        combined_contrasts=[combined_contrasts contrasts];

        curr_err=plot(contrasts, ma,'o','color',col(h,:),'Markersize',2);
        X0=[25 2];
        options=optimset('MaxFunEvals',10^100);
        options=optimset('MaxIter',10^100);
        X= fminsearch('weibull',X0,options,contrasts,ma);
        err=ERR;
        ma2=1-0.5*exp(-(contrasts./X(1)).^X(2));
        %%% testing real data vs. the fit, i.e. determining goodness of fit
        [residual1, xmle, chisq]=Max_like_estimates_weibull(n_s,ma,ma2);
        p_flank_no_flank_GoF(4,h)=1-chi2cdf(residual1,length(contrasts)-3);
        if X(1)>100
            X(1)=100;
        end

        thresholds(4,h)=X(1);
        slopes(4,h)=X(2);
        R=1-0.5*exp(-(xvals./X(1)).^X(2))
        plot(xvals,R,':', 'color',col(h,:));
        outtext=sprintf('thres: %.2f slope: %.2f p: %.5f',thresholds(4,h), slopes(4,h), p_flank_no_flank_GoF(4,h));
        text(40,0.5+0.1*(h-1),outtext,'Fontsize', [6],'color',col(h,:));

    end
    set(test,'XLim',[0 max_cont],'YLim',[0.4 1]);
    titletext=sprintf('flankers, wavelength: 3');
    title(titletext)
    ylabel('ROC value','FontSize',8)
    xlabel('contrast [%]','FontSize',8)
    prl=gca;
    set(prl,'FontSize',[8]);
    %%%%%%%%%%%%%% if data from same drug conditions are combined, then test
    %%%%%%%%%%%%%% for drug effect on fits
    if combine_same_drug==1
        X0=[25 2];
        options=optimset('MaxFunEvals',10^100);
        options=optimset('MaxIter',10^100);
        X3= fminsearch('weibull',X0,options,combined_contrasts,combined_data);
        err=ERR;
        ma_xpect2=1-0.5*exp(-(combined_contrasts./X3(1)).^X3(2))
        [residual2, xmle, chisq]=Max_like_estimates_weibull(combined_ntrials,combined_data,ma_xpect2);
        abs(residual2-residual1)

        p_drug(1,4)=1-chi2cdf(abs(residual2-residual1),1);
        outtext=sprintf('p drug: %.5f', p_drug(1,4));
        text(40,0.5+0.1*((h+1)-1),outtext,'Fontsize', [6],'color',col(2,:));
    end
elseif num_cond==6
    %%%%%% ROC curves from the Flanker with 3*times wavelength separation conditions
    test=subplot(num_cond/2,2,5);
    hold on;
    maxval=0
    combined_data=[];
    combined_ntrials=[];
    combined_contrasts=[];

    for h=1:num_comp
        ma=ROCs(h,17:24);
        n_s=new_ntrials(h,17:24);
        combined_data=[combined_data ma];
        combined_ntrials=[combined_ntrials n_s];
        combined_contrasts=[combined_contrasts contrasts];

        curr_err=plot(contrasts, ma,'o','color',col(h,:),'Markersize',2);
        X0=[25 2];
        options=optimset('MaxFunEvals',10^100);
        options=optimset('MaxIter',10^100);
        X= fminsearch('weibull',X0,options,contrasts,ma);
        err=ERR;
        ma2=1-0.5*exp(-(contrasts./X(1)).^X(2));
        %%% testing real data vs. the fit, i.e. determining goodness of fit
        [residual1, xmle, chisq]=Max_like_estimates_weibull(n_s,ma,ma2);
        p_flank_no_flank_GoF(3,h)=1-chi2cdf(residual1,length(contrasts)-3);
        if X(1)>100
            X(1)=100;
        end

        thresholds(3,h)=X(1);
        slopes(3,h)=X(2);
        R=1-0.5*exp(-(xvals./X(1)).^X(2))
        plot(xvals,R,':', 'color',col(h,:));
        outtext=sprintf('thres: %.2f slope: %.2f p: %.5f',thresholds(3,h), slopes(3,h), p_flank_no_flank_GoF(3,h));
        text(40,0.5+0.1*(h-1),outtext,'Fontsize', [6],'color',col(h,:));

    end
    set(test,'XLim',[0 max_cont],'YLim',[0.4 1]);
    titletext=sprintf('flankers, wavelength: 3');
    title(titletext)
    ylabel('ROC value','FontSize',8)
    xlabel('contrast [%]','FontSize',8)
    prl=gca;
    set(prl,'FontSize',[8]);
    %%%%%%%%%%%%%% if data from same drug conditions are combined, then test
    %%%%%%%%%%%%%% for drug effect on fits
    if combine_same_drug==1
        X0=[25 2];
        options=optimset('MaxFunEvals',10^100);
        options=optimset('MaxIter',10^100);
        X3= fminsearch('weibull',X0,options,combined_contrasts,combined_data);
        err=ERR;
        ma_xpect2=1-0.5*exp(-(combined_contrasts./X3(1)).^X3(2))
        [residual2, xmle, chisq]=Max_like_estimates_weibull(combined_ntrials,combined_data,ma_xpect2);
        abs(residual2-residual1)

        p_drug(1,3)=1-chi2cdf(abs(residual2-residual1),1);
        outtext=sprintf('p drug: %.5f', p_drug(1,3));
        text(40,0.5+0.1*((h+1)-1),outtext,'Fontsize', [6],'color',col(2,:));
    end
end
test=subplot(num_cond/2-1,2,2);
hold on
%%%%%%%%%%%%%%%%% This currently compares the thresholds and fits for no
%%%%%%%%%%%%%%%%% flanker vs flanker,
for h=1:num_comp
    ma=ROCs(h,1:8);
    n_s=new_ntrials(h,1:8);
    curr_err=plot(contrasts, ma,'o','color',col(h,:),'Markersize',2);
    X0=[25 2];
    options=optimset('MaxFunEvals',10^100);
    options=optimset('MaxIter',10^100);
    X1= fminsearch('weibull',X0,options,contrasts,ma);
    err=ERR;
    R=1-0.5*exp(-(xvals./X1(1)).^X1(2))
    plot(xvals,R,'-', 'color',col(h,:));
    hold on
    ma2=ROCs(h,17:24);
    n_s2=new_ntrials(h,17:24);
    curr_err=plot(contrasts, ma2,'s','color',col(h,:),'Markersize',4);
    X0=[25 2];
    options=optimset('MaxFunEvals',10^100);
    options=optimset('MaxIter',10^100);
    X2= fminsearch('weibull',X0,options,contrasts,ma2);
    err=ERR;
    R=1-0.5*exp(-(xvals./X2(1)).^X2(2))
    plot(xvals,R,':', 'color',col(h,:));
    %%%%%%%%%%%%%%%%%% now fit to combined values %%%%%%%%%%%%%
    %%%%%%%%%%%%% get residual for single fit %%%%%%%%%%%%%%%%
    ma_xpect=1-0.5*exp(-(contrasts./X(1)).^X(2))
    [residual1, xmle, chisq]=Max_like_estimates_weibull(n_s,ma,ma_xpect);

    %%%% fit to combined data %%%%%
    X0=[25 2];
    options=optimset('MaxFunEvals',10^100);
    options=optimset('MaxIter',10^100);
    X3= fminsearch('weibull',X0,options,[contrasts contrasts],[ma ma2]);
    err=ERR;

    ma_xpect2=1-0.5*exp(-([contrasts contrasts]./X3(1)).^X3(2))
    [residual2, xmle, chisq]=Max_like_estimates_weibull([n_s; n_s2],[ma ma2],ma_xpect2);
    abs(residual2-residual1)
    if num_cond==6
        p_no_flank_flank(2,h)=1-chi2cdf(abs(residual2-residual1),1);
        outtext=sprintf('p: %.5f', p_no_flank_flank(2,h));
    elseif  num_cond==8
        p_no_flank_flank(3,h)=1-chi2cdf(abs(residual2-residual1),1);
        outtext=sprintf('p: %.5f', p_no_flank_flank(3,h));
    end
    text(40,0.5+0.1*(h-1),outtext,'Fontsize', [6],'color',col(h,:));
end
set(test,'XLim',[0 max_cont],'YLim',[0.4 1]);
titletext=sprintf('no flankers vs flankers wavelength 3');
title(titletext)
ylabel('ROC value','FontSize',8)
xlabel('contrast [%]','FontSize',8)
prl=gca;
set(prl,'FontSize',[8]);

if num_cond==8
    test=subplot(num_cond/2-1,2,4);
    hold on
    for h=1:num_comp
        ma=ROCs(h,1:8);
        n_s=new_ntrials(h,1:8);
        curr_err=plot(contrasts, ma,'o','color',col(h,:),'Markersize',2);
        X0=[25 2];
        options=optimset('MaxFunEvals',10^100);
        options=optimset('MaxIter',10^100);
        X1= fminsearch('weibull',X0,options,contrasts,ma);
        err=ERR;
        R=1-0.5*exp(-(xvals./X1(1)).^X1(2))
        plot(xvals,R,'-', 'color',col(h,:));
        hold on
        ma2=ROCs(h,25:32);
        n_s2=new_ntrials(h,25:32);
        curr_err=plot(contrasts, ma2,'s','color',col(h,:),'Markersize',4);
        X0=[25 2];
        options=optimset('MaxFunEvals',10^100);
        options=optimset('MaxIter',10^100);
        X2= fminsearch('weibull',X0,options,contrasts,ma2);
        err=ERR;
        R=1-0.5*exp(-(xvals./X2(1)).^X2(2))
        plot(xvals,R,':', 'color',col(h,:));
        %%%%%%%%%%%%%%%%%% now fit to combined values %%%%%%%%%%%%%
        %%%%%%%%%%%%% get residual for single fit %%%%%%%%%%%%%%%%
        ma_xpect=1-0.5*exp(-(contrasts./X(1)).^X(2))
        [residual1, xmle, chisq]=Max_like_estimates_weibull(n_s,ma,ma_xpect);

        %%%% fit to combined data %%%%%
        X0=[25 2];
        options=optimset('MaxFunEvals',10^100);
        options=optimset('MaxIter',10^100);
        X3= fminsearch('weibull',X0,options,[contrasts contrasts],[ma ma2]);
        err=ERR;
        prl_ma1=1-0.5*exp(-(contrasts./X1(1)).^X1(2))
        R=1-0.5*exp(-(xvals./X3(1)).^X3(2))
        %plot(xvals,R,':', 'color','g');

        ma_xpect2=1-0.5*exp(-([contrasts contrasts]./X3(1)).^X3(2))
        [residual2, xmle, chisq]=Max_like_estimates_weibull([n_s; n_s2],[ma ma2],ma_xpect2);
        p_no_flank_flank(2,h)=1-chi2cdf(abs(residual2-residual1),1)
        outtext=sprintf('p: %.5f', p_no_flank_flank(2,h));
        text(40,0.5+0.1*(h-1),outtext,'Fontsize', [6],'color',col(h,:));
    end
    set(test,'XLim',[0 max_cont],'YLim',[0.4 1]);
    titletext=sprintf('no flankers vs flankers wavelength 2');
    title(titletext)
    ylabel('ROC value','FontSize',8)
    xlabel('contrast [%]','FontSize',8)
    prl=gca;
    set(prl,'FontSize',[8]);
    new_ntrials=new_ntrials;
    test=subplot(num_cond/2-1,2,6);
    hold on
    for h=1:num_comp
        ma=ROCs(h,1:8);
        n_s=new_ntrials(h,1:8);
        curr_err=plot(contrasts, ma,'o','color',col(h,:),'Markersize',2);
        X0=[25 2];
        options=optimset('MaxFunEvals',10^100);
        options=optimset('MaxIter',10^100);
        X1= fminsearch('weibull',X0,options,contrasts,ma);
        err=ERR;
        R=1-0.5*exp(-(xvals./X1(1)).^X1(2))
        plot(xvals,R,'-', 'color',col(h,:));
        hold on
        ma2=ROCs(h,9:16);
        n_s2=new_ntrials(h,9:16);
        curr_err=plot(contrasts, ma2,'s','color',col(h,:),'Markersize',4);
        X0=[25 2];
        options=optimset('MaxFunEvals',10^100);
        options=optimset('MaxIter',10^100);
        X2= fminsearch('weibull',X0,options,contrasts,ma2);
        err=ERR;
        R=1-0.5*exp(-(xvals./X2(1)).^X2(2))
        plot(xvals,R,':', 'color',col(h,:));
        %%%%%%%%%%%%%%%%%% now fit to combined values %%%%%%%%%%%%%
        %%%%%%%%%%%%% get residual for single fit %%%%%%%%%%%%%%%%
        ma_xpect=1-0.5*exp(-(contrasts./X(1)).^X(2))
        [residual1, xmle, chisq]=Max_like_estimates_weibull(n_s,ma,ma_xpect);

        %%%% fit to combined data %%%%%
        X0=[25 2];
        options=optimset('MaxFunEvals',10^100);
        options=optimset('MaxIter',10^100);
        X3= fminsearch('weibull',X0,options,[contrasts contrasts],[ma ma2]);
        err=ERR;
        prl_ma1=1-0.5*exp(-(contrasts./X1(1)).^X1(2))
        R=1-0.5*exp(-(xvals./X3(1)).^X3(2))
        %plot(xvals,R,':', 'color','g');

        ma_xpect2=1-0.5*exp(-([contrasts contrasts]./X3(1)).^X3(2))
        [residual2, xmle, chisq]=Max_like_estimates_weibull([n_s; n_s2],[ma ma2],ma_xpect2);
        p_no_flank_flank(1,h)=1-chi2cdf(abs(residual2-residual1),1)
        outtext=sprintf('p: %.5f', p_no_flank_flank(1,h));
        text(40,0.5+0.1*(h-1),outtext,'Fontsize', [6],'color',col(h,:));
    end
    set(test,'XLim',[0 max_cont],'YLim',[0.4 1]);
    titletext=sprintf('no flankers vs flankers wavelength 1');
    title(titletext)
    ylabel('ROC value','FontSize',8)
    xlabel('contrast [%]','FontSize',8)
    prl=gca;
    set(prl,'FontSize',[8]);
    new_ntrials=new_ntrials;
elseif num_cond==6
    test=subplot(num_cond/2-1,2,4);
    hold on
    for h=1:num_comp
        ma=ROCs(h,1:8);
        n_s=new_ntrials(h,1:8);
        curr_err=plot(contrasts, ma,'o','color',col(h,:),'Markersize',2);
        X0=[25 2];
        options=optimset('MaxFunEvals',10^100);
        options=optimset('MaxIter',10^100);
        X1= fminsearch('weibull',X0,options,contrasts,ma);
        err=ERR;
        R=1-0.5*exp(-(xvals./X1(1)).^X1(2))
        plot(xvals,R,'-', 'color',col(h,:));
        hold on
        ma2=ROCs(h,9:16);
        n_s2=new_ntrials(h,9:16);
        curr_err=plot(contrasts, ma2,'s','color',col(h,:),'Markersize',4);
        X0=[25 2];
        options=optimset('MaxFunEvals',10^100);
        options=optimset('MaxIter',10^100);
        X2= fminsearch('weibull',X0,options,contrasts,ma2);
        err=ERR;
        R=1-0.5*exp(-(xvals./X2(1)).^X2(2))
        plot(xvals,R,':', 'color',col(h,:));
        %%%%%%%%%%%%%%%%%% now fit to combined values %%%%%%%%%%%%%
        %%%%%%%%%%%%% get residual for single fit %%%%%%%%%%%%%%%%
        ma_xpect=1-0.5*exp(-(contrasts./X(1)).^X(2))
        [residual1, xmle, chisq]=Max_like_estimates_weibull(n_s,ma,ma_xpect);

        %%%% fit to combined data %%%%%
        X0=[25 2];
        options=optimset('MaxFunEvals',10^100);
        options=optimset('MaxIter',10^100);
        X3= fminsearch('weibull',X0,options,[contrasts contrasts],[ma ma2]);
        err=ERR;
        prl_ma1=1-0.5*exp(-(contrasts./X1(1)).^X1(2))
        R=1-0.5*exp(-(xvals./X3(1)).^X3(2))
        %plot(xvals,R,':', 'color','g');

        ma_xpect2=1-0.5*exp(-([contrasts contrasts]./X3(1)).^X3(2))
        [residual2, xmle, chisq]=Max_like_estimates_weibull([n_s; n_s2],[ma ma2],ma_xpect2);
        p_no_flank_flank(1,h)=1-chi2cdf(abs(residual2-residual1),1)
        outtext=sprintf('p: %.5f', p_no_flank_flank(1,h));
        text(40,0.5+0.1*(h-1),outtext,'Fontsize', [6],'color',col(h,:));
    end
    set(test,'XLim',[0 max_cont],'YLim',[0.4 1]);
    titletext=sprintf('no flankers vs flankers wavelength 1');
    title(titletext)
    ylabel('ROC value','FontSize',8)
    xlabel('contrast [%]','FontSize',8)
    prl=gca;
    set(prl,'FontSize',[8]);
    new_ntrials=new_ntrials;
end