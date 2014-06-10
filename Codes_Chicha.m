%%%%%save lfps raws

clear
animals{1} ='blanco';
animals{2} = 'jack';
areas= ['v4';'v1'];
%areas= ['v4'];

for nm=1:size(animals,2)
    for xc=1:size(areas,1)
        area=areas(xc,:);
        animal=animals{nm};
        if strcmpi(animal,'blanco')
            if strncmp(area,'v4',2)
                channels= [1 2 3 4 7 12 13 14 18 20 22 24 33 34 36 37 38 40 42 49 50 51 52 53 54 55 57 59 60];
            elseif strncmp(area,'v1',2)
                channels=[8 9 10 11 15 17 19 21 23 25 26 27 28 29 31 44 45 46 48 61 62 63 64];
            end
        elseif strcmpi(animal,'jack')
            if strncmp(area,'v4',2)
                channels=[1:5 6 8 10 24 35 37 39 40 41 49 50 52:54 56];
            elseif strncmp(area,'v1',2)
                channels=[7 9 11 12 13 14 15 16 17 18 19 20 21 22 23 25 26 27 28 29 30 31 32 51 55];
            end
        end
        cons = [8:1:14];
        N = 2048;
        sessions = main_raw_sessions_final(animal,area,1,0);
        
        Lfp = zeros(30,length(cons),120,length(sessions),N);
        Sps = zeros(30,length(cons),120,length(sessions),N);
        
        Ltes = zeros(length(cons),120,length(sessions));
        nts = zeros(length(cons),length(sessions));
        
        for nj=1:length(sessions)
            
            session= sessions(nj);
            
            loadText=['load /home/chamanthi/Thiele_data/pl_LFP/',animal,'/',num2str(session),'/SGL_trial_LFP_-512_1536_ch',num2str(channels(6)),'.mat LFP_data event_arr'];
            eval(loadText)
            loadText =['load /home/chamanthi/Thiele_data/pl_corr_art_trials_copy/',animal,'/',num2str(session),'_corrtrialartifact.mat'];
            eval(loadText)
            
            
            for k =1:length(cons)
                
                lfps00 = squeeze(LFP_data(cons(k),1:event_arr(1,cons(k)),2:2049));
                te = squeeze(LFP_data(cons(k),1:event_arr(1,cons(k)),1));
                tes =[];
                %%%eliminate the trials corresponding to the ones removed
                %%%for the spikes
                for kk=1:length(te)
                    if isempty(find(removeTrialsTimestamps(:,1)==te(kk)))==1
                        tes=[tes kk];
                    end
                end
                if length(tes)>120
                   tes = tes(1:120); 
                end
                Ltes(k,1:length(tes),nj) = tes;
                nts(k,nj)= length(tes);
            end
            
            
            ntes = zeros(length(cons),120);
            for bi=1:length(channels)
                channel = channels(bi);
                [nm xc nj bi]
                loadText=['load /home/chamanthi/Thiele_data/pl_LFP/',animal,'/',num2str(session),'/SGL_trial_LFP_-512_1536_ch',num2str(channel),'.mat LFP_data event_arr'];
                eval(loadText)
                loadText=['load /home/chamanthi/Thiele_data/spikeData/',animal,'/',num2str(channel),'_',num2str(session),'_30.mat'];
                eval(loadText)
                
                for k =1:length(cons)
                    
                    lfps00 = squeeze(LFP_data(cons(k),1:event_arr(1,cons(k)),2:2049));
                    ttes = Ltes(k,1:nts(k,nj),nj);
                    
                    epochs=[1 2 3 4];
                    for n=1:length(ttes)
                        
                        st1 =[];
                        for kk=1:length(epochs)
                            st0= 512+round(matarray{cons(k),kk}{n});
                            st1= [st1 st0];
                        end
                        st =intersect(st1,[1:2048]);
                        Sps(bi,k,n,nj,st)=1;
                    end
                    %%%%%eliminate the trials that have saturations or
                    %%%%%strange trends
                    lfps0= lfps00(ttes,:);
                    mst =0;
                    sst = zeros(size(lfps0,1),1);
                    for jk = 1:size(lfps0,1)
                        sst(jk) = std(squeeze(lfps0(jk,:)));
                        mst = mst +sst(jk);
                    end
                    mst = mst/size(lfps0,1);
                    %sttt =std(lfps0,[],2);
                    %mst = mean(sttt);
                    te =find(sst<1.5*mst);
                    
                    if length(te)>120
                        te = te(1:120);
                    end
                    nl = length(te);
                    ntes(k,te) = ntes(k,te)+1;
                    lfps =lfps0(te,:);
                    for n=1:nl
                       s =squeeze(lfps(n,:));
                       Lfp(bi,k,te(n),nj,:)=squeeze(Lfp(bi,k,te(n),nj,:))+ s';
                    end
                end
            end
            
        end
        tl = sprintf('save /home/chamanthi/Thiele_data/datacorrel14/datalfpsraw1%s_%s  Ltes Sps nts Lfp',animal,area);
        eval(tl);
        
    end
end

%%The average is generated during the run of the discriminability analysis, is the variable mlfp


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%discriminability in quartiles with averaged LFP across channels


clear
animals{1} ='blanco';
animals{2} = 'jack';
areas= ['v4';'v1'];
%areas= ['v4'];
nb = 300;
for nm=1:size(animals,2)
    for xc=1:size(areas,1)
        area=areas(xc,:);
        animal=animals{nm};
        if strcmpi(animal,'blanco')
            if strncmp(area,'v4',2)
                channels= [1 2 3 4 7 12 13 14 18 20 22 24 33 34 36 37 38 40 42 49 50 51 52 53 54 55 57 59 60];
            elseif strncmp(area,'v1',2)
                channels=[8 9 10 11 15 17 19 21 23 25 26 27 28 29 31 44 45 46 48 61 62 63 64];
            end
        elseif strcmpi(animal,'jack')
            if strncmp(area,'v4',2)
                channels=[1:5 6 8 10 24 35 37 39 40 41 49 50 52:54 56];
            elseif strncmp(area,'v1',2)
                channels=[7 9 11 12 13 14 15 16 17 18 19 20 21 22 23 25 26 27 28 29 30 31 32 51 55];
            end
        end
        sessions = main_raw_sessions_final(animal,area,1,0);
        
        nch(nm,xc) = length(channels);
        cons = [1:1:7];
        N = 2048;
        fqs = [1 5 8 13 18 30 60];
        fqs2 = [4 8 12 18 23 35 65];
        Fs = 1/(0.001);
        hFs = Fs/2;
        nf =length(fqs);
        nsu = 19;
        t0s = [512 1536];
        wissps = [100 500];
        q = 3;
        th = 20;
        
        tl = sprintf('load /home/chamanthi/Thiele_data/datacorrel14/datalfpsraw2%s_%s',animal,area);
        eval(tl);
        numsessionssel = zeros(length(cons),1);
        
        selsess = zeros(length(cons),length(sessions));
        mauroc = zeros(nf,length(cons),length(t0s),length(wissps),2,2);
        muauroc = zeros(nf,length(cons),length(t0s),length(wissps),nsu,2,2);
 
        mSps = squeeze(sum(Sps(1:length(channels),:,:,:,:),1));
        
        
        for k =1:length(cons)
            [nm xc k]
            
            cx =0;
            for nj=1:length(sessions)
                nu = nts(k,nj);   
                lfps = squeeze(Lfp(1:length(channels),k,1:nu,nj,:));
                slfps = sum(lfps,3);
                mlfp = zeros(nu,N);
                seltr = zeros(nu,1);
                for rt =1:nu
                  fd = find(abs(slfps(:,rt))>0);
                  if length(fd)>0
                    seltr(rt)=1;
                    %%%the averae is done only across the "healthy"
                    %%%channels
                    mlfp(rt,:)= mean(lfps(fd,rt,:),1);
                  end
                end
                nu = sum(seltr);
                num = round(nu/q);
                if mod(num,2)>0
                   num = num+1; 
                end
            
                if num >=th
                    cx = cx+1;
                    selsess(k,nj)=1;
                    ttr = find(seltr==1);
                    sps = squeeze(mSps(k,ttr,nj,:));
                    
                    tyy= zeros(nsu,num,2);
                    for nn = 1:nsu
                        tyy(nn,:,1) = randperm(num);
                        tyy(nn,:,2) = randperm(num);
                    end
                    for j = 1:length(fqs)
                        
                        lband = fqs(j);
                        uband = fqs2(j);
                        [B,A] = butter(3,[lband/hFs uband/hFs]);
                        tpw = zeros(nu,N);
                        for n=1:nu
                            
                            s =squeeze(mlfp(ttr(n),:));
                            t = filtfilt(B,A,s-mean(s));
                            
                            tpw(n,:)=abs(hilbert(t-mean(t)))';
                            
                        end
                        for kk = 1:length(t0s)
                            tt = t0s(kk);
                            %npw = squeeze(sum(tpw(:,1:tt),2));
                            %npw = squeeze(sum(tpw(:,tt-400:tt),2));
                            npw = squeeze(sum(tpw(:,tt:tt+400),2));
                            
                            [bn bv]= sort(npw);
                            sps2 = sps(bv(nu-num+1:nu),:);
                            sps1 = sps(bv(1:num),:);
                            
                            %%lower quartile
                            tt = 513;
                            tt2 = 1537;
                            for k2=1:length(wissps)
                                
                                nsps1l = squeeze(sum(sps1(:,tt:tt+wissps(k2)),2));
                                nsps2l = squeeze(sum(sps1(:,tt2:tt2+wissps(k2)),2));
                                nspsl = [nsps1l nsps2l];
                                
                                nsps1u = squeeze(sum(sps2(:,tt:tt+wissps(k2)),2));
                                nsps2u = squeeze(sum(sps2(:,tt2:tt2+wissps(k2)),2));
                                nspsu = [nsps1u nsps2u];
                                
                                dt = nsps2l-nsps1l;
                                tdt = find(dt>0);
                                tdt2 = find(dt==0);
                                 
                                %mauroc(njj,j,k,kk,k2,1,2) = mauroc(njj,j,k,kk,k2,1,2) +0.5-length(tdt)/num-length(tdt2)/(2*num);
                                mauroc(j,k,kk,k2,1,2) = mauroc(j,k,kk,k2,1,2) +0.5-length(tdt)/(num-length(tdt2));
                                
                                maxr = max([max(nsps1l) max(nsps2l) max(nsps1u) max(nsps2u)]);
                                minr = min([min(nsps1l) min(nsps2l) min(nsps1u) min(nsps2u)]);
                                ib = minr-(maxr-minr)/nb:(maxr-minr)/nb:maxr+(maxr-minr)/nb;

                                ttd = find(abs(squeeze(nspsl(:,1))-squeeze(nspsl(:,2)))>0);
                                lttd = length(ttd);

                                n1 = histc(squeeze(nspsl(ttd,1)),ib)/lttd;
                                n2 = histc(squeeze(nspsl(ttd,2)),ib)/lttd;
                                mk = length(n1);
                                ui =0;
                                for jj = 1:mk
                                    ui = ui + n1(jj)*sum(n2(jj+1:mk))+n1(jj)*n2(jj)/2;
                                end
                                
                                mauroc(j,k,kk,k2,1,1) = mauroc(j,k,kk,k2,1,1) + 0.5-ui;
                                
                                for nn = 1:nsu
                                    usps1 = zeros(num,1);
                                    usps1(1:num/2) = nsps1l(tyy(nn,1:num/2,1));
                                    usps1(num/2+1:num) = nsps1u(tyy(nn,1:num/2,2));
                                    usps2 = zeros(num,1);
                                    usps2(1:num/2) = nsps2l(tyy(nn,1:num/2,1));
                                    usps2(num/2+1:num) = nsps2u(tyy(nn,1:num/2,2));
                                    
                                    unsps = [usps1 usps2];
                                    dt = usps2-usps1;
                                    tdt = find(dt>0);
                                    tdt2 = find(dt==0);
                                    
                                    muauroc(j,k,kk,k2,nn,1,2) = muauroc(j,k,kk,k2,nn,1,2) +  0.5-length(tdt)/(num-length(tdt2));
                                
                                    ttd = find(abs(squeeze(unsps(:,1))-squeeze(unsps(:,2)))>0);
                                    lttd = length(ttd);

                                    n1 = histc(squeeze(unsps(ttd,1)),ib)/lttd;
                                    n2 = histc(squeeze(unsps(ttd,2)),ib)/lttd;
                                    mk = length(n1);
                                    ui =0;
                                    for jj = 1:mk
                                        ui = ui + n1(jj)*sum(n2(jj+1:mk))+n1(jj)*n2(jj)/2;
                                    end
                                    
                                    muauroc(j,k,kk,k2,nn,1,1) = muauroc(j,k,kk,k2,nn,1,1) +0.5-ui;
                                end
                                
                                
                                %%upper quartile
                                
                                dt = nsps2u-nsps1u;
                                tdt = find(dt>0);
                                tdt2 = find(dt==0);
                                 
                                mauroc(j,k,kk,k2,2,2) = mauroc(j,k,kk,k2,2,2) +0.5-length(tdt)/(num-length(tdt2));
                                
                                ttd = find(abs(squeeze(nspsu(:,1))-squeeze(nspsu(:,2)))>0);
                                lttd = length(ttd);

                                n1 = histc(squeeze(nspsu(ttd,1)),ib)/lttd;
                                n2 = histc(squeeze(nspsu(ttd,2)),ib)/lttd;
                                mk = length(n1);
                                ui =0;
                                for jj = 1:mk
                                    ui = ui + n1(jj)*sum(n2(jj+1:mk))+n1(jj)*n2(jj)/2;
                                end
                                
                                mauroc(j,k,kk,k2,2,1) = mauroc(j,k,kk,k2,2,1) + 0.5-ui;
                                
                                for nn = 1:nsu
                                    usps1 = zeros(num,1);
                                    usps1(1:num/2) = nsps1l(tyy(nn,num/2+1:num,1));
                                    usps1(num/2+1:num) = nsps1u(tyy(nn,num/2+1:num,2));
                                    usps2 = zeros(num,1);
                                    usps2(1:num/2) = nsps2l(tyy(nn,num/2+1:num,1));
                                    usps2(num/2+1:num) = nsps2u(tyy(nn,num/2+1:num,2));
                                    
                                    unsps = [usps1 usps2];
                                    dt = usps2-usps1;
                                    tdt = find(dt>0);
                                    tdt2 = find(dt==0);
                                     
                                    muauroc(j,k,kk,k2,nn,2,2) = muauroc(j,k,kk,k2,nn,2,2) +0.5-length(tdt)/(num-length(tdt2));
                                
                                ttd = find(abs(squeeze(unsps(:,1))-squeeze(unsps(:,2)))>0);
                                lttd = length(ttd);

                                n1 = histc(squeeze(unsps(ttd,1)),ib)/lttd;
                                n2 = histc(squeeze(unsps(ttd,2)),ib)/lttd;

                                    mk = length(n1);
                                    ui =0;
                                    for jj = 1:mk
                                        ui = ui + n1(jj)*sum(n2(jj+1:mk))+n1(jj)*n2(jj)/2;
                                    end
                                     
                                    muauroc(j,k,kk,k2,nn,2,1) = muauroc(j,k,kk,k2,nn,2,1) +0.5-ui;
                                end
                            end
                        end
                    end
                end
                
             end
            
                if cx>0
                    numsessionssel(k)= cx;
                    mauroc(:,k,:,:,:,:) = mauroc(:,k,:,:,:,:)/cx;
                    muauroc(:,k,:,:,:,:,:) = muauroc(:,k,:,:,:,:,:)/cx;
                end
        end
        tl = sprintf('save /home/chamanthi/Thiele_data/datacorrel25/datacorrel%s_%s2  numsessionssel mauroc muauroc selsess',animal,area);
        eval(tl);
    end
end


