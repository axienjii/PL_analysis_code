function plot_all_session_rasters2
%Written by Xing 20/08/12.
%Routine to plot rasters across all sessions, for a given channel.


animal = 'jack';
animal = 'blanco';
area = 'v1';
area = 'v4';
channels=[1 2 3 18 33 34 36 37 38 40 42 49 50 51 52 53 54 55 57 58 59 60];
channels=1;
readMat=0;
readNSE=1;
lineCounter=0;
duration=[529 529 529 529 400];
setThreshold=tic
sessions=main_raw_sessions(animal, area);
sessions=307;
for chNum=1:length(channels)
    sessionLabels=[];
    meanSpontanRate=[];
    mean2Rate=[];
    mean3Rate=[];
    mean4Rate=[];
    mean5Rate=[];
    channel=channels(chNum);
    figtest=figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.2, 0.04, 0.5, 0.9]);
    set(figtest,'PaperUnits','centimeters','PaperType','A4','PaperOrientation', 'portrait', 'PaperPosition', [0.63452 0.63452 21 28.41]);
    for sessionCount=1:length(sessions)
        overlap=[];
        numSpikes1=[];
        numSpikes2=[];
        numSpikes3=[];
        numSpikes4=[];
        numSpikes5=[];
        session=sessions(sessionCount);
        if readMat==1
            if session==333
                %to remove empty cells from end of matarray for session 333:
                loadText=['load F:\blanco\mat_files\',num2str(channel),'_333.mat'];
                eval(loadText)
                for cond=1:size(matarray,1)
                    for i=1:length(matarray{cond,1})
                        if isempty(matarray{cond,1}{i})
                            for epoch=1:5
                                matarray{cond,epoch}=matarray{cond,epoch}(1:i-1)
                            end
                            break
                        end
                    end
                end
                saveText=['save F:\blanco\mat_files\',num2str(channel),'_333.mat matarray'];
                eval(saveText)
            end
            currdir=cd
            cd 'F:\blanco\mat_files'
            if round(channel)~=channel
                filename2=[num2str(ch),'_',num2str(10*(channel-round(channel)))];
            else
                filename2=[num2str(channel)];
            end
            loadtext=sprintf('load %s_%s.mat matarray',filename2,num2str(session))
            eval(loadtext);
            chdirtext=sprintf('cd ''%s''',currdir);
            eval(chdirtext);
            cd
            %cd to folder containing corresponding .ncs file for that channel and
            %session.
            for cond=1:14
                %subplot(1,5,epoch);
                for epoch=1:5
                    for n=1:length(matarray{cond,epoch})
                        for m=1:length(matarray{cond,epoch}{n})%raster plots
                            plot([matarray{cond,epoch}{n}(m) matarray{cond,epoch}{n}(m)],[n+lineCounter-0.5 n+lineCounter+0.5],'k');hold on
                            if matarray{cond,epoch}{n}(m)>0&&epoch==1
                                overlap=[overlap matarray{cond,epoch}{n}(m)];
                            end
                            if epoch==1
                                numSpikes1=[numSpikes1 length(matarray{cond,epoch}{n})];%list the number of spikes per trial during spontaneous period
                            end
                        end
                    end
                end
                lineCounter=lineCounter+length(matarray{cond,epoch});%add up number of trials across sessions
            end
            meanSpontanRate=[meanSpontanRate mean(numSpikes1)/529*1000];%list of spontan rate during each session, across trials and conditions
            sessionLabels=[sessionLabels;session lineCounter];
        elseif readNSE==1
%             step=1;
%             folder=['F:\blanco\before_300113\v4_1_sorted_spikes\4\',num2str(session)];%monitor artifact present
            step=2;
            step=3;
            folder=['I:\pl_spnorm_nse\blanco\',num2str(session)];%monitor artifact removed, spontaneous activity matched across sessions            
            listing=dir(folder);
            if size(listing,1)>2
                for q=3:size(listing,1)
                    fileName=listing(q,1).name;
                    ind=find(fileName=='_');
                    if length(ind)==2
                        if step==1
                            ch=fileName(ind(1)+1:ind(2)-1);
                        elseif step==2
                            if strcmp(fileType,'spk__spnorm.nse')
                                ch=fileName(ind(1)+1:ind(2)-1);
                            else
                                ch=[];
                            end
                            elseif step==3
                                fileType=[fileName(1:ind(1)),fileName(ind(2):end)];
                                if strcmp(fileType,'SpikeCh__.NSE')
                                    ch=fileName(ind(1)+1:ind(2)-1);
                                end
                        end
                        if strcmp(ch,num2str(channel))
                            nse_fname=[folder,'\',fileName];
                            file_of_int = session_metadata(session, animal);
                            nev_fname = nev_finder(animal, area, session);
                            vals = psycho_EV_vals(nev_fname, file_of_int);
                            alltp_sptimes = make_alltp_sptimesmf_keepvoltage(nse_fname, vals);
                            numTrialsPerSess=0;
                            for n=1:size(alltp_sptimes{1},2)
                                temp=alltp_sptimes{1}(:,n);
                                rasters1=temp(~isnan(temp));
                                temp=alltp_sptimes{2}(:,n);
                                rasters2=temp(~isnan(temp));
                                temp=alltp_sptimes{3}(:,n);
                                rasters3=temp(~isnan(temp));
                                temp=alltp_sptimes{4}(:,n);
                                rasters4=temp(~isnan(temp));
                                temp=alltp_sptimes{5}(:,n);
                                rasters5=temp(~isnan(temp));
                                if ~isempty(rasters1)||~isempty(rasters2)||~isempty(rasters3)||~isempty(rasters4)||~isempty(rasters5)%make sure it's not an empty trial
                                    numTrialsPerSess=numTrialsPerSess+1;
                                    rasters1=-rasters1;
                                    numSpikes1=[numSpikes1 length(rasters1)];%list the number of spikes per trial during spontaneous period
                                    numSpikes2=[numSpikes2 length(rasters2)];%list the number of spikes per trial during spontaneous period
                                    rasters3=rasters3+529;
                                    rasters3=rasters3(rasters3<1058);%due to variable length of sample-test interval
                                    numSpikes3=[numSpikes3 length(rasters3)];%list the number of spikes per trial during spontaneous period
                                    rasters4=rasters4+529*2;
                                    numSpikes4=[numSpikes4 length(rasters4)];%list the number of spikes per trial during spontaneous period
                                    rasters5=rasters5+529*3;
                                    numSpikes5=[numSpikes5 length(rasters5)];%list the number of spikes per trial during spontaneous period
                                    for m=1:length(rasters1)%raster plots
                                        plot([rasters1(m) rasters1(m)],[numTrialsPerSess+lineCounter-0.5 numTrialsPerSess+lineCounter+0.5],'k');hold on
                                    end
                                    for m=1:length(rasters2)%raster plots
                                        plot([rasters2(m) rasters2(m)],[numTrialsPerSess+lineCounter-0.5 numTrialsPerSess+lineCounter+0.5],'k');hold on
                                    end
                                    for m=1:length(rasters3)%raster plots
                                        plot([rasters3(m) rasters3(m)],[numTrialsPerSess+lineCounter-0.5 numTrialsPerSess+lineCounter+0.5],'k');hold on
                                    end
                                    for m=1:length(rasters4)%raster plots
                                        plot([rasters4(m) rasters4(m)],[numTrialsPerSess+lineCounter-0.5 numTrialsPerSess+lineCounter+0.5],'k');hold on
                                    end
                                    for m=1:length(rasters5)%raster plots
                                        plot([rasters5(m) rasters5(m)],[numTrialsPerSess+lineCounter-0.5 numTrialsPerSess+lineCounter+0.5],'k');hold on
                                    end
                                end
                            end
                            lineCounter=lineCounter+numTrialsPerSess;%add up number of trials across sessions
                            meanSpontanRate=[meanSpontanRate mean(numSpikes1)/529*1000];%list of spontan rate during each session, across trials and conditions
                            mean2Rate=[mean2Rate mean(numSpikes2)/529*1000];%list of spontan rate during each session, across trials and conditions
                            mean3Rate=[mean3Rate mean(numSpikes3)/529*1000];%list of spontan rate during each session, across trials and conditions
                            mean4Rate=[mean4Rate mean(numSpikes4)/529*1000];%list of spontan rate during each session, across trials and conditions
                            mean5Rate=[mean5Rate mean(numSpikes5)/400*1000];%list of spontan rate during each session, across trials and conditions
                            sessionLabels=[sessionLabels;session lineCounter];
                        end
                    end
                end
            end
            for timePeriod=1:5
                tally=sum(~isnan(alltp_sptimes{timePeriod}),1);
                zeroInd=find(tally~=0);
                if ~isempty(zeroInd)
                    tally=tally(1:zeroInd(end));
                end
                meanTally(timePeriod)=mean(tally)/duration(timePeriod)*1000;
                SDTally(timePeriod)=std(tally);
            end
        end
    end
    mean2Rate
    mean3Rate
    mean4Rate
    mean5Rate
    ylabel(sessionLabels);
    set(gca,'XTick',[-529 0 529 1058 1587 1987]);
    set(gca,'XTickLabel',[-529 0 529 1058 1587 1987]);
    set(gca,'YTick',sessionLabels(:,2));
    set(gca,'YTickLabel',sessionLabels(:,1));
    for i=1:length(meanSpontanRate)
        ptext=sprintf('%.2f     ',meanSpontanRate(i));
        text('Position',[-900 sessionLabels(i,2)],'FontSize',9,'String',ptext);
    end
    xlim([-529 1987]);
    filename1=['F:\blanco\auto_spontan_analysis\ch',num2str(channel),'_act_artifacts_removed_sorted3'];
    printtext=sprintf('print -dpng %s_act',filename1);
    eval(printtext);
    toc(setThreshold)
end

% overlap=[];
% for epoch=1:5
%     for cond=1:14
%         for n=1:length(matarray{cond,epoch})
%             for m=1:length(matarray{cond,epoch}{n})%raster plots
%                 if matarray{cond,epoch}{n}(m)>0&&epoch==1
%                     overlap=[overlap matarray{cond,epoch}{n}(m)];
%                 end
%             end
%         end
%     end
% end