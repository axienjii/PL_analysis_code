function PL_LFP_list_xing(animal,area)
% if nargin<3 || isempty(channels)
%     channels = main_channels(animal,area);
% end
% if nargin<4 || isempty(sessions)
%     sessions = main_raw_sessions_final(animal,area,1,0);
% end
% for j=1:length(sessions)
%     for i=1:length(channels)
%         perceptual_learning_LFP_xing(animal,area,channels(i),sessions(j), 1,[1],1) %%%..
%     end
% end

%split session 355:
%Run on 355.1 and 355.2.
animal='blanco';
area='v1_1';
channels = main_channels(animal,area);
[dummy,testContrasts]=session_metadata(355.1,animal);
for i=1:length(channels)
    loadText=['load M:\Xing\pl_LFP\blanco\355.1\SGL_trial_LFP_-512_1536_ch',num2str(channels(i)),'.mat'];
    eval(loadText)
    LFPdata1=LFP_data;
    event_arr1=event_arr;
    loadText=['load M:\Xing\pl_LFP\blanco\355.2\SGL_trial_LFP_-512_1536_ch',num2str(channels(i)),'.mat'];
    eval(loadText)
    LFPdata2=LFP_data;
    event_arr2=event_arr;
    event_arr=event_arr1+event_arr2;
    LFP_data=zeros(length(testContrasts(1,:)),max(event_arr),2049); %%% downsampled to 1 kHz
    for h=1:length(testContrasts(1,:))
        LFP_data(h,1:event_arr(1,h),:)=[squeeze(LFPdata1(h,1:event_arr1(1,h),:));squeeze(LFPdata2(h,1:event_arr2(1,h),:))];
    end
    saveText=['save M:\Xing\pl_LFP\blanco\355\SGL_trial_LFP_-512_1536_ch',num2str(channels(i)),'.mat LFP_data event_arr'];
    eval(saveText)
end


