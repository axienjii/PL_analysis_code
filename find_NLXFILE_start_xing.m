function [on_time, off_time]=find_NLXFILE_start_xing(file_of_int,eventsFilePath)

if nargin<3
[TimeStamps, EventIDs, Nttls, Extras, EventStrings, Header] = Nlx2MatEV(eventsFilePath, [1 1 1 1 1],1,1,1);
else
    TimeStamps=varargin{1};
    event_arr_NLX=varargin{2};
    EventStrings=varargin{3};
end
on_time=[];
off_time=[];
onfile=sprintf('%s on',file_of_int);
offfile=sprintf('%s off',file_of_int);
onfile2=sprintf('%s on',file_of_int);
offfile2=sprintf('%s off',file_of_int);
loc_on=find(strcmp(EventStrings,onfile));
if isempty(loc_on)
    loc_on=TimeStamps(1);
else
    on_time=TimeStamps(loc_on(1));
end

loc_off=find(strcmp(EventStrings,offfile));
if ~isempty(loc_off)
    off_time=TimeStamps(loc_off(1));
else
    off_time=TimeStamps(end);
end
% for j=1:length(EventStrings)
%     prl=EventStrings(j);
%     x=strfind(prl,'on');
%     if ~isempty(x{1});
%         prl=EventStrings(j);
%         k=strcmp(prl,onfile);
%         k2=strcmp(prl,onfile2);
%         if k==1 || k2==1
%             on_time=TimeStamps(j);
%         end
%         
%     end
%     k=strcmp(prl,offfile); 
%     k2=strcmp(prl,offfile2); 
%     if k==1 || k2==1
%         off_time=TimeStamps(j);
%     end
%     x=strfind(prl,'off');
%     if ~isempty(x{1});
%         prl=EventStrings(j);
%     end
% end
if isempty(on_time)
    on_time=TimeStamps(1);
end
if isempty(off_time)
    off_time=TimeStamps(end);
end
