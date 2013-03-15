function RC=STRF_script(CTXpath,NEVpath,NSEpath,NCSpath)
% function RC=STRF_script(CTXpath,NEVpath,NCSpath)

% NSEpath = 'D:\';
CTXname = 'STRF072S';% dont change, identifier for paradigm parameters
NSEname = 'Sc1';% tag for spike data
NCSname = 'MUA';% tag for csc data

NEV = NLX_LoadNEV(NEVpath,'FULL',[],[]);
NEV = NLX_getEventType(NEV,'USER');

CTXpar = ctx_ParameterSwitch(CTXname);
NLXpar = nlx_ParameterSwitch(CTXname);
nameInd=find(CTXpath=='\');
ctxWholeName=CTXpath(nameInd(end)+1:end);
onTime=[];offTime=[];
for i=1:length(NEV.Eventstring)
    if strcmp([ctxWholeName,' on'],NEV.Eventstring(i))
        onTime=[onTime NEV.TimeStamps(i)];
    end
    if strcmp([ctxWholeName,' off'],NEV.Eventstring(i))
        offTime=[offTime NEV.TimeStamps(i)];
    end
end
if length(onTime)>1||length(offTime)>1
    error('More than one encode for onTime and/or offTime found in event file!!')
end
NLXTime = [onTime offTime];

[SPK,P] = MakeSPK_2(NEVpath,NLXTime,NLXpar,CTXpath,CTXpar, ...
    'TrialEventWin',{'NLX_SUBJECT_START' 'NLX_TRIAL_END'},...
    'TrialEventOffset',[-500000 500000],...% in microsec
    'IncludeTrialEvents',{'NLX_READ_DATA'},...% these events are mandatory for each trial
    'ExcludeTrialEvents',{},...% occurence of these events exclude trial
    'NSEpath',{NSEpath},...
    'NSEname',{NSEname},...
    'NCSpath',{NCSpath},...
    'NCSname',{NCSname},...
    'ErrorTrials',[],...
    'ConvertEOGX',15.0*2/4096,...
    'ConvertEOGY',11.3*2/4096);



% % % % % ChanName = 'Sc1.01';
ChanName = NCSname;

% BasePar = {'NLX_SUBJECT_START' 'NLX_STIM_ON' []};% defines baseline window between two events
BasePar = {'NLX_SUBJECT_START' [0] [300]};% defines baseline window as time window around reference event
% BasePar = 0;% defines baseline as MapTau = 0 ms;
MapTau = []; %%%[100];
MapTauWin = []; %%%[-50 0];


[RC,Stim,M] = STRF_rate(SPK,CTXname,ChanName,true, ...
    'BasePar',BasePar, ...
    'MapTau',MapTau, ...
    'MapTauWin',MapTauWin);
