% **** Batch file ****
function batch_identify_trialcorr_artifact_scott_arrays(animal, area, sessions, icopy, ncopies, istart)

if nargin<3 || isempty(sessions)
    sessions = main_raw_sessions_final(animal,area,[],1);
end
if nargin<5
    icopy = 1;
end
if nargin<6
    ncopies = 1;
end
if nargin<7
    istart = 1;
end

verbose = 10;

for iter = icopy:ncopies:length(sessions)
    
    if iter<istart
        continue;
    end
    
    session = sessions(iter);
    if sum(session==[355 405 435])==0
        if verbose; fprintf('%d(%d): %d/%d Processing session %.1f\n',...
                icopy, ncopies, iter, length(sessions), session);
        end
        input_fname=['F:\PL\pl_corr_art_trials\',animal,'\',num2str(session),'_corrtrialartifact.mat'];
        loadText=['load ',input_fname,' removeTrialsTimestamps rlist'];
        eval(loadText)
        output_fname=['F:\PL\pl_corr_art_trials\scott_copy\',animal,'\',num2str(session),'_corrtrialartifact.mat'];
        saveText=['save ',output_fname,' removeTrialsTimestamps'];
        eval(saveText)
    end
end

load('F:\PL\pl_corr_art_trials\scott_copy\blanco\355.1_corrtrialartifact.mat')
save 'F:\PL\pl_corr_art_trials\scott_copy\blanco\355.1_corrtrialartifact.mat' removeTrialsTimestamps
load('F:\PL\pl_corr_art_trials\scott_copy\blanco\355.2_corrtrialartifact.mat')
save 'F:\PL\pl_corr_art_trials\scott_copy\blanco\355.2_corrtrialartifact.mat' removeTrialsTimestamps
load('F:\PL\pl_corr_art_trials\scott_copy\blanco\405.1_corrtrialartifact.mat')
save 'F:\PL\pl_corr_art_trials\scott_copy\blanco\405.1_corrtrialartifact.mat' removeTrialsTimestamps
load('F:\PL\pl_corr_art_trials\scott_copy\blanco\405.2_corrtrialartifact.mat')
save 'F:\PL\pl_corr_art_trials\scott_copy\blanco\405.2_corrtrialartifact.mat' removeTrialsTimestamps
load('F:\PL\pl_corr_art_trials\scott_copy\blanco\435.1_corrtrialartifact.mat')
save 'F:\PL\pl_corr_art_trials\scott_copy\blanco\435.1_corrtrialartifact.mat' removeTrialsTimestamps
load('F:\PL\pl_corr_art_trials\scott_copy\blanco\435.2_corrtrialartifact.mat')
save 'F:\PL\pl_corr_art_trials\scott_copy\blanco\435.2_corrtrialartifact.mat' removeTrialsTimestamps
