function fname = trial_corr_artifact_finder(animal, session)

animal = lower(animal);

if ispc    
    trialcorr_artifact_filename=[num2str(session),'_corrtrialartifact.mat'];
    fname = fullfile('F:','PL','pl_corr_art_trials',animal,trialcorr_artifact_filename);
end
