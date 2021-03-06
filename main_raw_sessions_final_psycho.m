function sessions = main_raw_sessions_final_psycho(animal, area, expt_type, horz)
% A good set of raw sessions to analyse
% For all of these,
% sample = 30;
% 14 conditions
% Data is collected on consecutive days

% Default values
if nargin<3
    expt_type = 1;
end

animal = lower(animal);
area = lower(area);

if ispc
    switch [animal,area]
        case ['blanco','v4_0_1']
            sessions=304;
        case ['blanco','v4_0_2']
            sessions=305;
        case ['blanco','v4_0_3']
            sessions=306;
        case ['blanco','v4']
            %                 sessions = 306:342;
            % We know there is missing data for this, so just use the
            % sessions which actually exist.
            %                 sessions = intersect(sessions,unique(floor(sessions_available(animal,area))));
            % sessions is really
            sessions = [307 308 311 312 313 314 316 317 318 320 321 322 323 324 327 328 329:341];
            if horz==1
            sessions = [307 308 311 312 313 314 316 317 318 320 321 322 323 324 327 328 329:342];
            end
            %                 sessions = 329:341;
            %                 conditions = [10 15 20 25 27 28 29 31 32 33 35 40 50 60];
        case ['blanco','v4_1']
            sessions = [307 308 311 312 313 314 316 317 318 320 321 322 323 324 327 328 329:341];
            if horz==1
            sessions = [307 308 311 312 313 314 316 317 318 320 321 322 323 324 327 328 329:342];
            end
            
        case ['blanco','v1']
            % NB: session 355 is in two parts; 355_1 and 355_2
            % and the session metadata for these is as 355.1 and 355.2
            sessions = 343:359;
            %                 conditions = [5 10 15 20 22 25 28 32 35 40 45 50 60 90];
            
          case ['blanco','v1_1']
            sessions = 343:359;
            
        case ['blanco','v4_2']
            sessions = [360 361 362 363 364];%have to restore and play back 362 from tape
            
        case ['blanco','v1_2']
            sessions = [388:401 403:422 431:434 435 436:443 451:459]; %no flankers: 388 to 422. flankers: 431:452 no flankers: 453 to 459
               
        case ['blanco','v1_2_1']
            sessions = [388:401 403:422]; %no flankers: 388 to 422. flankers: 431:452 no flankers: 453 to 459
               
        case ['blanco','v1_2_2']
            sessions = [431:434 435 436:443 451:452]; %no flankers: 388 to 422. flankers: 431:452 no flankers: 453 to 459
               
        case ['blanco','v1_2_3']
            sessions = 453:459; %no flankers: 388 to 422. flankers: 431:452 no flankers: 453 to 459
            
        case ['jack','v4_0_1']
            sessions=22;
            
        case ['jack','v4_0_2']
            sessions=23;
            
        case ['jack','v4']
            % All present and correct
            sessions = [24 25 27:49];
            if horz==1
                sessions = [24 25 27:50];
            end
            %                 conditions = [10 15 20 25 27 28 29 31 32 33 35 40 50 60];
            
        case ['jack','v4_1']
            % All present and correct
            sessions = [24 25 27:49];
            if horz==1
                sessions = [24 25 27:50];
            end
            
        case ['jack','v1']
            % All present and correct
            sessions = 51:72;
            %                 conditions = [5 10 15 20 22 25 28 32 35 40 45 50 60 90];
        
        case ['jack','v1_1']
            % All present and correct
            sessions = 51:72;
            %                 conditions = [5 10 15 20 22 25 28 32 35 40 45 50 60 90];
            
        case ['jack','v4_2']
            % All present and correct
            sessions = 73:77;
            %                 conditions = [10 15 20 25 27 28 29 31 32 33 35 40 50 60];
            
        case ['jack','v1_2']
            sessions = 78:119; % 
            
        case ['jack','v1_2_1']
            sessions = 78:93; % 
            
        case ['jack','v1_2_2']
            sessions = 94:115; % 
            
        case ['jack','v1_2_3']
            sessions = 116:119; % 
            
        case ['jack','v1_3']
            sessions = []; % BU/PL task
            
        case ['jack','v1_4']
            sessions = 140:182; % control task with stimuli at Blanco's V1 location
            
        case ['jack','v1_4_1']
            sessions = 140:162; % control task with stimuli at Blanco's V1 location
            
        case ['jack','v1_4_2']
            sessions = 163:184; % control task with stimuli at Blanco's V1 location
            
        case ['jack','v1_4_3']
            sessions = 185:189; % control task with stimuli at Blanco's V1 location
            
    end
else
    switch expt_type
        case 0
            % Pre-trianing
            
            switch [animal,area]
                case ['blanco','v4']
                    sessions = []; % Not empty!
                case ['blanco','v1']
                    sessions = [];
                    
                case ['jack','v4']
                    sessions = 1:23;
                case ['jack','v1']
                    sessions = [];
            end
            
        case 1
            % Main experiment
            
            switch [animal,area]
                case ['blanco','v4']
                    %                 sessions = 306:342;
                    % We know there is missing data for this, so just use the
                    % sessions which actually exist.
                    %                 sessions = intersect(sessions,unique(floor(sessions_available(animal,area))));
                    % sessions is really
                    sessions = [307 308 , 311 , 313 314 , 317 318 , 320 321 , 329:341];
                    %                 sessions = 329:341;
                    %                 conditions = [10 15 20 25 27 28 29 31 32 33 35 40 50 60];
                    
                case ['blanco','v1']
                    % NB: session 355 is in two parts; 355_1 and 355_2
                    % and the session metadata for these is as 355.1 and 355.2
                    sessions = [343:354 355.1 355.2 356:359];
                    %                 conditions = [5 10 15 20 22 25 28 32 35 40 45 50 60 90];
                    
                case ['blanco','v4_2']
                    sessions = [360 361 363 364];%have to restore and play back 362 from tape
                    
                case ['jack','v4']
                    % All present and correct
                    sessions = [24 25 27:49];
                    % Exclude session 26
                    %                 conditions = [10 15 20 25 27 28 29 31 32 33 35 40 50 60];
                    
                case ['jack','v1']
                    % All present and correct
                    sessions = 51:72;
                    %                 conditions = [5 10 15 20 22 25 28 32 35 40 45 50 60 90];
            end
            
        case 2
            % Roving experiment
            switch [animal,area]
                case ['blanco','v1_2']
                    sessions = [388:410 412:422]; %do 411 when data available
                case ['jack','v1_2']
                    sessions = []; % Not empty!
            end
            
        otherwise
            error('Main sessions not defined for expt_type %s',num2str(expt_type));
    end
end
end