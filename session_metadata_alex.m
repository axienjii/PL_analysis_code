function [tag_on, conditions, sample_contrasts, expt_type, rotated, area, rovingType,tag_off] = session_metadata(session_num, animal)
% Look up function for the metadata involved in different sessions for the
% animals, such as which test contrasts were used on that day

% sample_contrast = 30;

% version2 taken from CSC_Ext.m

% Only 2 sessions have rotated stimuli. Default, not rotated.
rotated = 0;
rovingType=0;
sample_contrasts=30;


% -------------------------------------------------------------------------
% For blanco
% -------------------------------------------------------------------------
if strcmpi(animal,'blanco')

if session_num==304
% REC_PEN    = [304];
% V_Channels = [1 2 3 4 7 12 13 14 15 18 20 22 24 33 34 36 37 38 40 42 49 50 51 52 53 54 55 57 59 60]; 
conditions = [5 10 20 25 35 40 60 90];
tag_on     = '21653090.1';
tag_off    = '21653090.1';
area       = 'v4';
expt_type  = 1;

elseif session_num==305
% REC_PEN    = [305];
% V_Channels = [1 2 3 4 7 12 13 14 15 18 20 22 24 33 34 36 37 38 40 42 49 50 51 52 53 54 55 57 59 60]; 
conditions = [10 15 20 25 27 29 31 33 35 40 50 60];
tag_on     = '216136.2';
tag_off    = '216136.2';
area       = 'v4';
expt_type  = 1;

elseif session_num>=306 && session_num<=312
% REC_PEN    = [306 307 308 311 312];
% REC_PEN    = [306];
% V_Channels = [1 2 3 4 7 12 13 14 15 18 20 22 24 33 34 36 37 38 40 42 49 50 51 52 53 54 55 57 59 60]; 
conditions = [10 15 20 25 27 28 29 31 32 33 35 40 50 60];
tag_on     = '21613614.1';
tag_off    = '21613614.1';
expt_type  = 1;
area       = 'v4';

elseif session_num==313
% REC_PEN    = [313];
% V_Channels = [1 2 3 4 7 12 13 14 15 18 20 22 24 33 34 36 37 38 40 42 49 50 51 52 53 54 55 57 59 60]; 
conditions = [10 15 20 25 27 28 29 31 32 33 35 40 50 60];
tag_on     = '21613614.2';
tag_off    = '21613614.2';
expt_type  = 1;
area       = 'v4';

elseif session_num>=314 && session_num<=341
% REC_PEN    = [314 316 317 318 320 321 322 323 324 327 328 329 ...
%              330 331 332 333 334 335 336 337 338 339 340 341 342];
% V_Channels = [1 2 3 4 7 12 13 14 15 18 20 22 24 33 34 36 37 38 40 42 49 50 51 52 53 54 55 57 59 60]; 
conditions = [10 15 20 25 27 28 29 31 32 33 35 40 50 60];
tag_on     = '21613614.1';
tag_off    = '21613614.1';
expt_type  = 1;
area       = 'v4';

elseif session_num==342
% **** NOTE: session 342 has rotated stimuli! ****
% V_Channels = [1 2 3 4 7 12 13 14 15 18 20 22 24 33 34 36 37 38 40 42 49 50 51 52 53 54 55 57 59 60]; 
conditions = [10 15 20 25 27 28 29 31 32 33 35 40 50 60];
tag_on     = '21613614.1';
tag_off    = '21613614.1';
expt_type  = 1;
rotated    = 1;
area       = 'v4';

elseif session_num>=343 && session_num<=355.1
% REC_PEN    = [343 344 345 346 347 348 349 350 351 352 353 354 355.1];
% V_Channels = [8 9 10 11 15 17 19 21 23 25 26 27 28 29 31 44 45 46 48 61 62 63 64];
conditions = [5 10 15 20 22 25 28 32 35 40 45 50 60 90];
tag_on     = '2353914.1';
tag_off    = '2353914.1';
expt_type  = 1;
area       = 'v1';

elseif session_num==355.2
% REC_PEN    = [355.2];
% V_Channels = [8 9 10 11 15 17 19 21 23 25 26 27 28 29 31 44 45 46 48 61 62 63 64];
conditions = [5 10 15 20 22 25 28 32 35 40 45 50 60 90];
tag_on     = '2353914.2';
tag_off    = '2353914.2';
expt_type  = 1;
area       = 'v1';

elseif session_num>=356 && session_num<=359
% REC_PEN    = [356 357 358 359];
% V_Channels = [8 9 10 11 15 17 19 21 23 25 26 27 28 29 31 44 45 46 48 61 62 63 64];
conditions = [5 10 15 20 22 25 28 32 35 40 45 50 60 90];
tag_on     = '2353914.1';
tag_off    = '2353914.1';
expt_type  = 1;
area       = 'v1';

elseif session_num==360
% REC_PEN    = [360];
% V_Channels = [1 2 3 4 7 12 13 14 15 18 20 22 24 33 34 36 37 38 40 42 49 50 51 52 53 54 55 57 59 60]; 
conditions = [10 15 20 25 27 28 29 31 32 33 35 40 50 60];
tag_on     = '21613614.2';
tag_off    = '21613614.2';
expt_type  = 1;
area       = 'v4';

elseif session_num>=361 && session_num<=364
% REC_PEN    = [361 362 363 364];
% V_Channels = [1 2 3 4 7 12 13 14 15 18 20 22 24 33 34 36 37 38 40 42 49 50 51 52 53 54 55 57 59 60]; 
conditions = [10 15 20 25 27 28 29 31 32 33 35 40 50 60];
tag_on     = '21613614.1';
tag_off    = '21613614.1';
expt_type  = 1;
area       = 'v4';

elseif session_num>=374 && session_num<=375%bottom-up attention task, V4
conditions = 1:8;
tag_on     = 'BU_SAC.1';
tag_off    = 'BU_SAC.2';
expt_type  = 1;
area       = 'v4';

elseif session_num>=377 && session_num<=379%bottom-up attention task, V1
conditions = [1:8];
tag_on     = 'BU_SAC.1';
tag_off    = 'BU_SAC.2';
expt_type  = 1;
area       = 'v1';

elseif session_num>=380 && session_num<=384
% REC_PEN    = [380 381 382 383 384];
% V_Channels = [8 9 10 11 15 17 19 21 23 25 26 27 28 29 31 44 45 46 48 61 62 63 64];
conditions = [5 10 15 20 22 25 28 32 35 40 45 50 60 90];
tag_on     = '2353914.1';
tag_off    = '2353914.1';
expt_type  = 2;

elseif session_num>=385 && session_num<=387
% REC_PEN    = [385 386 387];
% V_Channels = [8 9 10 11 15 17 19 21 23 25 26 27 28 29 31 44 45 46 48 61 62 63 64];
conditions = [5 10 15 22 25 28 32 35 38 45 60 90];
tag_on     = '2253912.1';
tag_off    = '2253912.1';
expt_type  = 2;

elseif session_num>=388 && session_num<=405.1
% REC_PEN    = [388 389 390 391 392 393 394 395 396 397 398 399 400 ...
%              401 403 404 405.1]; add 402 for LFP
% V_Channels = [8 9 10 11 15 17 19 21 23 25 26 27 28 29 31 44 45 46 48 61 62 63 64];
conditions = [  5 10 15 22 25 28 32 35 38 45 60 90 ...
                5 10 12 15 18 22 25 28 35 45 60 90 ...
                5 10 15 25 32 35 38 42 45 50 60 90 ];
tag_on     = '225234912.1';
tag_off    = '225234912.1';
expt_type  = 2;
rovingType=1;%1st time without flankers

elseif session_num==405.2||session_num==4052
% REC_PEN    = [405.2];
% V_Channels = [8 9 10 11 15 17 19 21 23 25 26 27 28 29 31 44 45 46 48 61 62 63 64];
conditions = [  5 10 15 22 25 28 32 35 38 45 60 90 ...
                5 10 12 15 18 22 25 28 35 45 60 90 ...
                5 10 15 25 32 35 38 42 45 50 60 90 ];
tag_on     = '225234912.2';
tag_off    = '225234912.2';
expt_type  = 2;
rovingType=1;%1st time without flankers

elseif session_num>=406 && session_num<=414
% REC_PEN    = [406 407 408 409 410 411 412 413 414];
% V_Channels = [8 9 10 11 15 17 19 21 23 25 26 27 28 29 31 44 45 46 48 61 62 63 64];
conditions = [	5 10 15 22 25 28 32 35 38 45 60 90 ...
                5 10 12 15 18 22 25 28 35 45 60 90 ...
                5 10 15 25 32 35 38 42 45 50 60 90 ];
tag_on     = '225234912.1';
tag_off    = '225234912.1';
expt_type  = 2;
rovingType=1;%1st time without flankers

elseif session_num==415
% REC_PEN    = [415];
% V_Channels = [8 9 10 11 15 17 19 21 23 25 26 27 28 29 31 44 45 46 48 61 62 63 64];
conditions = [  5 10 15 22 25 28 32 35 38 45 60 90 ...
                5 10 12 15 18 22 25 28 35 45 60 90 ...
                5 10 15 25 32 35 38 42 45 50 60 90 ];
tag_on     = '225234912.2';
tag_off    = '225234912.2';
expt_type  = 2;
rovingType=1;%1st time without flankers

elseif session_num>=416 && session_num<=422
% REC_PEN    = [416 417 418 419 420 421 422];
% V_Channels = [8 9 10 11 15 17 19 21 23 25 26 27 28 29 31 44 45 46 48 61 62 63 64];
conditions = [  5 10 15 22 25 28 32 35 38 45 60 90 ...
                5 10 12 15 18 22 25 28 35 45 60 90 ...
                5 10 15 25 32 35 38 42 45 50 60 90 ];
tag_on     = '225234912.1';
tag_off    = '225234912.1';
expt_type  = 2;
rovingType=1;%1st time without flankers

elseif session_num>=423 && session_num<=430
% REC_PEN    = [423 424 425 426 427 428 429 430];
% V_Channels = [8 9 10 11 15 17 19 21 23 25 26 27 28 29 31 44 45 46 48 61 62 63 64];
conditions = [  25 35 ...
                15 25 ...
                35 45 ];
tag_on     = '222346.1';
tag_off    = '222346.1';
expt_type  = 2;

elseif session_num>=431 && session_num<=435.1
% REC_PEN    = [431 432 433 434 435.1];
% V_Channels = [8 9 10 11 15 17 19 21 23 25 26 27 28 29 31 44 45 46 48 61 62 63 64];
conditions = [  5 10 15 22 25 28 32 35 38 45 60 90 ...
                5 10 12 15 18 22 25 28 35 45 60 90 ...
                5 10 15 25 32 35 38 42 45 50 60 90 ];
tag_on     = '2223412f.1';
tag_off    = '2223412f.1';
expt_type  = 2;
rovingType=2;%flankers added

elseif session_num==435.2||session_num==4352
% REC_PEN    = [435.2];
% V_Channels = [8 9 10 11 15 17 19 21 23 25 26 27 28 29 31 44 45 46 48 61 62 63 64];
conditions = [  5 10 15 22 25 28 32 35 38 45 60 90 ...
                5 10 12 15 18 22 25 28 35 45 60 90 ...
                5 10 15 25 32 35 38 42 45 50 60 90 ];
tag_on     = '2223412f.2';
tag_off    = '2223412f.2';
expt_type  = 2;

elseif session_num>=436 && session_num<=452
% REC_PEN    = [436 437 438 439 440 441 442 443 450 451 452];
% V_Channels = [8 9 10 11 15 17 19 21 23 25 26 27 28 29 31 44 45 46 48 61 62 63 64];
conditions = [  5 10 15 22 25 28 32 35 38 45 60 90 ...
                5 10 12 15 18 22 25 28 35 45 60 90 ...
                5 10 15 25 32 35 38 42 45 50 60 90 ];
tag_on     = '2223412f.1';
tag_off    = '2223412f.1';
expt_type  = 2;
rovingType=2;%flankers added

elseif session_num>=453 && session_num<=459
% REC_PEN    = [453 454 455 456 457 458 459];
% V_Channels = [8 9 10 11 15 17 19 21 23 25 26 27 28 29 31 44 45 46 48 61 62 63 64];
conditions = [  5 10 15 22 25 28 32 35 38 45 60 90 ...
                5 10 12 15 18 22 25 28 35 45 60 90 ...
                5 10 15 25 32 35 38 42 45 50 60 90 ];
tag_on     = '225234912.1';
tag_off    = '225234912.1';
expt_type  = 2;
rovingType=3;%flankers removed

else
   ME = MException('SessionMetadata:UnknownSession', ...
    'No known session metadata for %s, session %s.',animal,num2str(session_num));
   throw(ME);
end


% -------------------------------------------------------------------------
% For jack
% -------------------------------------------------------------------------
elseif strcmpi(animal,'jack')

if session_num<22
conditions = [10 90];
tag_on     = '2161392.1';
tag_off    = '2161392.1';
expt_type  = 1;
area       = 'v4';

elseif session_num==22
% REC_PEN    = [22];
% V_Channels = [1 2 3 4 5 6 8 10 24 35 37 39 40 41 49 50 52 53 54 56];
conditions = [5 10 20 25 35 40 60 90];
tag_on     = '2161398.1';
tag_off    = '2161398.1';
expt_type  = 1;
area       = 'v4';

elseif session_num==23
% REC_PEN    = [23];
% V_Channels = [1 2 3 4 5 6 8 10 24 35 37 39 40 41 49 50 52 53 54 56];
conditions = [10 15 20 25 27 29 31 33 35 40 50 60];
tag_on     = '21613612.1';
tag_off    = '21613612.1';
expt_type  = 1;
area       = 'v4';

elseif session_num>=24 && session_num<=49
% REC_PEN    = [24 25 27 28 29 30 31 32 33 34 35 36 37 38 39 40 ...
%              41 42 43 44 45 46 47 48 49 50];
% V_Channels = [1 2 3 4 5 6 8 10 24 35 37 39 40 41 49 50 52 53 54 56];
conditions = [10 15 20 25 27 28 29 31 32 33 35 40 50 60];
tag_on     = '21613614.1';
tag_off    = '21613614.1';
expt_type  = 1;
area       = 'v4';

elseif session_num==50
% **** NOTE: session 50 has rotated stimuli! ****
% V_Channels = [1 2 3 4 5 6 8 10 24 35 37 39 40 41 49 50 52 53 54 56];
conditions = [10 15 20 25 27 28 29 31 32 33 35 40 50 60];
tag_on     = '21613614.1';
tag_off    = '21613614.1';
expt_type  = 1;
rotated    = 1;
area       = 'v4';

elseif session_num>=51 && session_num<=72
% REC_PEN    = [51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72];
% V_Channels = [7 9 11 12 13 14 15 16 17 18 19 20 21 22 23 25 26 27 28 29 30 31 32 51 55];
conditions = [5 10 15 20 22 25 28 32 35 40 45 50 60 90];
tag_on     = '47553914.1';
tag_off    = '47553914.1';
expt_type  = 1;
area       = 'v1';

elseif session_num>=73 && session_num<=77
% **** NOTE: Sinewave function with circular rectification presented instead of a Gabor function! ****
% REC_PEN    = [73 74 75 76 77];
% V_Channels = [1 2 3 4 5 6 8 10 24 35 37 39 40 41 49 50 52 53 54 56];
conditions = [10 15 20 25 27 28 29 31 32 33 35 40 50 60];
tag_on     = '21613614.1';
tag_off    = '21613614.1';
expt_type  = 1;
area       = 'v4';

elseif session_num>=78 && session_num<=93
% REC_PEN    = [79 80];
% V_Channels = [7 9 11 12 13 14 15 16 17 18 19 20 21 22 23 25 26 27 28 29 30 31 32 51 55];
conditions = [  5 10 15 22 25 28 32 35 38 45 60 90 ...
                5 10 12 15 18 22 25 28 35 45 60 90 ...
                5 10 15 25 32 35 38 42 45 50 60 90 ];
tag_on     = '47523412.1';
tag_off    = '47523412.1';
expt_type  = 2;
area       = 'v1';
rovingType=1;%1st time without flankers

elseif session_num>=94 && session_num<=115
% REC_PEN    = [79 80];
% V_Channels = [7 9 11 12 13 14 15 16 17 18 19 20 21 22 23 25 26 27 28 29 30 31 32 51 55];
conditions = [  5 10 15 22 25 28 32 35 38 45 60 90 ...
                5 10 12 15 18 22 25 28 35 45 60 90 ...
                5 10 15 25 32 35 38 42 45 50 60 90 ];
tag_on     = '4723412f.1';
tag_off    = '4723412f.1';
expt_type  = 2;
area       = 'v1';
rovingType=2;%flankers added

elseif session_num>=116 && session_num<=120%116 to 119- V1_2 task with flankers removed
% REC_PEN    = [79 80];
% V_Channels = [7 9 11 12 13 14 15 16 17 18 19 20 21 22 23 25 26 27 28 29 30 31 32 51 55];
conditions = [  5 10 15 22 25 28 32 35 38 45 60 90 ...
                5 10 12 15 18 22 25 28 35 45 60 90 ...
                5 10 15 25 32 35 38 42 45 50 60 90 ];
tag_on     = '47523412.1';
tag_off    = '47523412.1';
expt_type  = 2;
area       = 'v1';
rovingType=3;%flankers removed

elseif session_num==121
conditions = 1:4;
tag_on     = 'bu_pasv4.1';
tag_off    = 'bu_pasv4.1';
expt_type  = 1;
area       = 'v4';

elseif session_num==122
conditions = 1:4;
tag_on     = 'bu_pasv1.1';
tag_off    = 'bu_pasv1.1';
expt_type  = 1;
area       = 'v1';

elseif session_num==123
conditions = 1:4;
tag_on     = 'bu_rfrig.1';
tag_off    = 'bu_rfrig.1';
expt_type  = 1;
area       = 'v4';

elseif session_num==124
conditions = 1:8;
tag_on     = 'bu_jack.1';
tag_off    = 'bu_jack.1';
expt_type  = 1;
area       = 'v4';

elseif session_num==125
conditions = 1:8;
tag_on     = 'bu_jack.1';
tag_off    = 'bu_jack.1';
expt_type  = 1;
area       = 'v1';

elseif session_num>=140 && session_num<=162%116 to 119- Jack control task at Blanco's V1 location pre-flankers
% REC_PEN    = [79 80];
% V_Channels = [7 9 11 12 13 14 15 16 17 18 19 20 21 22 23 25 26 27 28 29 30 31 32 51 55];
conditions = [  5 10 15 22 25 28 32 35 38 45 60 90 ...
                5 10 12 15 18 22 25 28 35 45 60 90 ...
                5 10 15 25 32 35 38 42 45 50 60 90 ];
tag_on     = '225234912.1';
tag_off    = '225234912.1';
expt_type  = 2;
area       = 'v1';
rovingType=1;%pre-flankers 

elseif session_num>=163 && session_num<=184%116 to 119- Jack control task at Blanco's V1 location with flankers added
% REC_PEN    = [79 80];
% V_Channels = [7 9 11 12 13 14 15 16 17 18 19 20 21 22 23 25 26 27 28 29 30 31 32 51 55];
conditions = [  5 10 15 22 25 28 32 35 38 45 60 90 ...
                5 10 12 15 18 22 25 28 35 45 60 90 ...
                5 10 15 25 32 35 38 42 45 50 60 90 ];
tag_on     = '225234912f.1';
tag_off    = '225234912f.1';
expt_type  = 2;
area       = 'v1';
rovingType=2;%flankers added

elseif session_num>=185 && session_num<=189%116 to 119- Jack control task at Blanco's V1 location with flankers removed
% REC_PEN    = [79 80];
% V_Channels = [7 9 11 12 13 14 15 16 17 18 19 20 21 22 23 25 26 27 28 29 30 31 32 51 55];
conditions = [  5 10 15 22 25 28 32 35 38 45 60 90 ...
                5 10 12 15 18 22 25 28 35 45 60 90 ...
                5 10 15 25 32 35 38 42 45 50 60 90 ];
tag_on     = '225234912.1';
tag_off    = '225234912.1';
expt_type  = 2;
area       = 'v1';
rovingType=3;%flankers removed

% Roving Task.
% If condition number is between  1 and 12, sample contrast is 30%
% If condition number is between 13 and 24, sample contrast is 20%
% If condition number is between 25 and 36, sample contrast is 40%


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
else
   ME = MException('SessionMetadata:UnknownSession', ...
    'No known session metadata for %s, session %s.',animal,num2str(session_num));
   throw(ME);
end

% Unfamiliar animal
else
   ME = MException('SessionMetadata:UnknownAnimal', ...
       'No known metadata for animal %s.',animal);
   throw(ME);
end


% Start of file tag is the same as tag_on and tag_off
if ~strcmp(tag_on,tag_off)
    warning('SessionMetadata:NonMatchingTags' , ...
        'On and off tags do not match for %s, session %s' ,...
        animal, num2str(session_num));
end
file_of_int = tag_on;

% Test contrasts indexed by condition number
test_contrasts = conditions;

if expt_type  == 2;
    area       = 'v1';
end

% Sample contrasts indexed by condition number
if expt_type==1
    sample_contrasts = repmat(30,1,length(conditions));
elseif expt_type==2
    ncond_ea = length(conditions)/3;
    sample_contrasts = [repmat(30,1,ncond_ea) repmat(20,1,ncond_ea) repmat(40,1,ncond_ea)];
else
    error('Uncertain about the sample contrasts.');
end


end

%% old version
%
% version1 taken from various files
% --- jack_SE_V4_1_roc_batch.m
% --- blanco_V1_params.txt
% --- blanco_V4_params.txt
% sample_contrast = 30;
% 
% % For jack
% if strcmpi(animal,'jack')
%     if session_num==22
%         file_of_int='2161398.1';
%         test_contrasts=[10 15 20 25 35 40 60 90];
%     elseif session_num<22
%         file_of_int='2161392.1';
%         test_contrasts=[10 90];
%     elseif session_num==23
%         file_of_int='21613912.1';%wrongly named in Cheetah encode- actual name should be 21613612.1
%         test_contrasts=[10 15 20 25 27 29 31 33 35 40 50 60];
%     elseif session_num>23&&session_num<51||session_num>72
%         file_of_int='21613614.1';
%         test_contrasts=[10 15 20 25 27 28 29 31 32 33 35 40 50 60];
%     elseif session_num>50&&session_num<73
%         file_of_int='47553914.1';
%         test_contrasts=[5 10 15 20 22 25 28 32 35 40 45 50 60 90];
%     else
%        ME = MException('SessionMetadata:UnknownSession', ...
%         'No known session metadata for %s, session %s.',animal,num2str(session_num));
%        throw(ME);
%     end
% % For blanco
% elseif strcmpi(animal,'blanco')
%     % v4_1
%     if session_num==313
%         file_of_int='21613614.2';
%         test_contrasts=[10 15 20 25 27 28 29 31 32 33 35 40 50 60];
%     elseif session_num>=306 && session_num<=342
%         file_of_int='21613614.1';
%         test_contrasts=[10 15 20 25 27 28 29 31 32 33 35 40 50 60];
%     % v1 non-roving
%     elseif session_num>=343 && session_num<=359
%         file_of_int='2353914.1';
%         test_contrasts=[5 10 15 20 22 25 28 32 35 40 45 50 60 90];
%     % v4_2
%     elseif session_num==360
%         file_of_int='21613614.2';
%         test_contrasts=[10 15 20 25 27 28 29 31 32 33 35 40 50 60];
%     elseif session_num>360 && session_num<=364
%         file_of_int='21613614.1';
%         test_contrasts=[10 15 20 25 27 28 29 31 32 33 35 40 50 60];
%     % v1 roving without flankers
%     elseif (session_num>=388 && session_num<=422) || ...
%             (session_num>=453 && session_num<=459)
%         file_of_int='225234912.1';
%         % Unsure of sample_contrast and test_contrasts
%         sample_contrast = NaN;
%     % v1 roving with flankers
%     elseif (session_num>=431 && session_num<=452)
%         file_of_int='2223412f.1';
%         % Unsure of sample_contrast and test_contrasts
%         sample_contrast = NaN;
%     % Unknown session for blanco
%     else
%        ME = MException('SessionMetadata:UnknownSession', ...
%         'No known session metadata for %s, session %s.',animal,num2str(session_num));
%        throw(ME);
%     end
% % Unfamiliar animal
% else
%    ME = MException('SessionMetadata:UnknownAnimal', ...
%        'No known metadata for animal %s.',animal);
%    throw(ME);
% end
