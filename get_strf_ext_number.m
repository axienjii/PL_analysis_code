function extNum=get_strf_ext_number(session)
%list of extension numbers for each of the STRF Cortex files, e.g. '.1' for
%'strf072s.1'
switch(session)
    case 307
        extNum=[2 4];%2:SF=[0.125 0.25 0.5] & 4:SF=[1 2 4]
    case 308
        extNum=[1 0];%1:SF=[0.125 0.25 0.5] & no data for higher SFs
    case 311
        extNum=[1 2];
    case 312
        extNum=[1 0];%or [1 4]
    case 313
        extNum=[1 0];%find the second strf072s file! (.2)
    case 314
        extNum=[1 2];
    case 316
        extNum=[1 1];%forgot to create new file for 2nd set
    case 317
        extNum=[1 2];
    case 318
        extNum=[1 3];
    case 320
        extNum=[1 2];
    case 321
        extNum=[1 2];
    case 322
        extNum=[1 2];%events file missing
    case 323
        extNum=[1 2];
    case 324
        extNum=[1 2];
    case 327
        extNum=[1 2];
    case 328
        extNum=[1 2];
    case 329
        extNum=[1 2];
    case 330
        extNum=[1 2];
    case 331
        extNum=[1 2];
    case 332
        extNum=[1 2];
    case 333
        extNum=[1 2];
    case 334
        extNum=[1 2];
    case 335
        extNum=[1 2];
    case 336
        extNum=[1 3];%the .1 file has 2 'on' encodes- the second encode is the right one
    case 337
        extNum=[1 2];
    case 338
        extNum=[1 2];
    case 339
        extNum=[1 2];
    case 340
        extNum=[1 2];
    case 341
        extNum=[1 2];
    case 342
        extNum=[1 2];
    case 360
        extNum=[1 2];
    case 361
        extNum=[2 3];
    case 362
        extNum=[1 2];
    case 363
        extNum=[1 2];
    case 364
        extNum=[1 2];
    case 343
        extNum=[2];
    case 344
        extNum=[2];
    case 345
        extNum=[1];
    case 346
        extNum=[1];
    case 347
        extNum=[1];
    case 348
        extNum=[1];
    case 349
        extNum=[1];
    case 350
        extNum=[1];
    case 351
        extNum=[1];
    case 352
        extNum=[1];
    case 353
        extNum=[1];
    case 354
        extNum=[1];
    case 355
        extNum=[1];
    case 356
        extNum=[1];
    case 357
        extNum=[1];
    case 358
        extNum=[1];
    case 359
        extNum=[1];
    otherwise
        extNum=1;
end