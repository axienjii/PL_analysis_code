function [sampleContrasts testContrasts]=area_metadata(area)

sampleContrasts=30;
if strcmp(area,'v4_0_1')
    testContrasts=[10 15 20 25 35 40 60 90];
elseif strcmp(area,'v4_0_2')
    testContrasts=[10 15 20 25 27 29 31 33 35 40 50 60];
elseif strcmp(area,'v4')||strcmp(area,'v4_1')||strcmp(area,'v4_2')||strcmp(area,'v4_0_3')
    testContrasts=[10 15 20 25 27 28 29 31 32 33 35 40 50 60];
elseif strcmp(area,'v1')||strcmp(area,'v1_1')
    testContrasts=[5 10 15 20 22 25 28 32 35 40 45 50 60 90];
elseif strncmp(area,'v1_2',4)||strncmp(area,'v1_4',4)
    sampleContrasts=[20 30 40];
    testContrasts=[5 10 12 15 18 22 25 28 35 45 60 90;
        5 10 15 22 25 28 32 35 38 45 60 90;
        5 10 15 25 32 35 38 42 45 50 60 90];
end