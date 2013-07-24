function [flankerType flankerSessions]=getFlankerType(animal,area)

if strcmp(animal,'jack')
    if strncmp(area,'v1_2',4)
        flankerType=1:3;
        flankerSessions=[{78:93} {94:115} {116:119}];
    elseif strcmp(area,'v1_4')
        flankerType=1:3;
        flankerSessions=[{140:162} {163:184} {185:189}];
    end
elseif strcmp(animal,'blanco')
    if strncmp(area,'v1_2',4)
        flankerType=1:3;
        flankerSessions=[{[388:401 403:422]} {[431:443 451 452]} {453:459}];
    end
end