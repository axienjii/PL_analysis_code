function remove_empty_333(session)
%remove empty cells from session 333 and 311
animal='blanco';
area='v4_1';
channels=main_channels(animal,area);
sessions=[311 333];
minTrials=[50 60];
for sessInd=1:length(sessions)
    for i=1:length(channels)
        loadText=['load F:\PL\spikeData\blanco\',num2str(channels(i)),'_',num2str(sessions(sessInd)),'_30.mat'];
        eval(loadText);
        for condInd=1:14
            emptyRows=[];
            for rowInd=1:size(matarray{condInd,2},1)
                if isempty(matarray{condInd,2}{rowInd})
                    emptyRows=[emptyRows rowInd];
                end
            end
            firstEmpty=emptyRows(diff(emptyRows)==1);%first empty cell, out of sequence of at least 2 consequtive indices
            firstEmpty=firstEmpty(firstEmpty>minTrials(sessInd));
            for columnInd=1:size(matarray,2)
                matarray{condInd,columnInd}=matarray{condInd,columnInd}(1:firstEmpty(1)-1);
            end
        end
        saveText=['save F:\PL\spikeData\blanco\',num2str(channels(i)),'_',num2str(sessions(sessInd)),'_30.mat matarray'];
        eval(saveText);
    end
end