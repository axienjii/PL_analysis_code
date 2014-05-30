function[time_Arr, value_Arr, SampFreq]  = read_nlx_CSC (name)


fid = fopen(name, 'rb');
if fid>0
    [head, count]=fread(fid, 16384, 'int8');
    
    % for j=1:50
    %      head(j)
    % end
    %
    % pause
    position1 = ftell(fid);
    stat1=fseek(fid,0,'eof');
    position2 = ftell(fid);
    
    %[A,count]=fscanf (fid,'%c', inf);
    numevents=(position2-position1)/1044;
    
    
    fclose(fid);
    fid = fopen(name, 'rb');
    [head, count]=fread(fid, 16384, 'int8');
    time_Arr=zeros(512,ceil(numevents));
    value_Arr=zeros(512,ceil(numevents));
    SampFreq=0;
    j=1;
    prl=[2:2:512];
    while ( ~feof (fid))
        timestamp=fread(fid, 1, 'int64');
        if (~feof (fid))
            ScNumber=fread(fid, 1, 'int32');
            SampFreq=fread(fid, 1, 'int32');
            NumValidSamp=fread(fid, 1, 'int32');
            Samp=fread(fid, 512, 'int16');
            prl=Samp;
            if NumValidSamp<512
                Samp(:)=-1;
            end
            us_per_value=[0:1:511]*(1000/SampFreq)*1000;
            time_Arr(:,j)=(us_per_value+timestamp)';
            %         size(Samp)
            %         size(value_Arr)
            value_Arr(:,j)=Samp;
            
            j=j+1;
            
        end;
    end;
    
    % prl=[1:2:size(time_Arr,1)];
    % SampFreq=round(SampFreq/2);
    % time_Arr=time_Arr(prl,:);
    % value_Arr=value_Arr(prl,:);
    time_Arr=time_Arr(:);
    value_Arr=value_Arr(:);
    fclose(fid);
else
    time_Arr=0;
    value_Arr=0;
    SampFreq=0;
end



