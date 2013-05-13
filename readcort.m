function[time_arr,event_arr,eog_arr,header,trial]  = readcort (name,read_eog)

tic;
trial=0;
name
fid = fopen(name, 'rb');
hd=zeros(1,13);
max_t_arr=0;
max_ev_arr=0;
max_eog_arr=0;
max_epp_arr=0;
while ( ~feof (fid))
   length = fread(fid, 1, 'uint16');
   if (isempty(length)~=1)
      hd(1,1:8)= (fread(fid, 8, 'uint16'))';
      hd(1,5)=hd(1,5)/4;
      if hd(1,5)>max_t_arr
          max_t_arr=hd(1,5);
      end
      if hd(1,6)>max_ev_arr
          max_ev_arr=hd(1,6);
      end
      if hd(1,7)>max_eog_arr
          max_eog_arr=hd(1,7);
      end
      if hd(1,8)>max_epp_arr
          max_epp_arr=hd(1,8);
      end

      hd(1,6)=hd(1,6)/2;
      hd(1,7)=hd(1,7)/2;
      hd(1,8)=hd(1,8)/2;
      hd(1,9:10) = (fread(fid, 2, 'uchar'))';
      hd(1,11:13)= (fread(fid, 3, 'uint16'))';
      time_arr  = (fread (fid,(hd(1,5)) , 'ulong'));
      event_arr = (fread (fid,(hd(1,6)), 'uint16'));
      eog_arr = fread (fid,(hd(1,7)), 'short');
      epp = fread (fid,(hd(1,8)), 'uint16');
      trial=trial+1;
   end; 
end;
toc
fclose(fid);
tic;
fid = fopen(name, 'rb');
time_arr =zeros(max_t_arr, trial);
event_arr=zeros(max_ev_arr, trial);
if (nargin == 2)
   eog_arr =zeros(max_eog_arr, trial)-5000;
   read_eog=1;
else
   eog_arr=[];
	read_eog=0;   
end;


header=zeros(13,trial);
trial=0;
while ( ~feof (fid))
   length = fread(fid, 1, 'uint16');
   if (isempty(length)~=1)
      hd(1,1:8)= (fread(fid, 8, 'uint16'))';
      hd(1,5)=hd(1,5)/4;
      hd(1,6)=hd(1,6)/2;
      hd(1,7)=hd(1,7)/2;
      hd(1,8)=hd(1,8)/2;
      hd(1,9:10) = (fread(fid, 2, 'uchar'))';
      hd(1,11:13)= (fread(fid, 3, 'uint16'))';
      %if (hd(1,13)==0)
      %   eog=hd(1,7)
      %   epp=hd(1,8)
      %   pause
      %end
      
      if hd(1,5)>0
         time_arr(1:(hd(1,5)),trial+1)   = (fread (fid,(hd(1,5)) , 'ulong')); 
      end;
      if hd(1,6)>0
         event_arr(1:(hd(1,6)),trial+1) = (fread (fid,(hd(1,6)), 'uint16'));
      end;
      
      %find(event_arr(:,trial+1)==1)
      if read_eog==1
          if hd(1,7)>0
              eog_arr(1:(hd(1,7)),trial+1) = fread (fid,(hd(1,7)), 'short');
          end
      else
          eog=fread (fid,(hd(1,7)), 'uint16');
      end
      epp = fread (fid,(hd(1,8)), 'uint16');
      header(1:13,trial+1)=hd(1,1:13)';
      trial=trial+1;
   end; 
end;
toc
