function [roc,vec1,vec2]= sglroc3_xing(arr1,arr2) 
%Modified from Alex's function from 13/02/97		   
    
maximum=[0 0];
maximum(1,1)=max(arr1);
maximum(1,2)=max(arr2);
maximum=max(maximum);
vec1=zeros(1,ceil(maximum)); 
vec2=zeros(1,ceil(maximum));
for i=1:ceil(maximum)
   vec1(1,i+1)=size(find(arr1<i),2)/size(arr1,2);
   vec2(1,i+1)=size(find(arr2<i),2)/size(arr2,2);
end;
roc=0;  
for i=1:ceil(maximum)-1    
   y=(vec2(i)+vec2(i+1))/2;    
   x=(vec1(i+1)-vec1(i));    
   roc=roc+(x*y); 
end  