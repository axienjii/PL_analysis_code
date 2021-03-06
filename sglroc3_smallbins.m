function [roc,vec1,vec2]= sglroc3_smallbins(arr1,arr2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	 
%%%%%%%%	Make ROC from normalized (to mean) data	%%%%%%%% %			
%AUTH: AT %			
%VERS: 1.0 %			
%DATE: 13/02/97 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%Modified by Xing to use smaller bin divisions
binMultiple=100;
sum_arr=size(arr1,2)+size(arr2,2);    
maximum=[0 0];
size_arr=[0 0];
size_arr(1,1)=size(arr1,2);
size_arr(1,2)=size(arr2,2);
max_size=max(size_arr);
maximum(1,1)=max(arr1);
maximum(1,2)=max(arr2);
maximum=max(maximum);
steps=binMultiple*maximum/(sum_arr/2);   
vec1=zeros(1,binMultiple*ceil(maximum+2)); 
vec2=zeros(1,binMultiple*ceil(maximum+2));
count=0; 
for i=1:binMultiple*maximum+2
   vec1(1,i)=size(find(arr1>maximum-(i-1)/binMultiple),2)/size(arr1,2);
   vec2(1,i)=size(find(arr2>maximum-(i-1)/binMultiple),2)/size(arr2,2);
   %(i-1)*steps
end;
choice_prob=0;  
for i=2:binMultiple*maximum+2    
   y=1-(vec1(i-1)+vec1(i))/2;    
   x=(vec2(i)-vec2(i-1));    
   choice_prob	= choice_prob+(x*y); 
end  
roc=1-choice_prob;
  	   