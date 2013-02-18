ses=[343:356 358:359];
ses=[307 308 , 311 , 313 314 , 317 318 , 320 321 , 329:342];
missin=[];present=0;
for i=1:length(ses)
    for j=1:size(rocValsMat,1)
        if rocValsMat{j,2}==ses(i)
            present=1;
        end
    end
    if present==0
        missin=[missin ses(i)];
    end
    present=0;
end