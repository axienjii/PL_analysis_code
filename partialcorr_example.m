function partialcorr_example
X = randn(100,4);
m = size(X,2);
R = eye(m);
for i = 1:m
    for j = 1:i-1
        k = setdiff(1:m,[i j])
        R(i,j) = partialcorr(X(:,i),X(:,j),X(:,k)); 
        R(j,i) = R(i,j)
    end
end