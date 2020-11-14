function ranMat = diagonallyDominantMaker(n)

% generate a diagonally dominant matrix: 
ranMat = rand(n);

for i = 1:n
    sum1 = sum(ranMat(:,i));
    sum2 = sum(ranMat(i,:));
    
    % check to see which one is more: 
    minimumDiagVal = max(sum1, sum2);
    
    added = rand;
    minimumDiagVal = added + minimumDiagVal;
    
    ranMat(i,i) = minimumDiagVal;
end

end
