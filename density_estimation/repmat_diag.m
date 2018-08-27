function BigA=repmat_diag(A,num)
ACell = repmat({A}, 1, num);
BigA = blkdiag(ACell{:});
end