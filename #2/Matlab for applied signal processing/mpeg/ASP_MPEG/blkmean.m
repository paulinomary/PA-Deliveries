function Ed = blkmean(D,s)

Ed = zeros(s);
for i=1:s(1),
    for j=1:s(2),
        Ed(i,j) = mean(mean(D(i:s(1):end,j:s(2):end)));
    end;
end;
