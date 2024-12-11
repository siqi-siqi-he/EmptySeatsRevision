function D = BtoD0(x)
    if any(x < 0)
        D = 1;
    else
        D = 2;
        for i = 1:length(x)
            D = D + x(i) * 2^(i-1);
        end
    end
end