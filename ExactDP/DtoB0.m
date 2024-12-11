function x = DtoB0(D,c)
    x = zeros(1,c);
    if D == 1
        return;
    end
    D = D - 2;
    for i = 1:length(x)
        x(i) = mod(D,2);
        D = floor(D / 2);
    end
end