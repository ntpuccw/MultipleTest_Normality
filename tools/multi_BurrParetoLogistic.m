function U = multi_BurrParetoLogistic(p, n, marginal)

    if marginal == "normal"
        a = 1;
        Yi = exprnd(1, n, p);
        X = repmat(gamrnd(a, 1, n, 1),1,p);
%         X = gamrnd(a, 1, n, p);
        U = norminv((1+Yi./X).^(-a)); % p.167
    end
    
    
    