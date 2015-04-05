% Data: a m-by-n matrix, m is the number of genes, n is the number of
% samples. survival time (Times) 
% d: censor information, 1 is uncensored, 0 is censored.
% S: normalized network information

function [bta] = NetworkRegularisedCox(Data, lamda, alpha, d, S)

% m: number of genes
% n: number of samples
[m,n] = size(Data);

%maximum interation for the outer loop
maxIter1 = 1000;

%maximum iteration for the inner loop (most of time it less than 10 iterations)
maxIter2 = 1000;

gamma = zeros(length(lamda),length(alpha),n);
bta = zeros(length(lamda),length(alpha),m);


 for nn = 1:length(alpha)
     
     % laplacian matrix 
     L = (1-alpha(nn))*(eye(m) - S) + alpha(nn)*eye(m);
     
     Z = Data'/L*Data;
     
     for mm = 1:length(lamda)
        
        gamma_outer = rand(n,1);
        lp = zeros(maxIter1,1);
        pena = zeros(maxIter1,1);
        
        C = Z*gamma_outer;

        M = zeros(n,n);
        for j = 1:n
            M(j,:) = C' - C(j); 
        end
        M = exp(M);

        tmp = zeros(n,1);
        for k = 1:n
            IX = find(d(1:k)==1);
            for p = 1:length(IX)
                tmp(k) = tmp(k) + 1/sum(M(k,IX(p):n));
            end
            %clear IX
        end

        lnh = zeros(n,1);
        for p = 1:n-1
            lnh(p) = -C(p) - log(1 + sum(M(p,p+1:n)));
        end
        lnh(n) = -C(n);

        lp_old = 0;

        for i = 1:maxIter1
            gamma_outer_old = gamma_outer;
    
            pena(i) = (1/2)*lamda(mm)*(gamma_outer'*Z*gamma_outer);
            lp(i) = sum(-tmp + d.*(lnh + C)) - pena(i);
     
            gamma_inner = gamma_outer;
    	    lp_inner_old = 0;
            
            %Newton-Raphson
            for t = 1:maxIter2
                gamma_inner_old = gamma_inner;
                dta = d -tmp;
                firstD = Z*dta - lamda(mm)*Z*gamma_inner;
                D = diag(tmp);
        
                secondD = -Z*D*Z'- lamda(mm)*Z;
                gamma_inner = gamma_inner - secondD\firstD;
                tmp = tmp.*exp(Z*gamma_inner-Z*gamma_inner_old);
                
                pena_inner = (1/2)*lamda(mm)*(gamma_inner'*Z*gamma_inner);
                lp_inner = sum(-tmp + d.*(lnh + C)) - pena_inner;
        
                if (max(abs(gamma_inner-gamma_inner_old)) < 1e-8||abs(lp_inner-lp_inner_old)<1e-8)
                    break;
                end

                lp_inner_old = lp_inner;

                if t == 1000
                    disp('The Newton-Raphson is not converge')
                end
        
            end
            gamma_outer = gamma_inner;
    
            if (max(abs(gamma_outer-gamma_outer_old)) < 1e-6||abs(lp(i)-lp_old)<1e-6)
                break;
            end
            %update the baseline hazard
    	    lp_old = lp(i);

            C = Z*gamma_outer;
    
            M = zeros(n,n);
            for j = 1:n
                M(j,:) = C' - C(j); 
            end
            M = exp(M);

            tmp = zeros(n,1);
            for k = 1:n
                IX = find(d(1:k)==1);
                for p = 1:length(IX)
                    tmp(k) = tmp(k) + 1/sum(M(k,IX(p):n));
                end
                %clear IX
            end
            lnh = zeros(n,1);
            for p = 1:n-1
                lnh(p) = -C(p) - log(1 + sum(M(p,p+1:n)));
            end
            lnh(n) = -C(n);
        end
        
        gamma(mm,nn,:) = gamma_outer;
        bta(mm,nn,:) = L\Data*gamma_outer;
        
        
     end
     %clear L Z
end
