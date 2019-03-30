%% 
% This will solve a Krueger and Uhlig (2006) style model with idiosyncratic
% income uncertainty and one-sided commitment. 

%% Parameters and Housekeeping

clear 
clc
tic
V.par.sigma             = 2;
V.par.gamma             = 0.975;
V.par.markov            = [V.par.gamma,1-V.par.gamma;1-V.par.gamma,V.par.gamma];
V.par.end               = [0.7,1.3];
V.par.beta              = 0.96;
V.par.AvgEarn           = sum(MarkovChains.StationaryDistribution(V.par.markov).*V.par.end);
if V.par.sigma == 1
    V.fun.U             = @(c) log(c);
else
    V.fun.U             = @(c) (c.^(1-V.par.sigma)-1)./(1-V.par.sigma);
end
autV                    = solve(@(L,H) V.fun.U(V.par.end)'-[L;H]+V.par.beta*sum(V.par.markov.*[L,H],2));
autV                    = double([autV.L,autV.H]);

%% Outer Bisection - R

R_bisect                    = [1,1/V.par.beta];
thr_outer                   = [Inf,Inf];
while thr_outer(1) > eps && thr_outer(2) > eps
    
    %% Inner Bisection - Cn
    
    R                       = mean(R_bisect);
    C_bisect                = V.par.end;
    thr_inner               = [Inf,Inf];
    while thr_inner(1) > eps && thr_inner(2) > eps
        
        % To find n...
        Cn                  = mean(C_bisect);
        c                   = (((V.par.beta*R).^(1/V.par.sigma)).^(0:200)).*Cn;
        c                   = fliplr([c(1:(find(c<V.par.end(1)==1,1)-1)),V.par.end(1)])';
        Us                  = V.fun.U(c);
        v                   = zeros(length(c),1);
        
        % To find vn...
        for k = 1:length(c) 
            if k == 1
                v(k,:)          = autV(1);
            elseif k == length(c)
                v(k,:)          = Us(k) + V.par.beta*(V.par.gamma*autV(2)+(1-V.par.gamma)*v(k-1));
            else
                v(k,:)          = Us(k)+V.par.beta*(V.par.gamma*v(k-1)+(1-V.par.gamma)*autV(2));
            end
        end
        
        % Bisection step...
        F = v(end) - autV(2);
        if F > 0; C_bisect(2)     = Cn;
        elseif F < 0 ; C_bisect(1)     = Cn;
        end
        thr_inner           = double([norm(Cn-mean(C_bisect)),norm(F)]);
    end
    
    PTM                     = zeros(length(c),length(c));
    PTM(1,1)                = V.par.gamma;
    PTM(end,end)            = V.par.gamma;
    PTM(1:(end-1),end)      = 1-V.par.gamma;
    for j = 2:length(c)
        PTM(j,j-1)          =   V.par.gamma;
        if j == length(c)
            PTM(j,j-1)          = 1-V.par.gamma;
        end
    end
    
    % Market clearing...
    stdis                   = MarkovChains.StationaryDistribution(PTM);%((ones(1,length(c))/length(c))*PTM^1000)';
    AvgEarn                 = sum(c.*stdis);
    F                       = AvgEarn - V.par.AvgEarn;
    
    %Bisection step...
    if F > 0 
        R_bisect(2) = R;
    elseif F < 0
        R_bisect(1) = R;
    end
    thr_outer           = double([norm(R-mean(R_bisect)),norm(AvgEarn-1)]);
end

%% Pretty Results

out     = array2table([c,MarkovChains.StationaryDistribution(PTM)],'VariableNames',{'ConsumptionLevel','StationaryProbability'},'RowNames',strseq('C',[1:length(c)]));
disp(out)
plot(c,'LineWidth',3)
xlim([1,length(c)]); xlabel('Cycle'); ylabel('Consumption'); title('Consumption Levels'); ylim(V.par.end);
disp('Market-clearing interest rate:'); disp('=============================='); disp(R);
toc