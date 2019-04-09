%% Housekeeping

clear 
clc
close all
% This requires some Markov functions...
addpath('C:\Users\curta\OneDrive - The University of Western Ontario\Github\NumericalMethods\Markov Processes')

%% Exogenous Parameters

% Object Dimensions (N)
N.y                     = 200;
N.b                     = 500;

% Model Parameters (P)
P.beta                  = 0.953;
P.sigma                 = 2.0d0;
P.rf                    = 1.017;
P.rho                   = 0.945;
P.shocksd               = 0.025;
P.theta                 = 0.282;
P.yhat                  = 0.969;
P.blo                   = -0.35;
P.bhi                   = 0.2;
P.dispinc               = 50;               % How often to display convergence results..

% Grid Objects (G)
[G.y,P.markov]      = Rouwenhorst(N.y,P.shocksd^2,0.5*(1+P.rho));
G.y                 = sort(exp(G.y));
G.b                 = sort([linspace(P.blo,P.bhi,N.b),0])';

% The first item here is the ad hoc default cost, the second is an
% alternative for cost that is a fixed amount...

G.h                 = G.y.*(G.y <= P.yhat )+ P.yhat.*(G.y > P.yhat);
%P.deltaaut          = 0.995;
%G.h                 = G.y*P.deltaaut;
N.b                 = size(G.b,1);
P.zerob             = find(G.b == 0);

% Functions (F)
F.u                 = @(c) (c.^(1-P.sigma))./(1-P.sigma);
F.cr                = @(y,q,bp) y + G.b - (q.*bp)';
F.cd                = G.h;
F.evr               = @(vo) P.markov*vo';
F.evd               = @(vo,vd) P.theta*P.markov*vo(P.zerob,:)' + (1-P.theta)*P.markov*vd;

%% VFI
itr                 = 1;
thr                 = [Inf,Inf];
timerr              = cputime;
while max(thr) > 1e-04

    if itr == 1
        tic
    end
    
    % For repayers...
    if itr == 1
        C                   = G.y'+G.b;
        v.r1                = F.u(G.y'+G.b);
        v.r1(C<=0)           = -realmax;
    else
        EvR                 = F.evr(v.o0)';
        for i = 1:N.y
            C               = F.cr(G.y(i),v.q(:,i),G.b);
            U               = F.u(C);
            U(C<0)          = -realmax;
            [v.r1(:,i),...
                po.r1(:,i)] = max(U+P.beta*EvR(:,i)',[],2);
        end
        thr(1)              = norm(v.r0-v.r1);
        clear U C EvR
    end
    
    % For defaulters...
    C                   = G.h;
    U                   = F.u(C);
    U(C<0)              = -realmax;
    if itr == 1
        v.d1            = U;
        clear U C EvD
    else
        EvD             = F.evd(v.o0,v.d0);
        for i = 1:N.y
            v.d1        = U+P.beta*EvD;
        end
        thr(2)          = norm(v.d0-v.d1);
            clear U C EvD
    end
    
    % Options, decision rules, and prices...
    for y = 1:N.y
        for b = 1:N.b
            v.o1(b,y)           = max(v.r1(b,y),v.d1(y));
            v.R(b,y)            = v.r1(b,y) > v.d1(y);
        end
    end
    
    % Updating...
    v.q                         = (P.markov*v.R'/P.rf)';
    v.r0                        = v.r1;
    v.d0                        = v.d1;
    v.o0                        = v.o1;
    if mod(itr,P.dispinc) == 0
        disp(['[Itr. #',sprintf('%04d',itr),' in ',num2str(round(toc)),' second(s)]  [Repay Convergence: (',sprintf('%12.11f',thr(1)),')] [Default Convergence: (',sprintf('%12.11f',thr(2)),')]'])
        tic
    end
    itr                         = itr + 1;
    clear y b 
end
disp(['[Itr. #',sprintf('%04d',itr),' in ',num2str(round(toc)),' second(s)]  [Repay Convergence: (',sprintf('%12.11f',thr(1)),')] [Default Convergence: (',sprintf('%12.11f',thr(2)),')]'])

%% Plots

% Generate moments of the Markov process and find the index of relevant
% quantities...
% P.moments                       = MarkovMoments(P.markov,G.y);
% P.meanY                         = P.moments.Mean;
% P.sdY                           = (P.moments.Variance)^(1/2);
% [~,meannum]                     = min(abs(G.y-P.meanY));
% [~,sdhinum]                     = min(abs(G.y-P.meanY-0.05));%P.sdY));
% [~,sdlonum]                     = min(abs(G.y-P.meanY+0.05));%P.sdY));

% % Bond Price Schedule by Shock...
% figure('DefaultAxesFontSize',15)
% hold on
% Blev    = [round(N.b/5),round(N.b/5)*2,round(N.b/5)*3];
% Bval    = round(G.b(Blev),3);
% plot(G.y,v.q(Blev(1),:),'LineWidth',4)
% plot(G.y,v.q(Blev(2),:),'LineWidth',4)
% plot(G.y,v.q(Blev(3),:),'LineWidth',4)
% legend(strcat('B''=',num2str(Bval(1))),strcat('B''=',num2str(Bval(2))),strcat('B''=',num2str(Bval(3))))
% xlabel('Income Shock $y$','Interpreter','latex')
% ylabel('Bond Price $Q(B'',Y)$','Interpreter','latex')
% 
% % Probability of Default by B'...
% figure('DefaultAxesFontSize',15)
% hold on
% plot(G.b,1-P.rf*v.q(:,meannum),'LineWidth',4)
% plot(G.b,1-P.rf*v.q(:,sdhinum),'LineWidth',4)
% plot(G.b,1-P.rf*v.q(:,sdlonum),'LineWidth',4)
% legend('Mean Shock','High Shock','Low Shock')
% xlabel('New Debt $B''$','Interpreter','latex')
% ylabel('Probability of Default $\mu(B'',Y)$','Interpreter','latex')

% Policy Functions - Figure 4
% figure('DefaultAxesFontSize',15)
% hold on
% plot(G.b,G.b(po.r1(:,meannum)),'LineWidth',4)
% plot(G.b,G.b(po.r1(:,sdhinum)),'LineWidth',4)
% plot(G.b,G.b(po.r1(:,sdlonum)),'LineWidth',4)
% refline(1,0)
% xlabel('$B$','Interpreter','latex')
% ylabel('$B''$','Interpreter','latex')
% title('Savings Function $B''(B,y)$, $\delta = 0.995$','Interpreter','latex')
% legend('Mean Shock','1SD Up Shock','1SD Down Shock','Location','southeast')
% xlim([P.blo,0]); ylim([P.blo,0])

% Value Functions - Figure 4
% figure('DefaultAxesFontSize',15)
% hold on
% plot(G.b,repelem(v.d1(meannum),N.b,1),'r--','LineWidth',2)
% plot(G.b,v.r1(:,meannum),'r-','LineWidth',2)
% plot(G.b,repelem(v.d1(sdhinum),N.b,1),'b--','LineWidth',2)
% plot(G.b,v.r1(:,sdhinum),'b-','LineWidth',2)
% xlim([P.blo,P.bhi])
% xlabel('$B$','Interpreter','latex')
% ylabel('Utility','Interpreter','latex')
% title('Value Functions $V^d(y)$ and $V^{nd}(B,y)$, $\delta = 0.995$','Interpreter','latex')
% legend('$V^{d}(B,y^{mean})$','$V^{nd}(y^{mean})$','$V^{d}(B,y^{hi})$','$V^{nd}(y^{hi})$','Location','southeast','Interpreter','latex')

% Bond Price Schedule - Figure 3
% figure('DefaultAxesFontSize',15)
% hold on
% plot(G.b,v.q(:,meannum),'LineWidth',4)
% plot(G.b,v.q(:,sdhinum),'LineWidth',4)
% plot(G.b,v.q(:,sdlonum),'LineWidth',4)
% xlabel('$B''$','Interpreter','latex')
% ylabel('$Q(B'',A)$','Interpreter','latex')
% title('Bond Price Schedule $Q(B'',y)$, $\delta = 0.995$','Interpreter','latex')
% xlim([P.blo,0])
% legend('Mean Shock','1SD Up Shock','1SD Down Shock','Location','northwest')

% Simulation...
% clear D B BI C
% P.T                     = 1000000;
% StateShocks             = MarkovSimulate(P.markov,P.T);
% B(1)                    = 0;
% BI(1)                   = P.zerob;
% 
% for t = 1:P.T
%     if t == 1
%         % Default Check...
%         D(t)            = 1 - v.R(BI(t),StateShocks(t));
%         % If better to default...
%         if D(t) == 1
%             BI(t+1)     = NaN;
%             B(t+1)      = NaN;
%             C(t)        = G.h(StateShocks(t));
%         % If better to repay...
%         else
%             BI(t+1)     = po.r1(BI(t),StateShocks(t));
%             B(t+1)      = G.b(BI(t+1));
%             C(t)        = G.y(StateShocks(t)) + B(t) - B(t+1)*v.q(BI(t),StateShocks(t));
%         end
%     else
%         % Default last period check...
%         if D(t-1) == 0
%             % Default Check...
%             D(t)            = 1 - v.R(BI(t),StateShocks(t));
%             % If better to default...
%             if D(t) == 1
%                 BI(t+1)     = NaN;
%                 B(t+1)      = NaN;
%                 C(t)        = G.h(StateShocks(t));
%             % If better to repay...
%             else
%                 BI(t+1)     = po.r1(BI(t),StateShocks(t));
%                 B(t+1)      = G.b(BI(t+1));
%                 C(t)        = G.y(StateShocks(t)) + B(t) - B(t+1)*v.q(BI(t),StateShocks(t));
%             end
%         else
%             % Draw probability of reentering
%             Reenter = rand;
%             % Reenter
%             if Reenter <= P.theta
%                 D(t)        = 0;
%                 BI(t+1)     = P.zerob;
%                 B(t+1)      = 0;
%                 C(t)        = G.y(StateShocks(t));
%             % Don't reenter
%             else
%                 D(t)        = NaN;
%                 BI(t+1)     = NaN;
%                 B(t+1)      = NaN;
%                 C(t)        = G.h(StateShocks(t));
%             end
%         end
%     end
% end
% 
% Results = array2table([C',D',B(1:(end-1))',G.y(StateShocks)],'VariableNames',{'C','D','B','Y'});
% % Simulation Results
% figure('DefaultAxesFontSize',15)
% scatter3(Results.Y(Results.B < -0.001),Results.B(Results.B < -0.001),Results.D(Results.B < -0.001))
% xlabel('Income Shock'); ylabel('Debt'); 
% zticks([0,1]); set(gca,'ztick',[0,1],'zticklabel',{'Repay','Default'})

