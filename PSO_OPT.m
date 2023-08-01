clc;clear;close all;

%optimization options
%{
hybridopts = optimoptions('fminunc', 'Display','iter',...
    'Algorithm','quasi-newton');

options = optimoptions(@particleswarm,...
    'MaxIterations',50000,...
   'Display','iter',...
    PlotFcn='pswplotbestf');
%}
hybridopts = optimoptions('fmincon','OptimalityTolerance',1e-1000);
options = optimoptions('ga','HybridFcn',{'fmincon',hybridopts},...
    'MaxStallGenerations',200,'MaxGenerations',10000,'Display','iter',...
    'PlotFcn','gaplotbestf');
%[x,fval] = ga(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,options)

% Define lower and upper bounds for the variables
lb = [];
ub = [];
x0=6;
[x,fval,exitflag,output] = ga(@power_consumption,x0,[],[],[],[],lb,ub,@constraints,options);  
return

 function [f, df] = power_consumption(x)
        % x = [Pt, R, G_loss, P_min, eta, alpha]
        Pt = x(1);
        R = x(2);
        G_loss = x(3);
        P_min = x(4);
        BW = 1e6; % bandwidth (Hz)
        S = 1e6; % signal power (W)
        N = 1e9; % noise power (W)
        R_min = BW / (log2(1 + S/N)); % minimum data rate (bps)
        
        eta = x(5);
        alpha = x(6);

        % Constants
        c = 299792458; % speed of light (m/s)
        lambda = c / (R * log2(1 + S/N)); % wavelength (m)
        D = lambda / alpha; % antenna diameter (m)
        fc = c / lambda; % carrier frequency (Hz)
        G = 10*log10((eta*pi*D^2)/lambda^2) + G_loss; % antenna gain (dBi)

        % Objective function
        f = Pt;

        % Gradient of the objective function
        df = [1; 0; 0; 0; 0; 0]; % gradient of Pt w.r.t. x
    end

    function [c, ceq, dc, dceq] = constraints(x)
        % x = [Pt, R, G_loss, P_min, eta, alpha]
        Pt = x(1);
        R = x(2);
        G_loss = x(3);
        P_min = x(4);
        eta = x(5);
        alpha = x(6);
        c = 299792458; % speed of light (m/s)
         
        % Data rate constraint
        BW = 1e6; % bandwidth (Hz)
        S = 1e-6; % signal power (W)
        N = 1e-9; % noise power (W)
        R_min = BW / (log2(1 + S/N)); % minimum data rate (bps)
        lambda = c / (R * log2(1 + S/N)); % wavelength (m)
        fc = c / lambda; % carrier frequency (Hz)
        D = lambda / alpha; % antenna diameter (m)
        G = 10*log10((eta*pi*D^2)/lambda^2) + G_loss;
        ceq(1) = R - R_min; % constraint
        
        % Antenna gain constraint
        G_min = 10*log10((0.55*pi*D^2)/lambda^2); % minimum antenna gain (dBi)
        ceq(2) = G - G_min; % constraint

        % Transmit power constraint
        c(1) = P_min - Pt; % constraint

        % Antenna diameter constraint
        ceq(3) = D - lambda/alpha; % constraint

        % Carrier frequency constraint
        ceq(4) = fc - c/lambda; % constraint

    end
