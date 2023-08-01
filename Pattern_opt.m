clc;clear;close all;

%optimization options
options = optimoptions('patternsearch',...
    'Algorithm','classic',...
    'Display','iter',...
    PlotFcn='psplotmeshsize');

% Define lower and upper bounds for the variables
lb = [];
ub = [];
x0=[1 1 1 1 1 1];
[x,fval,exitflag,output] = patternsearch(@power_consumption,x0,...
    [],[],[],[],lb,ub,...
    @constraints,options);
return

   function [f] = power_consumption(x)
        Pt = x(1);
        % Objective function
        f = Pt;
    
    end

    function [c, ceq, dc, dceq] = constraints(x)
        % x = [Pt, R, G_loss, P_min, eta, alpha]
        Pt = x(1);
        R = x(2);
        G_loss = x(3);
        P_min = 100;
        eta = x(5);
        alpha = x(6);
        Sl = 299792458; % speed of light (m/s)
         
        % Data rate constraint
        BW = 1e6; % bandwidth (Hz)
        S = 1e-6; % signal power (W)
        N = 1e-9; % noise power (W)
        R_min = BW / (log2(1 + S/N)); % minimum data rate (bps)
        lambda = Sl / (R * log2(1 + S/N)); % wavelength (m)
        fc = Sl / lambda; % carrier frequency (Hz)
        D = lambda / alpha; % antenna diameter (m)
        G = 10*log10((eta*pi*D^2)/lambda^2) + G_loss;
       
        c(3) = R - R_min; % constraint
        
        % Antenna gain constraint
        G_min = 10*log10((0.55*pi*D^2)/lambda^2); % minimum antenna gain (dBi)
        c(2) = G - G_min; % constraint

        % Transmit power constraint
        c(1) = Pt - P_min; % constraint

        % Antenna diameter constraint
        ceq(1) = D - lambda/alpha; % constraint

        % Carrier frequency constraint
        ceq(2) = fc - Sl/lambda; % constraint
        

    end