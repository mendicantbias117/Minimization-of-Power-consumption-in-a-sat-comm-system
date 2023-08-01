clc; clear; close;

% Optimization options
     options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'Hessian', 'bfgs', ...
    'ste'
    'MaxIterations', 500, ...
    'Display', 'iter',...
    PlotFcn= 'optimplotfval');

% Run the optimization
x0 = [1 1 1 1 1 1]';
[x, fval, exitflag, output, lambda, grad, hessian] = fmincon(@power_consumption, x0, [], [], [], [], [], [], @constraints, options);
return
    function [f, df] = power_consumption(x)
        % x = [Pt, R, G_loss, P_min, eta, alpha]
        Pt = x(1);
        R = x(2);
        G_loss = x(3);
        P_min = x(4);
        BW = 1e6; % bandwidth (Hz)
        S = 1e-6; % signal power (W)
        N = 1e-9; % noise power (W)
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
        f = Pt-P_min;

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
   
   
% Evaluate the objective function and constraints at the optimal point
%[f, df] = power_consumption(x);
%[c, ceq, dc, dceq] = constraints(x);
% Define the ranges of Pt and R
%Pt_range = linspace(1e-3, 1e-1, 100);
%R_range = linspace(1e6, 1e8, 100);

% Evaluate the objective function at each combination of Pt and R
%[X,Y] = meshgrid(Pt_range, R_range);
%Z = zeros(length(Pt_range), length(R_range));
%for i = 1:length(Pt_range)
 %   for j = 1:length(R_range)
  %      Z(i,j) = power_consumption([Pt_range(i), R_range(j), x(3:end)]);
  %  end
%end

% Plot the contour plot
%contourf(X, Y, Z');
%xlabel('Transmit power (W)');
%ylabel('Data rate (bps)');
%colorbar;


