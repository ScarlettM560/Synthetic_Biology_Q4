%% Parameters
                 
q = 0.5002;              % example value (must satisfy p+q=1)
p = 1-q; 
c = 1/24;                % set your desired time constant
r = q/p;                 % r > 1 by assumption
N = 1000;                % number of points4
%====================================================
%% Calcualting the time to extinction from N
%====================================================
% First use recurrsion
T = extinction(p,q,c,r,N);

%% Now try the analytic solution
% Note the analytic solution is no quicker if you are calculating all
% soltions up to n. If is quicker if you are calculating the time to
% extinction for a 'one-off' very large value of n.
T_an  = extinction_an(p,q,c,r,N);


T_noise  = extinction_noise(p,q,c,r,N, 0.0005);

%% Plot a comparison of the two method for calcualting time to extinction
close all
%p1 = compare_n(T,T_an,T_noise, r,N);
p1 = compare(T,T_an,r,N);

%
%====================================================
%% Now calculate the time for going from 1 to N 
%====================================================
%% Again, start with the recurrence method
%
M = time_to_i(p,q,c,r,N);
%% Now try the analytic solution
M_an = time_to_i_an(p,q,c,r,N);


%% Plot a comparison of the two method for calcualting time to n
p2 = compare2(M,M_an,r,N)

%% Time to go extinct when p = q

q = 0.500;  
p = 1-q; 
c = 1/24;               % set your desired constant
r = q/p;                % r > 1 by assumption
N = 1000;               %number of points 
n1 = 10^16                % carrying capacity

%====================================================
%% Calcualting the time to extinction from N with carry capacity n
%====================================================
%
T1 = unbiased_extinction(p,r,c,N,n1);
T1_an = unbiased_extinction_an(p,r,c,N,n1);
n2 = 10^4
T2 = unbiased_extinction(p,r,c,N,n2);
T2_an = unbiased_extinction_an(p,r,c,N,n2);
p1 = compare3(T1,T1_an,T2, T2_an, r,N, n1, n2);


%%
function T = extinction(p,q,c,r,N)
%====================================================
%% Calcualting the time to extinction from N
%====================================================

%% Use the recurrence relationship first
% Preallocate
Delta = zeros(N,1);
T = zeros(N,1);

%% Initial value of Delta_1 ensuring T_i stays finite
Delta(1) = -(c/p) * log(1 - 1/r);

%% Recurrence for Delta_i
for i = 1:N-1
    Delta(i+1) = r*Delta(i) - c/(p*i);
end

%% Build T_i = sum of deltas
T = cumsum(Delta);
end


%%
function T = extinction_noise(p,q,c,r,N,noise)
%====================================================
%% Calcualting the time to extinction from N
%====================================================

%% Use the recurrence relationship first
% Preallocate
Delta = zeros(N,1);
T = zeros(N,1);
np=(p+(0.5-rand())*noise);
r=(1-np)/np;

%% Initial value of Delta_1 ensuring T_i stays finite
Delta(1) = -(c/p) * log(1 - 1/r);

%% Recurrence for Delta_i
for i = 1:N-1
    Delta(i+1) = r*Delta(i) - c/(np*i);
end

%% Build T_i = sum of deltas
T = cumsum(Delta);
end



%% Analtytic solution of time to extinction r> 1
function T_an = extinction_an(p,q,c,r,N)
T_an = zeros(N,1);
for k = 2:N
    term1 = integral(@(t) (1-t.^(k-1))./(1-t),0,1/r);
    term2 = log(k-1)+0.57721 + 1./(2*(k-1)) -1./(12*(k-1).^3);
    term3 = (r.^k -1)*log(r/(r-1));
    T_an(k) = -c/(p*(r-1))*((r.^k)*term1 - term2 -term3);
end
end

function p1 = compare(T,T_an,r,N)
figure('Units','centimeters','Position',[2 2 16 12]);
hold on
plot(1:N,log10(T_an),'r', 'LineWidth', 1.6);
plot(1:N, log10(T), 'k--', 'LineWidth', 1.6);  


xlabel('\it{n}', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Log(Time) (days)', 'FontSize', 14, 'FontWeight', 'bold');
s = "Time to extinction \it {q/p} = " + r;
title(s, 'FontSize', 15, 'FontWeight', 'bold');

grid on;
set(gca, ...
    'FontSize', 13, ...
    'LineWidth', 1.2, ...
    'Box', 'on', ...
    'TickDir', 'out');
legend({'Analytic','Recurrence' }, ...
    'Location', 'northwest', ...
    'FontSize', 12, ...
    'Box', 'off');

% Improve figure appearance
set(gcf, 'Color', 'w');     % white background for publication
p1 = 1;
% Uncomment to export as vector PDF for publication:
% print(gcf, '-dpdf', 'logTi_plot.pdf');
hold off
%% 
end


function p1 = compare_n(T,T_an,T_n, r,N)
figure('Units','centimeters','Position',[2 2 16 12]);
hold on
plot(1:N,(T_an),'r', 'LineWidth', 1.6);
plot(1:N, (T), 'k--', 'LineWidth', 1.6); 
plot(1:N, (T_n), 'b', 'LineWidth', 1.6); 


xlabel('\it{n}', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Log(Time) (days)', 'FontSize', 14, 'FontWeight', 'bold');
s = "Time to extinction \it {q/p} = " + r;
title(s, 'FontSize', 15, 'FontWeight', 'bold');

grid on;
set(gca, ...
    'FontSize', 13, ...
    'LineWidth', 1.2, ...
    'Box', 'on', ...
    'TickDir', 'out');
legend({'Analytic','Recurrence', 'Noisy' }, ...
    'Location', 'northwest', ...
    'FontSize', 12, ...
    'Box', 'off');

% Improve figure appearance
set(gcf, 'Color', 'w');     % white background for publication
p1 = 1;
% Uncomment to export as vector PDF for publication:
% print(gcf, '-dpdf', 'logTi_plot.pdf');
hold off
%% 
end

function M = time_to_i(p,q,c,r,N)
%% Preallocate
Delta = zeros(N,1);
M = zeros(N,1);

%% Initial value of Delta_1 ensuring that you cant go from 0  to 1
Delta(1) = 0;%% Preallocat

%% Recurrence for Delta_i
for i = 2:N
    Delta(i) = r*Delta(i-1) + c/(p*(i-1));
end

%% Build T_i = sum of deltas
M = cumsum(Delta);
end

function M_an = time_to_i_an(p,q,c,r,N);

%%  analytic solution
% Note the analytic solution is no quicker if you are calculating all
% soltions up to n. If is quicker if you are calculating the time to
% extinction for a 'one-ff' very large value of n.
M_an = zeros(N,1);
for k = 2:N
    term1 = integral(@(t) (1-t.^(k-1))./(1-t),0,1/r);
    term2 = log(k-1)+0.57721 + 1./(2*(k-1)) -1./(12*(k-1).^3);
    
    M_an(k) = c/(p*(r-1))*((r.^k)*term1 - term2 );
end

end

function p2 = compare2(M,M_an,r,N)
figure('Units','centimeters','Position',[2 2 16 12]);


hold on
plot(1:N, log10(M), 'r-', 'LineWidth', 1.6);
plot(1:N, log10(M_an), 'k--', 'LineWidth', 1.6);


xlabel('\it{n}', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Log(Time) (days)', 'FontSize', 14, 'FontWeight', 'bold');
s = "Time to go from 1 to \it{n  q/p} = " + r;
title(s, 'FontSize', 15, 'FontWeight', 'bold');

grid on;

set(gca, ...
    'FontSize', 13, ...
    'LineWidth', 1.2, ...
    'Box', 'on', ...
    'TickDir', 'out');
legend({'Analytic','Recurrence' }, ...
    'Location', 'northwest', ...
    'FontSize', 12, ...
    'Box', 'off');

%% Improve figure appearance
set(gcf, 'Color', 'w');     % white background for publication
hold off
p2 = 2;
end

function T = unbiased_extinction(p,r,c,N,n)
%% Use the recurrence relationship first
% Preallocate
Delta = zeros(N,1);
T = zeros(N,1);
H_n = log(n-1)+0.57721 + 1./(2*(n-1)) -1./(12*(n-1).^3);


%% Initial value of Delta_1 ensuring T_i stays finite
Delta(1) = 1/N + 2*c*H_n;

%% Recurrence for Delta_i
for i = 1:N-1
    Delta(i+1) = r*Delta(i) - c/(p*i);
end

%% Build T_i = sum of deltas
T = cumsum(Delta);
end

function T = unbiased_extinction_an(p,r,c,N,n)
%% Use the recurrence relationship first
% Preallocate
NN = [2:N];
T = zeros(N,1);
H_n = log(n-1)+0.57721 + 1./(2*(n-1)) -1./(12*(n-1).^3);
H_i = log(NN-1)+0.57721 + 1./(2*(NN-1)) -1./(12*(NN-1).^3);
T(2:N) = 2*c.*NN.*(H_n-H_i)+ NN./n + 2*c.*(H_i - 1);
end

function p1 = compare3(T,T_an,T2, T2_an,r,N, n1, n2)
figure('Units','centimeters','Position',[2 2 16 12]);
hold on
plot(1:N,log10(T_an),'r', 'LineWidth', 1.6);
plot(1:N, log10(T), 'r--', 'LineWidth', 1.6);  
plot(1:N,log10(T2_an),'b', 'LineWidth', 1.6);
plot(1:N, log10(T2), 'b--', 'LineWidth', 1.6);  
% semilogy(1:N,T_an,'r', 'LineWidth', 1.6);
% semilogy(1:N, T, 'r--', 'LineWidth', 1.6);  
% semilogy(1:N,T2_an,'b', 'LineWidth', 1.6);
% semilogy(1:N, T2, 'b--', 'LineWidth', 1.6);
xlabel('\it{i}', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Log(Time) (days)', 'FontSize', 14, 'FontWeight', 'bold');
s = "Time to extinction {\it{q/p}} = " + r + " with carrying capacity n";
title(s, 'FontSize', 15, 'FontWeight', 'bold');

grid on;
set(gca, ...
    'FontSize', 13, ...
    'LineWidth', 1.2, ...
    'Box', 'on', ...
    'TickDir', 'out');
     s1 = "Analytic {\it{n}} = 10^{16}";
     s2 = "Recurrence {\it{n}} = 10^{16}";
     s3 = "Analytic {\it{n}} = 10^{4}";
     s4 = "Recurrence {\it{n}} = 10^{4}";
legend({s1,s2,s3,s4 }, ...
    'Location', 'southeast', ...
    'FontSize', 12, ...
    'Box', 'off');

% Improve figure appearance
set(gcf, 'Color', 'w');     % white background for publication
p1 = 1;
% Uncomment to export as vector PDF for publication:
% print(gcf, '-dpdf', 'logTi_plot.pdf');
hold off
%% 
end
