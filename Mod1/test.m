y = y_max_l;%[1.0, 0.8, 0.64, 0.512, 0.4096];  % example data
% t = (0:length(y)-1) * dt;  % create time vector
t = t_max_l;
dt = t(2)-t(1);  % time step
ln_y = log(y);  % natural log of data

p = polyfit(t, ln_y, 1);  % linear fit: ln(y) = -alpha * t + ln(A)
alpha = -p(1);  % negative slope is the damping factor

printf("Estimated damping factor (alpha): %f\n", alpha);

% Example estimated parameters
A = 1;           % initial amplitude (you can adjust this)
N = 20;          % number of points to plot

t = (0:N-1) * dt;         % time vector
y_exp = A * exp(-alpha * t);  % exponential decay

% Plot
figure;
plot(t, y_exp, 'r-', 'LineWidth', 2);
xlabel('Time');
ylabel('Amplitude');
title('Exponential Decay: y(t) = Ae^{-\alpha t}');
grid on;
hold on;
y_original = y_max_l;%[1.0, 0.8, 0.64, 0.512, 0.4096];  % your actual data
plot((0:length(y_original)-1)*dt, y_original, 'bo-', 'LineWidth', 2);
legend('Estimated Exponential', 'Original Data');