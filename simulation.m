clc
close all
clear all

% Projeto 2 - Comunicacoes Moveis - PPGEE 2025/2
% Autoria: Pedro Henrique Dornelas Almeida

%% Rayleigh Channel Simulation

points = 100;
N = 15;  % simulation points
samples = 1e7 ;

gamma_th_db = linspace(-30, 30, points);   % para ficar igualmente distribuido no gráfico em db
gamma_th = 10 .^ (gamma_th_db ./ 10);      % conversao para linear

gamma_th_sim_db = linspace(-30, 30, N);         % para ficar igualmente distribuido no gráfico em db
gamma_th_sim = 10 .^ (gamma_th_sim_db ./ 10);   % conversao para linear

gamma_bar_db = [-20, 0, 20];
gamma_bar = 10 .^ (gamma_bar_db ./ 10);


% Processing Outage Probability
pout = zeros(points, length(gamma_bar));
pout_sim = zeros(N, length(gamma_bar));
for i = 1 : length(gamma_bar)
    pout(:, i) = pout_rayleigh(gamma_th, gamma_bar(i));
    [pout_sim(:, i)] = pout_sim_rayleigh(gamma_th_sim, gamma_bar(i), samples);
end

% Plot Outage Probability
figure(1)
color = ['b', 'r', 'g'];
for i = 1 : length(gamma_bar)
    semilogy(gamma_th_db, pout(:, i), 'Color', color(i), 'Linewidth', 1.5)
    hold on
    semilogy(gamma_th_sim_db, pout_sim(:, i), 'x', 'Color', 'k', 'Linewidth', 1.5)
    hold on
end
grid on
ylabel('OP', 'Interpreter', 'Latex')
xlabel('$\gamma_{th}$ (dB)', 'Interpreter', 'Latex')
legend("$\bar{\gamma_{s}} = " + num2str(gamma_bar_db(1)) + "$dB",...
       "",...
       "$\bar{\gamma_{s}} = " + num2str(gamma_bar_db(2)) + "$dB",...
       "",...
       "$\bar{\gamma_{s}} = " + num2str(gamma_bar_db(3)) + "$dB",...
       "Simulated",...
       'Interpreter', 'Latex', 'Location', 'southeast')
title("Rayleigh", "Interpreter", "Latex")
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 14;


%% Rice Channel Simulation

points = 100;
N = 15;  % simulation points

gamma_th_db = linspace(-30, 30, points);   % para ficar igualmente distribuido no gráfico em db
gamma_th = 10 .^ (gamma_th_db ./ 10);      % conversao para linear

gamma_th_sim_db = linspace(-30, 30, N);         % para ficar igualmente distribuido no gráfico em db
gamma_th_sim = 10 .^ (gamma_th_sim_db ./ 10);   % conversao para linear

gamma_bar_db = [-20, 0, 20];
gamma_bar = 10 .^ (gamma_bar_db ./ 10);

kr = [0.1, 1, 10];

fig_number = 2;
for j = 1 : length(kr)
    % Processing Outage Probability
    pout = zeros(points, length(gamma_bar));
    pout_sim = zeros(N, length(gamma_bar));
    for i = 1 : length(gamma_bar)
        pout(:, i) = pout_rice(gamma_th, gamma_bar(i), kr(j));
        [pout_sim(:, i)] = pout_sim_rice(gamma_th_sim, gamma_bar(i), samples, kr(j));
    end

    % Plot Outage Probability
    figure(fig_number)
    color = ['b', 'r', 'g'];
    for i = 1 : length(gamma_bar)
        semilogy(gamma_th_db, pout(:, i), 'Color', color(i), 'Linewidth', 1.5)
        hold on
        semilogy(gamma_th_sim_db, pout_sim(:, i), 'x', 'Color', 'k', 'Linewidth', 1.5)
        hold on
    end
    grid on
    ylabel('OP', 'Interpreter', 'Latex')
    xlabel('$\gamma_{th}$ (dB)', 'Interpreter', 'Latex')
    legend("$\bar{\gamma_{s}} = " + num2str(gamma_bar_db(1)) + "$dB",...
        "",...
        "$\bar{\gamma_{s}} = " + num2str(gamma_bar_db(2)) + "$dB",...
        "",...
        "$\bar{\gamma_{s}} = " + num2str(gamma_bar_db(3)) + "$dB",...
        "Simulated",...
        'Interpreter', 'Latex', 'Location', 'southeast')
    title("Rice ($K_r=" + num2str(kr(j)) + "$)", 'Interpreter', 'Latex')
    ylim([1e-6, 1e0])
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = 14;

    fig_number = fig_number + 1;
end


%% SEP Analisys
M = [4, 16, 64];

gamma_bar_db = linspace(-10, 30, points);
gamma_bar = 10 .^ (gamma_bar_db ./ 10);

gamma_bar_sim_db = linspace(-10, 30, N);
gamma_bar_sim = 10 .^ (gamma_bar_sim_db ./ 10);

sep = zeros(points, length(M));
sep_awgn = zeros(points, length(M));

for i = 1 : length(M)
    [sep(:, i), sep_awgn(:, i)] = sep_rayleigh(gamma_bar, M(i));
    [sep_sim(:, i), sep_sim_awgn(:, i)] = sep_sim_rayleigh(gamma_bar_sim, samples, M(i));
end

% Plot SEP
figure(fig_number)
color = ['b', 'r', 'g'];
for i = 1 : length(M)
    semilogy(gamma_bar_db, sep(:, i), 'Color', color(i), 'Linewidth', 1.5)
    hold on
    semilogy(gamma_bar_db, sep_awgn(:, i), '--', 'Color', color(i), 'Linewidth', 1.5)
    hold on
    semilogy(gamma_bar_sim_db, sep_sim(:, i), 'x', 'Color', 'k', 'Linewidth', 1.5)
    hold on
    semilogy(gamma_bar_sim_db, sep_sim_awgn(:, i), 'x', 'Color', 'k', 'Linewidth', 1.5)
    hold on
end
grid on
ylabel('SEP', 'Interpreter', 'Latex')
xlabel('$\bar{\gamma_{s}}$ (dB)', 'Interpreter', 'Latex')
ylim([1e-6, 1e0])
legend("$M = " + num2str(M(1)) + "$ Rayleigh",...
       "$M = " + num2str(M(1)) + "$ AWGN",...
       "",...
       "",...
       "$M = " + num2str(M(2)) + "$ Rayleigh",...
       "$M = " + num2str(M(2)) + "$ AWGN",...
       "",...
       "",...
       "$M = " + num2str(M(3)) + "$ Rayleigh",...
       "$M = " + num2str(M(3)) + "$ AWGN",...
       "Simulated",...
       "",...
       'Interpreter', 'Latex', 'Location', 'southwest')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 14;

%% Auxiliar Functions
function [pout] = pout_rayleigh(gamma_th, gamma_bar)
    pout = 1 - exp(-(gamma_th ./ gamma_bar));
end

function [pout] = pout_rice(gamma_th, gamma_bar, kr)
    pout = 1 - marcumq(sqrt(2*kr) , sqrt(2 * ((kr + 1) ./ gamma_bar) .* gamma_th));
end

function [pout] = pout_sim_rayleigh(gamma_th, gamma_bar, N)
    % pout = ones(1, length(gamma_th));
    for i = 1 : length(gamma_th)
        % normalized channel
        h = (1/sqrt(2)) * (randn(N, 1) + 1j*randn(N, 1));
        beta = abs(h);
        % gamma_s
        gamma_s = gamma_bar .* (beta .^ 2);
        % P(gamma_s < gamma_th)
        pout(i) = mean(gamma_s < gamma_th(i));
    end
end

function [pout] = pout_sim_rice(gamma_th, gamma_bar, N, kr)
    % pout = ones(1, length(gamma_th));
    for i = 1 : length(gamma_th)
        % line of sight component
        h_los = sqrt(kr/(kr+1));
        % non-line of sight component
        h_nlos = sqrt(1/(kr+1)) * (randn(N, 1) + 1j*randn(N, 1)) / sqrt(2);
        % normalized channel
        h = h_los + h_nlos;
        beta = abs(h);
        % gamma_s
        gamma_s = gamma_bar .* (beta .^ 2);
        % P(gamma_s < gamma_th)
        pout(i) = mean(gamma_s < gamma_th(i));
    end
end

function [sep, sep_awgn] = sep_rayleigh(gamma_bar, M)
    % sep
    cm = sqrt( (1.5 .* gamma_bar) ./ (M - 1 + (1.5 .* gamma_bar)) );
    cte = (sqrt(M) - 1) / sqrt(M);
    sep = (2 .* cte .* (1 - cm)) - ((cte .^ 2) .* (1 - ((4/pi) .* cm .* atan(1./cm))));
    % sep awgn
    q = qfunc(sqrt( (3 ./ (M - 1)) .* gamma_bar ));
    A = 4 * (1 - 1/sqrt(M));
    sep_awgn = A .* q - (A.^2 ./ 4) .* (q.^2);
end

function [sep, sep_awgn] = sep_sim_rayleigh(gamma_bar, N, M)
    % bits
    bits = randi([0, M-1], N, 1);
    % conversion to QAM symbols
    symbols = qammod(bits, M, 'UnitAveragePower', true);
    for i = 1 : length(gamma_bar)
        % normalized channel
        h = (1/sqrt(2)) * (randn(N, 1) + 1j*randn(N, 1));
        beta = abs(h);
        % normalized noise power, Es = 1
        sigma = 1 ./ gamma_bar(i);
        % awgn noise
        noise = sqrt(sigma/2) .* (randn(N, 1) + 1j*randn(N, 1));
        % received signal
        received = (beta .* symbols) + noise;
        % awgn
        received_awgn = symbols + noise;
        % equalization
        equalized = received ./ beta;
        % demodulation
        received_bits = qamdemod(equalized, M, 'UnitAveragePower', true);
        % awgn demodulation
        received_bits_awgn = qamdemod(received_awgn, M, 'UnitAveragePower', true);
        % SEP
        sep(i) = mean(bits ~= received_bits);
        % SEP awgn
        sep_awgn(i) = mean(bits ~= received_bits_awgn);
    end
end

