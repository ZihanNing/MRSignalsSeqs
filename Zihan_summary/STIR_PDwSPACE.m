clc;
clear;

%% Tissue definitions
% Add more tissues by expanding these arrays
tissue_names = {'fat', 'blood', 'myocardio'};
T1s = [187, 1121, 700];    % ms
T2s = [92, 263, 58];       % ms

% === Sequence settings ===
TI = 110;          % Inversion time [ms]
ETL = 30;          % Echo train length
FA_excite = 90;    % Excitation flip angle [deg]
FA_refoc = 175;    % Refocusing flip angle [deg]
ESP = 4.36;        % Echo spacing [ms]
TR = 1400;         % Repetition time [ms]
TD = TR - TI - (ESP/2 + ESP * ETL);  % Delay to complete TR
niT = 30;          % Number of TRs to simulate (to reach steady-state)

% RF train (in radians)
rf_train = ones(1, ETL+1);
rf_train(1) = FA_excite;
rf_train(2:end) = FA_refoc;
rf_train = deg2rad(rf_train);

% Preallocate output
Mz_all = {};
Mxy_all = {};

%% Simulation loop for each tissue
for t = 1:length(tissue_names)
    tissue = tissue_names{t};
    T1 = T1s(t);
    T2 = T2s(t);

    disp(['Simulating: ', tissue]);

    % Initial magnetization
    Q = epg_m0(5);
    M = {}; Mxy = {};
    [M{end+1}, Mxy{end+1}] = epg2mag(Q);

    for n = 1:niT
        lastM = M{end}; lastMxy = Mxy{end};
        M = {lastM}; Mxy = {lastMxy};

        Q = epg_rf(Q, pi, 0);           % non-selective inversion of STIR
        Q = epg_grad(Q, 1);             [M{end+1}, Mxy{end+1}] = epg2mag(Q); % crusher behind STIR
        Q = epg_relax(Q, T1, T2, TI);   [M{end+1}, Mxy{end+1}] = epg2mag(Q); % TI

        % Initial excitation
        Q = epg_rf(Q, rf_train(1), pi/2);    % 90° pulse
        Q = epg_relax(epg_grad(Q,1),T1,T2,ESP/2);

        % TSE echo train
        for i = 2:length(rf_train)
            Q = epg_rf(Q, rf_train(i), 0);         % 180° refocusing
            Q = epg_relax(epg_grad(Q,1),T1,T2,ESP/2);
            [M{end+1}, Mxy{end+1}] = epg2mag(Q);   % read echo signal
            Q = epg_relax(epg_grad(Q,1),T1,T2,ESP/2);
        end

        Q = epg_grad(Q, 1);             [M{end+1}, Mxy{end+1}] = epg2mag(Q); % rewinder after readout
        Q = epg_relax(Q, T1, T2, TD);   [M{end+1}, Mxy{end+1}] = epg2mag(Q); % dead time for recovery after readout
    end

    % Store results for plotting
    Mz_all{t} = cellfun(@(m) m(3), M);
    Mxy_all{t} = cell2mat(Mxy);
end

%% Plotting results
figure;
subplot(3,1,1); hold on;
for t = 1:length(tissue_names)
    plot(Mz_all{t}, '-o','DisplayName', tissue_names{t});
end
yline(0, '--r'); legend; ylabel('M_z'); title('M_z over STIR-TSE block');

subplot(3,1,2); hold on;
for t = 1:length(tissue_names)
    plot(abs(Mxy_all{t}), '-o','DisplayName', tissue_names{t});
end
yline(0, '--r'); legend; ylabel('|M_{xy}|'); title('|M_{xy}| over STIR-TSE block');

subplot(3,1,3); hold on;
for t = 1:length(tissue_names)
    plot(angle(Mxy_all{t}), '-o','DisplayName', tissue_names{t});
end
yline(0, '--r'); legend; ylabel('phase(M_{xy})'); title('Phase of M_{xy} over STIR-TSE block');
