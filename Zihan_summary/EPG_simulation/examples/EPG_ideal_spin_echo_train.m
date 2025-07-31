% Lecture 4, Example 02
% 
% Simple image reconstruction in matlab


%% Simulate
disp('-- Spin Echo with 180x pulses --')
T1 = 195; T2=95; T=10;     % Chosen for E1=0.95, E2=0.9

Q = epg_m0(2); disp('Equilibrium: Q ='); disp(Q); 
disp('M ='); M = epg_FZ2spins(Q); disp(M); disp('combined M =');[M0,Mxy0] = epg_combM_convfrmQ(M); disp(M0);

Q = epg_rf(Q,pi/2,pi/2); disp('Excitation-90: Q = '); disp(Q); 
disp('M ='); M = epg_FZ2spins(Q); disp(M); disp('combined M =');[M1,Mxy1] = epg_combM_convfrmQ(M); disp(M1);

Q = epg_relax(epg_grad(Q,1),T1,T2,T); disp('Pre-180: Q ='); disp(Q); 
disp('M ='); M = epg_FZ2spins(Q); disp(M); disp('combined M ='); [M3,Mxy3] = epg_combM_convfrmQ(M); disp(M3);

Q = epg_rf(Q,pi,0); disp('Post-180: Q ='); disp(Q); 
disp('M ='); M = epg_FZ2spins(Q); disp(M); disp('combined M ='); [M4,Mxy4] = epg_combM_convfrmQ(M); disp(M4);

Q = epg_relax(epg_grad(Q,1),T1,T2,T); disp('1st Spin Echo: Q ='); disp(Q); 
disp('M ='); M = epg_FZ2spins(Q); disp(M); disp('combined M ='); [M5,Mxy5] = epg_combM_convfrmQ(M); disp(M5);

Q = epg_relax(epg_grad(Q,1),T1,T2,T); disp('Pre-180: Q ='); disp(Q); 
disp('M ='); M = epg_FZ2spins(Q); disp(M); disp('combined M ='); [M6,Mxy6] = epg_combM_convfrmQ(M); disp(M6);

Q = epg_rf(Q,pi,0); disp('Post-180: Q ='); disp(Q); 
disp('M ='); M = epg_FZ2spins(Q); disp(M); disp('combined M ='); [M6,Mxy6] = epg_combM_convfrmQ(M); disp(M6);

Q = epg_relax(epg_grad(Q,1),T1,T2,T); disp('2nd Spin Echo: Q ='); disp(Q); 
disp('M ='); M = epg_FZ2spins(Q); disp(M); disp('combined M ='); [M7,Mxy7] = epg_combM_convfrmQ(M); disp(M7);

%% show in M


