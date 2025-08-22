% EPG simulation for demonstrating stimulated echo train
%
% modified by Zihan Ning @ King's College London
% 22-Aug-2025


%% Simulate
disp('-- Stimulation of stimulated echo --')
% T1 = 195; T2=95; T=10;     % Chosen for E1=0.95, E2=0.9
% we ignore relaxation in this demonstration

% define refocusing angle(s)
rf_train = [60,60,60]; % in deg
rf_train = deg2rad(rf_train); % in rad

M_whole = zeros(3,7);
Mxy_whole = zeros(1,7);

Q = epg_m0(5); disp('Equilibrium: Q ='); disp(Q); 
disp('M ='); M = epg_FZ2spins(Q); disp(M); disp('combined M =');[M_whole(:,1),Mxy_whole(1)] = epg_combM_convfrmQ(M); disp(M_whole(:,1));

Q = epg_rf(Q,rf_train(1),0); disp('Excitation-60: Q = '); disp(Q); 
disp('M ='); M = epg_FZ2spins(Q); disp(M); disp('combined M =');[M_whole(:,2),Mxy_whole(2)] = epg_combM_convfrmQ(M); disp(M_whole(:,2));

Q = epg_grad(Q,1); disp('First gradient: Q = '); disp(Q); 
disp('M ='); M = epg_FZ2spins(Q); disp(M); disp('combined M =');[M_whole(:,3),Mxy_whole(3)] = epg_combM_convfrmQ(M); disp(M_whole(:,3));

Q = epg_rf(Q,rf_train(2),0); disp('Excitation-60: Q = '); disp(Q); 
disp('M ='); M = epg_FZ2spins(Q); disp(M); disp('combined M =');[M_whole(:,4),Mxy_whole(4)] = epg_combM_convfrmQ(M); disp(M_whole(:,4));

Q = epg_grad(Q,1); disp('2nd gradient [echo]: Q = '); disp(Q); 
disp('M ='); M = epg_FZ2spins(Q); disp(M); disp('combined M =');[M_whole(:,5),Mxy_whole(5)] = epg_combM_convfrmQ(M); disp(M_whole(:,5));

Q = epg_rf(Q,rf_train(3),0); disp('Excitation-60: Q = '); disp(Q); 
disp('M ='); M = epg_FZ2spins(Q); disp(M); disp('combined M =');[M_whole(:,6),Mxy_whole(6)] = epg_combM_convfrmQ(M); disp(M_whole(:,6));

Q = epg_grad(Q,1); disp('3rd gradient [echo]: Q = '); disp(Q); 
disp('M ='); M = epg_FZ2spins(Q); disp(M); disp('combined M =');[M_whole(:,7),Mxy_whole(7)] = epg_combM_convfrmQ(M); disp(M_whole(:,7));

%% show in M
figure;
subplot(3,1,1); 
    plot(1:size(M_whole,2),M_whole(3,:),'o-'),...
        hold; plot(1:size(M_whole,2),zeros(size(M_whole,2)),'--r'),...
        xlabel('states, 2,4,6=RF, 5,7=echo'),ylabel('relative M'),title(['FA = ',num2str(rad2deg(rf_train(1))),'; M_z'])
subplot(3,1,2); 
    plot(1:length(Mxy_whole),abs(Mxy_whole),'o-'),...
    hold; plot(1:size(M_whole,2),zeros(size(M_whole,2)),'--r'),...    
    xlabel('states, 2,4,6=RF, 5,7=echo'),ylabel('relative M'),title('|M_{xy}|')
subplot(3,1,3);
    plot(1:length(Mxy_whole),angle(Mxy_whole),'o-'),...
    hold; plot(1:size(M_whole,2),zeros(size(M_whole,2)),'--r'),...
    xlabel('states, 2,4,6=RF, 5,7=echo'),ylabel('relative M'),title('phs(M_{xy})')




