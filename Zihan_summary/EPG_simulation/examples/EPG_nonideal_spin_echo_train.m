% EPG simulation for non-ideal spin echo
% train
%
% modified by Zihan Ning @ King's College London
% 21-Aug-2025


%% Simulate
disp('-- Spin Echo with non-180x pulses --')
T1 = 195; T2=95; T=10;     % Chosen for E1=0.95, E2=0.9

% define refocusing angle(s)
refocusing_angle = [60,60]; % in deg
refocusing_angle = deg2rad(refocusing_angle); % in rad


M_whole = zeros(3,9);
Mxy_whole = zeros(1,9);

Q = epg_m0(2); disp('Equilibrium: Q ='); disp(Q); 
disp('M ='); M = epg_FZ2spins(Q); disp(M); disp('combined M =');[M_whole(:,1),Mxy_whole(1)] = epg_combM_convfrmQ(M); disp(M_whole(:,1));

Q = epg_rf(Q,pi/2,pi/2); disp('Excitation-90: Q = '); disp(Q); 
disp('M ='); M = epg_FZ2spins(Q); disp(M); disp('combined M =');[M_whole(:,2),Mxy_whole(2)] = epg_combM_convfrmQ(M); disp(M_whole(:,2));

Q = epg_relax(epg_grad(Q,1),T1,T2,T); disp(['Pre-',num2str(rad2deg(refocusing_angle(1))),': Q =']); disp(Q); 
disp('M ='); M = epg_FZ2spins(Q); disp(M); disp('combined M ='); [M_whole(:,3),Mxy_whole(3)] = epg_combM_convfrmQ(M); disp(M_whole(:,3));

Q = epg_rf(Q,refocusing_angle(1),0); disp(['Post-',num2str(rad2deg(refocusing_angle(1))),': Q =']); disp(Q); 
disp('M ='); M = epg_FZ2spins(Q); disp(M); disp('combined M ='); [M_whole(:,4),Mxy_whole(4)] = epg_combM_convfrmQ(M); disp(M_whole(:,4));

Q = epg_relax(epg_grad(Q,1),T1,T2,T); disp('1st Spin Echo: Q ='); disp(Q); 
disp('M ='); M = epg_FZ2spins(Q); disp(M); disp('combined M ='); [M_whole(:,5),Mxy_whole(5)] = epg_combM_convfrmQ(M); disp(M_whole(:,5));

Q = epg_relax(epg_grad(Q,1),T1,T2,T); disp(['Pre-',num2str(rad2deg(refocusing_angle(2))),': Q =']); disp(Q); 
disp('M ='); M = epg_FZ2spins(Q); disp(M); disp('combined M ='); [M_whole(:,6),Mxy_whole(6)] = epg_combM_convfrmQ(M); disp(M_whole(:,6));

Q = epg_rf(Q,refocusing_angle(2),0); disp(['Post-',num2str(rad2deg(refocusing_angle(1))),': Q =']); disp(Q); 
disp('M ='); M = epg_FZ2spins(Q); disp(M); disp('combined M ='); [M_whole(:,7),Mxy_whole(7)] = epg_combM_convfrmQ(M); disp(M_whole(:,7));

Q = epg_relax(epg_grad(Q,1),T1,T2,T); disp('2nd Spin Echo: Q ='); disp(Q); 
disp('M ='); M = epg_FZ2spins(Q); disp(M); disp('combined M ='); [M_whole(:,8),Mxy_whole(8)] = epg_combM_convfrmQ(M); disp(M_whole(:,8));

%% show in M
figure;
subplot(3,1,1); 
    plot(1:size(M_whole,2),M_whole(3,:),'o-'),...
        hold; plot(1:size(M_whole,2),zeros(size(M_whole,2)),'--r'),...
        xlabel('states, 5=1st echo, 8=2nd echo'),ylabel('relative M'),title(['refocusing angle = ',num2str(rad2deg(refocusing_angle)),'; M_z'])
subplot(3,1,2); 
    plot(1:length(Mxy_whole),abs(Mxy_whole),'o-'),...
    hold; plot(1:size(M_whole,2),zeros(size(M_whole,2)),'--r'),...    
    xlabel('states, 5=1st echo, 8=2nd echo'),ylabel('relative M'),title('|M_{xy}|')
subplot(3,1,3);
    plot(1:length(Mxy_whole),angle(Mxy_whole),'o-'),...
    hold; plot(1:size(M_whole,2),zeros(size(M_whole,2)),'--r'),...
    xlabel('states, 5=1st echo, 8=2nd echo'),ylabel('relative M'),title('phs(M_{xy})')


