%%% Given the RF Pulse
%%% Simulate:
%%% - Slice profile (time & position)
%%% - Frequency profile (time & frequency)
%%%
%%% Zihan @ King's
%%% 17-Mar-2025

%% Define the RF pulse

dt = .004;		% ms, sample spacing
tip = 90;		% desired tip.	
t = [-5:dt:5];		% extended time period to get 100 Hz spectral res.
tplot = find(abs(t)<=1); % only work on central 2 ms
rf=0*t;			% allocate RF
RF_type = 'non-selective'; % ZN: could be 'selective' and 'non-selective'
switch RF_type
    case 'selective'
        disp('Selective RF pulse - shape sinc');
        rf(tplot) = msinc(length(tplot),3);	% put RF pulse in central part. % ZN: sinc RF pulse
    case 'non-selective'
        disp('Non-selective RF pulse - shape rect');
        rf(tplot) = 1;
    otherwise
        error('Unrecognize RF pulse!');
end

rf = rf/(sum(rf)*dt)*(tip/360)/42.58;	% scale for desired flip.

rfp = fftshift(fft(fftshift(rf)))*dt*42.58;		% RF profile, scaled
f = ([1:length(rfp)]-length(rfp)/2)/length(rfp)/dt;	% kHz

figure(1);		% Plot the RF and profile.
subplot(3,1,1);
switch RF_type
    case 'selective'
        plot(t(tplot),rf(tplot)); xlabel('time(ms)'); ylabel('B1 (mT)'); title('B1(t)')
    case 'non-selective'
        plot(t,rf); xlabel('time(ms)'); ylabel('B1 (mT)'); title('B1(t)')
end
subplot(3,1,2);
freqplot = find(abs(f)<5);			% Plot only range of freqs.
plot(f(freqplot),abs(rfp(freqplot))*360);	
title('Small Tip Approximation'); ylabel('Flip (deg)'); xlabel('Freq (kHz)');


%% Bloch simulation (slice profile)
switch RF_type
    case 'selective'
        Gz = 2.5;	% mT/m  (gamma/2pi*Gz ~ 100 kHz/m, 1kHz/cm; % slice-selective gradient
    case 'non-selective'
        Gz = 0; % no gradient applied if not selecting slice
end
off_f = 500e-3;		% kHz, off-resonance.
phi = 0; % deg, RF pulse, tan(phi)=y/x; 0 means RF y; -90 means RF x

% Simulate by time and position (slice profile)
pos = [-.05:.0001:.05];	% 	Positions to simulate
M = ones(3,length(tplot),length(pos)); % M in terms of position & time
M(1:2,:)=0;				% M=[0;0;1];

%plot Gz
subplot(3,1,3);
plot(pos*10,pos*Gz);
xlabel('Position(cm)'); ylabel('Delta_B0 (mT)'); title('Changed B0 along z axis')

% %plot on-resonance freq of every position
% subplot(4,1,4);
% gamma = 42.58; % kHz/mT
% on_freq = 2*pi*gamma*Gz/10; % kHz
% plot(pos*10,on_freq*pos);
% xlabel('Position(cm)'); ylabel('Delta_B0 (mT)'); title('Changed B0 along z axis')

% Simulate!!
% Note that we neglect relaxation during the RF.
for z = 1:length(pos)
  % Gradient rotation same for each interval
  Rgrad = zrot((42.58*pos(z)*Gz+off_f)*dt*360);	% zrot(a) = rotation matrix.
   
  for ti = 2:length(tplot)
    % Hard Pulse Approximation...
    alpha = rf(tplot(ti))*dt*42.58*360;		% RF rotation over interval
    Rrf = throt(alpha,phi);

    M(:,ti,z) = Rrf*Rgrad*M(:,ti-1,z);		% Apply RF, Gradient

  end
end

% plot
figure;
subplot(3,2,1);
plot(t(tplot),squeeze(M(1,:,ceil(length(pos)/2)))); 
ylim([-1 1])
xlabel('Time (ms)'); ylabel('M_x(t)'); title('M_x(t) in the centre of pos')

subplot(3,2,3);
plot(t(tplot),squeeze(M(2,:,ceil(length(pos)/2))));
ylim([-1 1])
xlabel('Time (ms)'); ylabel('M_y(t)'); title('M_y(t) in the centre of pos')

subplot(3,2,5);
plot(t(tplot),squeeze(M(3,:,ceil(length(pos)/2))));
ylim([-1 1])
xlabel('Time (ms)'); ylabel('M_z(t)'); title('M_z(t) in the centre of pos')

subplot(3,2,2)
plot(pos*100,squeeze(M(1,end,:)),'b--'); hold on;
plot(pos*100,squeeze(M(2,end,:)),'r-'); hold off;
xlabel('Position (cm)'); ylabel('M_x(z) and M_y(z)'); 
title('M_x(z) and M_y(z) at the end of time (final state)')
legend('M_x(z)','M_y(z)')

subplot(3,2,4)
plot(pos*100,abs(squeeze(M(1,end,:)+i*M(2,end,:))));
xlabel('Position (cm)'); ylabel('M_{xy}(z)');
title('M_{xy}(z) at the end of time (final state)')

subplot(3,2,6)
plot(pos*100,squeeze(M(3,end,:)));
xlabel('Position (cm)'); ylabel('M_y(z)');
title('M_z(z) at the end of time (final state)')

% At the end of the RF pulse (at the center pos)
Mxy = M(1,end,ceil(length(pos)/2)) + i*M(2,end,ceil(length(pos)/2));
Mz = M(3,end,ceil(length(pos)/2));
fprintf('=========== At the end of the RF pulse (at centre pos): \n');
fprintf('Expected RF pulse: flip angle %d deg- phi %d deg \n',tip,phi);
fprintf('Off-resonance freq: %d Hz\n',off_f*1000);
fprintf('--Simu Result:\n');
fprintf('Mag of the signal: abs(Mxy) = %.2f \n',abs(Mxy));
fprintf('Phs of the signal: angle(Mxy) = %.2f * pi \n',angle(Mxy)/pi);
fprintf('Residual Mz: abs(Mz) = %.2f \n',abs(Mz));
fprintf('Actual flip angle: %.1f deg \n',rad2deg(atan(abs(Mxy)/abs(Mz))));
fprintf('Phase shift: %.1f Hz \n',(phi+0.5*pi - angle(Mxy))/2*pi);

%% Bloch simulation (freq profile)
% Simulate by BW and frequency (freq profile)
BW = 2; % kHz
df = linspace(-BW,BW); % frequency to simulate
Mf = ones(3,length(tplot),length(df)); % M in terms of position & time
Mf(1:2,:)=0;		

% Simulate!!
% Note that we neglect relaxation during the RF.
for f = 1:length(df)
  % Gradient rotation same for each interval
  Rgrad = zrot(df(f)*dt*360);	% zrot(a) = rotation matrix.
   
  for ti = 2:length(tplot)
    % Hard Pulse Approximation...
    alpha = rf(tplot(ti))*dt*42.58*360;		% RF rotation over interval
    Rrf = throt(alpha,phi);

    Mf(:,ti,f) = Rrf*Rgrad*Mf(:,ti-1,f);		% Apply RF, Gradient

  end
end

% Plot
figure;
subplot(2,1,1)
plot(df,squeeze(Mf(1,end,:)),'b--'); hold on;
plot(df,squeeze(Mf(2,end,:)),'r-'); hold off;
xlabel('freq (kHz)'); ylabel('M_x(f) and M_y(f)'); 
title('M_x(f) and M_y(f) at the end of time (final state)')
legend('M_x(f)','M_y(f)')

subplot(2,1,2)
plot(df,abs(squeeze(Mf(1,end,:)+i*Mf(2,end,:))),'b--');hold on;
plot(df,squeeze(Mf(3,end,:)),'r-');hold off;
xlabel('freq (kHz)'); ylabel('M_{xy}(f)');
title('M_{xy}(f) at the end of time (final state)')
legend('M_xy(f)','M_z(f)')



