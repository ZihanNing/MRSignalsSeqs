% This will be an example for non-adiabatic RF pulse design

%% 1. Small FA with given RF duration & FA
% windowed vs non-windowed sinc function (windowed preferred as better
% slice profile)
gammabar = 42.58; % kHz/mT

% RF pulse design
% given parameters for design
flip = 60; % FA, deg
phi = 0; % phase of the RF, deg; tan(phi)=y/x; 0 means RF y; -90 means RF x
TBW = 8; % four pulse shape (sinc) cross 0 point for the sinc or msinc during the RF duration
Trf = 4; % RF duration, ms 

% Set the frame
% freq
BWplot = 1; % kHz -> BWrf = TBW/Trf = 1 kHz
df = linspace(-BWplot,BWplot);
% time
N = length(df);
IN = [-N/2:N/2-1]/N; % put RF pulse in the middle of the time axis
t = IN*Trf; dt = Trf/N;

% RF definition
% shape
RF_sinc_shape = sinc(IN * TBW); % sinc without windowing
RF_msinc_shape = hamming(N)' .* sinc(IN * TBW); % sinc with windowing, an alternative
% modulate by FA
RF_sinc = (flip*pi/180) * RF_sinc_shape/sum(RF_sinc_shape) / (2*pi*gammabar*dt); % FA = gamma*sum(dt*RF)dt; B1, mT
RF_msinc = (flip*pi/180) * RF_msinc_shape/sum(RF_msinc_shape) / (2*pi*gammabar*dt); % FA = gamma*sum(dt*RF)dt; B1, mT

% To plot
figure;
% RF pulse 
subplot(221),plot(IN*dt,RF_sinc),xlabel('Time [ms]'),ylabel('B1'),title('Sinc RF pulse [time]')
RF_sinc_freq = fftshift(fft(RF_sinc));
subplot(222),plot(df,abs(RF_sinc_freq)/max(abs(RF_sinc_freq))),xlabel('Freq [kHz]'),ylabel('B1'),ylim([0 1.2]),title('Sinc RF pulse [freq]')
subplot(223),plot(IN*dt,RF_msinc),xlabel('Time [ms]'),ylabel('B1'),title('mSinc RF pulse [time]')
RF_msinc_freq = fftshift(fft(RF_msinc));
subplot(224),plot(df,abs(RF_msinc_freq)/max(abs(RF_msinc_freq))),xlabel('Freq [kHz]'),ylabel('B1'),ylim([0 1.2]),title('mSinc RF pulse [freq]')

%% Frequency profile
Mf = ones(3,length(t),length(df)); % M in terms of position & time
Mf(1:2,:)=0;	

figure;

for typ = 1:2
    % select the RF for simulation
    if typ==1; rf = RF_sinc;else; rf = RF_msinc; end

    % Simulate!!
    % Note that we neglect relaxation during the RF.
    for f = 1:length(df)
      % Gradient rotation same for each interval
      Rgrad = zrot(df(f)*dt*360);	% zrot(a) = rotation matrix.

      for ti = 2:length(t)
        % Hard Pulse Approximation...
        alpha = rf(ti)*dt*gammabar*360;		% RF rotation over interval
        Rrf = throt(alpha,phi);

        Mf(:,ti,f) = Rrf*Rgrad*Mf(:,ti-1,f);		% Apply RF, Gradient

      end
    end
    
    % Plot
    subplot(2,1,typ)
    plot(df,abs(squeeze(Mf(1,end,:)+i*Mf(2,end,:))),'b--');hold on;
    plot(df,squeeze(Mf(3,end,:)),'r-');hold off;
    xlabel('freq (kHz)'); ylabel('M_{xy}(f)');
    if typ==1
        title('rf sinc:: M_{xy}(f) at the end of time (final state)')
    else
        title('rf msinc:: M_{xy}(f) at the end of time (final state)')
    end
    legend('M_xy(f)','M_z(f)')
end


%% Slice Profile
Gz = 2.5;	% mT/m  (gamma/2pi*Gz ~ 100 kHz/m, 1kHz/cm;
off_freq = 0;		% kHz, off-resonance.
pos = [-.05:.0001:.05];	% 	Positions to simulate

% RF rotation (excitation bloch simulation)
M = ones(3,length(IN),length(pos));
M(1:2,:)=0;				% M=[0;0;1];

for typ = 1:2
    % select the RF for simulation
    if typ==1; rf = RF_sinc;else; rf = RF_msinc; end

    for z = 1:length(pos)
      % Gradient rotation same for each interval
      Rgrad = zrot((gammabar*pos(z)*Gz+off_freq)*dt*360);	% zrot(a) = rotation matrix.

      for ti = 2:length(IN)
        % Hard Pulse Approximation...
        alpha = rf(ti)*dt*gammabar*360;		% RF rotation over interval
        Rrf = throt(alpha,phi);

        M(:,ti,z) = Rrf*Rgrad*M(:,ti-1,z);		% Apply RF, Gradient

      end
    end
    
    % To plot
    % plot
    figure;
    subplot(3,2,1);
    plot(t,squeeze(M(1,:,ceil(length(pos)/2)))); 
    ylim([-1 1])
    xlabel('Time (ms)'); ylabel('M_x(t)'); title('M_x(t) in the centre of pos')

    subplot(3,2,3);
    plot(t,squeeze(M(2,:,ceil(length(pos)/2))));
    ylim([-1 1])
    xlabel('Time (ms)'); ylabel('M_y(t)'); title('M_y(t) in the centre of pos')

    subplot(3,2,5);
    plot(t,squeeze(M(3,:,ceil(length(pos)/2))));
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
    fprintf('Expected RF pulse: flip angle %d deg- phi %d deg \n',flip,phi);
    fprintf('Off-resonance freq: %d Hz\n',off_freq*1000);
    fprintf('--Simu Result:\n');
    fprintf('Mag of the signal: abs(Mxy) = %.2f \n',abs(Mxy));
    fprintf('Phs of the signal: angle(Mxy) = %.2f * pi \n',angle(Mxy)/pi);
    fprintf('Residual Mz: abs(Mz) = %.2f \n',abs(Mz));
    fprintf('Actual flip angle: %.1f deg \n',rad2deg(atan(abs(Mxy)/abs(Mz))));
    fprintf('Phase shift: %.1f Hz \n',(phi+0.5*pi - angle(Mxy))/2*pi);
end





