% This will be an example for non-adiabatic RF pulse design

%% 1. Small FA with given RF duration & FA
% windowed vs non-windowed sinc function (windowed preferred as better
% slice profile)
gammabar = 42.58; % kHz/mT

% RF pulse design
% given parameters for design
flip = 60; % FA, deg
phi = 0; % phase of the RF, deg; tan(phi)=y/x; 0 means RF y; -90 means RF x
TBW = 2; % four pulse shape (sinc) cross 0 point for the sinc or msinc during the RF duration
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

% Gradient definition
Gz = 2.5;	% mT/m  (gamma/2pi*Gz ~ 100 kHz/m, 1kHz/cm;
off_freq = 0;		% kHz, off-resonance.
pos = [-.05:.0001:.05];	% 	Positions to simulate

% RF rotation (excitation bloch simulation)
M = ones(3,length(IN),length(pos));
M(1:2,:)=0;				% M=[0;0;1];
% select the RF for simulation
rf = RF_sinc;

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
figure;
% RF pulse 
if isequal(rf,RF_sinc) % sinc without windowing
    plot(IN*dt,RF_sinc),xlabel('Time [ms]'),ylabel('B1'),title('Sinc RF pulse [time]')
    RF_sinc_freq = fftshift(fft(RF_sinc));
    plot(df,abs(RF_sinc_freq)),xlabel('Freq [kHz]'),ylabel('B1'),title('Sinc RF pulse [freq]')
else
    plot(IN*dt,RF_msinc),xlabel('Time [ms]'),ylabel('B1'),title('mSinc RF pulse [time]')
    RF_msinc_freq = fftshift(fft(RF_msinc));
    plot(df,abs(RF_msinc_freq)/max(abs(RF_msinc_freq))),xlabel('Freq [kHz]'),ylabel('B1'),title('mSinc RF pulse [freq]')
end







gammarbar = 42.58; % kHz/mT

M0 = [0,0,1];
dt = 0.1; % ms
