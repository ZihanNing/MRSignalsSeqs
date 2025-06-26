% Lecture 3 - Example 9
%
% Example with real RF pulse and steady state.
% Same RF, now excitation/recovery sequence.
%
% Starts the same as prior example...
clear;
close all;
clc;

%% 0. setting the RF pulse (in GRE)
plotlevel = 3;		% Level of detail to plot (<2, do not plot about RF pulse)
dt = .004;		% ms, sample spacing
tip = 2;		% desired tip angle
phi = 0; % phase of the RF, deg; tan(phi)=y/x; 0 means RF y; -90 means RF x
t = [-5:dt:5];		% extended time period to get 100 Hz spectral res.
tplot = find(abs(t)<=1); % only work on central 2 ms
rf=0*t;			% allocate RF
rf(tplot) = msinc(length(tplot),3);	% put RF pulse in central part. (windowed sinc -> usually used)
rf = rf/(sum(rf)*dt)*(tip/360)/42.58;	% scale for desired flip.

% frequency domain
rff = fftshift(fft(fftshift(rf)));		% RF profile, without scaling 
f = ([1:length(rff)]-length(rff)/2)/length(rff)/dt;	% kHz

if (plotlevel>2)
  figure;		% Plot the RF and profile.
  subplot(2,1,1);
  plot(t(tplot),rf(tplot)); xlabel('time(ms)'); ylabel('B1 (mT)'); title('Time domain')
  subplot(2,1,2);
  freqplot = find(abs(f)<5);			% Plot only range of freqs.
  plot(f(freqplot),abs(rff(freqplot)));	
  title('Frequency domain'); ylabel('B1 (mT)'); xlabel('Freq (kHz)');
end

%% 1. frequency profile (flip angle in each freq)
% To note: 
% Here I show two ways for frequency profile computation
% a. small tip approximation (only for FA <= 30 deg)
% b. hard pulse approximation (regardless to FA range)

Mf = ones(3,length(t),length(f)); % M in terms of position & time
Mf(1:2,:)=0;	
% initailize flip angle in freq (freq profile)
tip_f_small = zeros(length(f)); % small tip approximation
tip_f_hard = zeros(length(f)); % hard pulse approximation

% A. SMALL TIP APPROXIMATION
rfp = fftshift(fft(fftshift(rf)))*dt*42.58;		% RF profile, scaled (small tip approximation)
% flip angle approximation
tip_f_small = abs(rfp)*360;
% profile
for fi = 1:length(f)
  % Gradient rotation same for each interval
  Rgrad = zrot(f(fi)*dt*360);	% zrot(a) = rotation matrix.
  alpha = tip_f_small(fi);
  Rrf = throt(alpha,phi);
  Mf(:,end,fi) = Rrf*Rgrad*Mf(:,end,fi);		% Apply RF, Gradient
end
Mf_small = Mf; 



% B. HARD PULSE APPROXIMATION
% profile
% Note that we neglect relaxation during the RF.
for fi = 1:length(f)
  % Gradient rotation same for each interval
  Rgrad = zrot(f(fi)*dt*360);	% zrot(a) = rotation matrix.

  for ti = 2:length(t)
    % Hard Pulse Approximation...
    alpha = rf(ti)*dt*42.58*360;		% RF rotation over interval
    Rrf = throt(alpha,phi);

    Mf(:,ti,fi) = Rrf*Rgrad*Mf(:,ti-1,fi);		% Apply RF, Gradient

  end
end
% flip angle computed based on hard pulse approximation
Mf_final = squeeze(Mf(:,end,:));
Mz = Mf_final(3,:);
Mmag_sq = sum(Mf_final.^2,1);
tip_f_hard = rad2deg(acos(Mz ./ sqrt(Mmag_sq)));
Mf_hard = Mf;

% plot
figure;
subplot(211)
freqplot = find(abs(f)<5);			% Plot only range of freqs.
plot(f(freqplot),tip_f_small(freqplot),'-r');	
hold on
plot(f(freqplot),tip_f_hard(freqplot),'-b');	
legend('small tip appro','hard pulse appro')
title(['Flip angle: ',num2str(tip),' (deg)']); ylabel('Flip (deg)'); xlabel('Freq (kHz)');

subplot(212)
plot(f(freqplot),abs(squeeze(Mf_small(1,end,freqplot)+i*Mf_small(2,end,freqplot))),'-r')
hold on
plot(f(freqplot),abs(squeeze(Mf_hard(1,end,freqplot)+i*Mf_hard(2,end,freqplot))),'-b')
legend('small tip appro','hard pulse appro')
title(['Freq profile: ',num2str(tip),' (deg)']); ylabel('M_{xy}'); xlabel('Freq (kHz)');

%% 2. Bloch Simulation of GRE sequence
% |-RF 30 deg---ACQ---RF 30 deg-....---|

TR = 5;			% ms
T1 = 500;		% ms
T2 = 50;		% ms

% Bloch simulation

Gz = 2.3;	% mT/m  (gamma/2pi*Gz ~ 100 kHz/m, 1kHz/cm);
df = 0;		% kHz, off-resonance.

pos = [-.05:.0001:.05];	% 	Positions to simulate

M = ones(3,length(tplot),length(pos));
M(1:2,:)=0;				% M=[1;0;0];

% -- Start right after RF pulse - here's where steady state will be!
% -- First set A,B for relaxation over TR

A = diag([exp(-(TR-2)/T2) exp(-(TR-2)/T2) exp(-(TR-2)/T1)]);
Bs = [0;0;1-exp(-(TR-2)/T1)];
As = [0 0 0;0 0 0;0 0 1]*A;	% Perfect Spoiling


% -- Relaxation during RF pulse:  A/B for each interval.

At = diag([exp(-(dt)/T2) exp(-(dt)/T2) exp(-(dt)/T1)]);
Bt = [0;0;1-exp(-(dt)/T1)];


% -- Hard pulse approximation, with relaxation.
for z = 1:length(pos)
  %tt = sprintf('Position %d of %d',z,length(pos)); disp(tt); % Show progress
  % Gradient rotation same for each interval
  Rgrad = zrot((42.58*pos(z)*Gz+df)*dt*360);	% zrot(a) = rotation matrix.
  A = As;	% Set to start value again!
  B = Bs;	% Set to start value again.
   
  for ti = 2:length(tplot)
    % Hard Pulse Approximation...
    alpha = rf(tplot(ti))*dt*42.58*360;		% RF rotation over interval
    Rrf = yrot(alpha);

    M(:,ti,z) = Rrf*Rgrad*M(:,ti-1,z);		% Apply RF, Gradient

    Ai = At*Rgrad*Rrf; % Relaxation during RF
    Bi = Bt;

    % Propagate A,B
    A = Ai*A;
    B = Ai*B+Bi;
   
  end

  MM(:,z) = inv(eye(3)-A)*B;	% Steady-state magnetization.
  MMe(:,z) = A*[0;0;1]+B;	% Magnetization after RF.  Note that
				% if we start with M=M0 and propagate
				% through 1 TR, there's no recovery before
				% the RF, so we just get the excitation profile.

end

Mxy = MM(1,:)+i*MM(2,:);
Mxye = MMe(1,:)+i*MMe(2,:);

figure;
plot(pos*100,abs(Mxy),'k-',pos*100,real(Mxy),'b-',pos*100,imag(Mxy),'r:');
grid on;
xlabel('Position (cm)');
ylabel('M_{xy}');
title('Excitation/Recovery Signal vs Position');
legend('abs(M_{xy})','real(M_{xy})','imag(M_{xy})')
setprops;

if (plotlevel>2)
  figure;
  plot(pos*100,180/pi*asin(abs(Mxye)));		% Plot of *ANGLE*
  grid on;
  xlabel('Position (cm)');
  ylabel('Flip Angle (deg)');
  title('Excitation Profile vs Position');
  setprops;
end

if (plotlevel>2)
  figure;
  fplot(@(x)sin(pi./180.*x).*(1-exp(-5./500))./(1-exp(-5./500).*cos(pi./180.*x)),[0,60]);
  %fplot('sin(pi/180*x)*(1-exp(-5/500))/(1-exp(-5/500)*cos(pi/180*x))',[0,60]);
  lplot('Signal','Flip Angle (deg)','Signal vs Flip Angle');
  setprops;
end

if (plotlevel < 4) disp('Try changing plotlevel for more/less detail.'); end;

