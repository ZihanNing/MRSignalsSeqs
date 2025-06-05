% exmaple of TBW with RF pulse
%%  Time domain
% Parameters
N = 100;
dt = 1;  % Sample interval (arbitrary units)
TBWs = [4, 8, 12, 16];

% Function handle for msinc
% f = @(TBW) msinc(N, TBW/2);  % ncyc = TBW / 2
f = @(TBW) ori_sinc(N, TBW/2);  % ncyc = TBW / 2

% Generate RF pulses
h1 = f(TBWs(1));
h2 = f(TBWs(2));
h3 = f(TBWs(3));
h4 = f(TBWs(4));

% Plot time domain
figure;
subplot(2,2,1), plot(h1), title('TBW = 4'), xlabel('Time'), ylabel('Amplitude')
subplot(2,2,2), plot(h2), title('TBW = 8'), xlabel('Time'), ylabel('Amplitude')
subplot(2,2,3), plot(h3), title('TBW = 12'), xlabel('Time'), ylabel('Amplitude')
subplot(2,2,4), plot(h4), title('TBW = 16'), xlabel('Time'), ylabel('Amplitude')


%% Frequency domain
% Frequency analysis
Fs = 1/dt;  % Sampling frequency
f_axis = linspace(-Fs/2, Fs/2, N);

% Compute FFTs
H1 = fftshift(fft(h1));
H2 = fftshift(fft(h2));
H3 = fftshift(fft(h3));
H4 = fftshift(fft(h4));

% Normalize magnitude
H1 = abs(H1) / max(abs(H1));
H2 = abs(H2) / max(abs(H2));
H3 = abs(H3) / max(abs(H3));
H4 = abs(H4) / max(abs(H4));

% Plot frequency domain
figure;
subplot(2,2,1), plot(f_axis, H1), ylim([0 1.5]),title('TBW = 4'), xlabel('Frequency'), ylabel('|H(f)|')
subplot(2,2,2), plot(f_axis, H2), ylim([0 1.5]),title('TBW = 8'), xlabel('Frequency'), ylabel('|H(f)|')
subplot(2,2,3), plot(f_axis, H3), ylim([0 1.5]),title('TBW = 12'), xlabel('Frequency'), ylabel('|H(f)|')
subplot(2,2,4), plot(f_axis, H4), ylim([0 1.5]),title('TBW = 16'), xlabel('Frequency'), ylabel('|H(f)|')

%% Support function
% Function to create a windowed sinc
function h = msinc(N, ncyc)
    x = (((0:N-1)+0.5) - (N/2 - 0.5)) / N * 2 * ncyc;
    hw = hamming(N);
    h = sinc(x(:)) .* hw(:);
    h = h';
end

function h = ori_sinc(N, ncyc)
    x = (((0:N-1)+0.5) - (N/2 - 0.5)) / N * 2 * ncyc;
    h = sinc(x(:));
    h = h';
end
