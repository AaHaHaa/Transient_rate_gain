function At_noise = shot_noise(Nt,dt,sim,cs)
%SHOT_NOISE It computes the shot noise included in the governing equation

h = 6.62607015e-34; % J*s

time_window = Nt*dt; % ps
f = ifftshift((-Nt/2:Nt/2-1)'/time_window,1); % THz
real_f = (f*cs+sim.f0)*1e12; % Hz

% I use analytical-signal representation for solving GMMNLSE, so the field
% covers only the positive frequencies.
real_f(real_f<0) = 0; % no noise at negative frequencies

noise_amplitude = sqrt(cs*h.*real_f/(time_window*1e-12));

%At_noise = fft(noise_amplitude.*randn(Nt,1).*exp(1i*2*pi*rand(Nt,1)),[],1);
At_noise = fft(noise_amplitude.*exp(1i*2*pi*rand(Nt,1)),[],1);

end