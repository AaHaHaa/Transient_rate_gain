function func = narrowband_scaledFT
%NARROWBAND_SCALEDFT A container for narrowband transformation following
%the scaled Fourier transform
%
%   convert(At,cs) computes the narrowband-transformed field
%   recover(transformed_At,cs) computes the inverse narrowband-transformed field; that is, it recovers the field

func.convert = @convert;
func.recover = @recover;

end

function transformed_At = convert(At,cs)
%CONVERT

% If At is all-zero, just downsample it by rate of cs
if all(abs(At) == 0)
    transformed_At = At(1:cs:end);
    return
end

Nt = size(At,1);

if cs < 1
    error('narrowband_scaledFTError:incorrect_cs',...
          'cs must be larger than or equal to 1.');
end
if mod(Nt,cs) ~= 0
    error('narrowband_scaledFTError:incorrect_cs',...
          'cs must be a divisor of Nt, the number of sampling points of At.');
end
if numel(At) ~= Nt
    error('narrowband_scaledFTError:AtError',...
          'Current implementation of narrowband transformation supports only a single-mode scalar field.');
end

[~,peak_position] = max(sum(abs(At).^2,2),[],1);
avg_center = sum((1:Nt)'.*abs(At).^2,[1,2])/sum(abs(At(:)).^2);
% If two positions are too far away, then it might indicate that the pulse
% is center at the left edge of the time window such that its other half is
% at the right edge of the window, due to periodic boundary condition under
% numerical DFT.
% If two positions are close, then we pick the avg_center as the pulse
% center.
if abs(peak_position - avg_center) > Nt/4
    fftshift_At = fftshift(At,1);
    avg_center = sum((1:Nt)'.*abs(fftshift_At).^2,[1,2])/sum(abs(fftshift_At(:)).^2);
    if avg_center >= floor(Nt/2)+1
        avg_center = avg_center - floor(Nt/2);
    else
        avg_center = avg_center + floor(Nt/2);
    end
end
pulse_center = avg_center;

phase_shift = ifftshift(2*pi*(1:Nt)'/Nt*pulse_center,1); % phase shift due to temporal offset

Aw = ifft(At,[],1);

phase = ifftshift(unwrap(angle(fftshift(Aw.*exp(-1i*phase_shift),1))),1);
phase = phase - phase(1,:,:);

transformed_phase = phase(1:cs:end)/cs;
transformed_Aw = sqrt(cs)*abs(Aw(1:cs:end)).*exp(1i*transformed_phase);

transformed_phase_shift = ifftshift(2*pi*(1:Nt/cs)'/(Nt/cs)*(pulse_center/cs),1); % phase shift due to temporal offset
transformed_At = fft(transformed_Aw.*exp(1i*transformed_phase_shift),[],1); % add the temporal offset back

end

function At = recover(transformed_At,cs)
%RECOVER

Nt = size(transformed_At,1);

% If transformed_At is all-zero, just upsample it by rate of cs
if all(abs(transformed_At) == 0)
    At = zeros(Nt*cs,1);
    return
end

if cs < 1
    error('narrowband_scaledFTError:incorrect_cs',...
          'cs must be larger than or equal to 1.');
end

[~,peak_position] = max(abs(transformed_At).^2,[],1);
avg_center = sum((1:Nt)'.*abs(transformed_At).^2)/sum(abs(transformed_At).^2);
% If two positions are too far away, then it might indicate that the pulse
% is center at the left edge of the time window such that its other half is
% at the right edge of the window, due to periodic boundary condition under
% numerical DFT.
% If two positions are close, then we pick the avg_center as the pulse
% center.
if abs(peak_position - avg_center) > Nt/4
    fftshift_transformed_At = fftshift(transformed_At);
    avg_center = sum((1:Nt)'.*abs(fftshift_transformed_At).^2)/sum(abs(fftshift_transformed_At).^2);
    if avg_center >= floor(Nt/2)+1
        avg_center = avg_center - floor(Nt/2);
    else
        avg_center = avg_center + floor(Nt/2);
    end
end
pulse_center = avg_center;

transformed_phase_shift = ifftshift(2*pi*(1:Nt)'/Nt*pulse_center,1); % phase shift due to temporal offset

transformed_Aw = ifft(transformed_At,[],1);

transformed_phase = ifftshift(unwrap(angle(fftshift(transformed_Aw.*exp(-1i*transformed_phase_shift),1))),1);
transformed_phase = transformed_phase - transformed_phase(1,:,:);

recover_idx = linspace(1,Nt+1,Nt*cs+1)';
recover_idx = recover_idx(1:end-1); % remove the last "Nt+1" point
if mod(Nt,2) == 1
    recover_idx = [recover_idx(end-ceil((cs-1)/2)+1:end)-Nt;recover_idx(1:end-ceil((cs-1)/2))];
end

phase = ifftshift(interp1((1:Nt)',fftshift(transformed_phase,1),recover_idx,'linear','extrap')*cs,1);
abs_Aw = ifftshift(interp1((1:Nt)',abs(fftshift(transformed_Aw,1)),recover_idx,'linear','extrap'),1);
phase_shift = ifftshift(2*pi*(1:Nt*cs)'/(Nt*cs)*(pulse_center*cs),1); % phase shift due to temporal offset
At = fft(abs_Aw.*exp(1i*(phase+phase_shift))/sqrt(cs),[],1);

end