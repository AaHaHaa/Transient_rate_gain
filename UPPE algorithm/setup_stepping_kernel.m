function [gpuDevice_Device,...
          cuda_SRSK,num_operations_SRSK] = setup_stepping_kernel(sim,Nt)
%SETUP_STEPPING_KERNEL It sets cuda for computing sums of SR and SK terms,
%and spontaneous Raman terms.

%% Use the specified GPU
% This needs to run at the beginning; otherwise, the already-stored values
% in GPU will be unavailable in a new GPU if the GPU device is switched.
try
    gpuDevice_Device = gpuDevice(sim.gpuDevice.Index); % use the specified GPU device
catch
    error('Please set the GPU you''re going to use by setting "sim.gpuDevice.Index".');
end

%% SR, SK
% Nonlinear term
SRSK_filename = 'GMMNLSE_nonlinear_sum';
num_operations_SRSK = 2;

% Kernal for computing the nonlinear term of UPPE
cuda_SRSK = setup_kernel_SRSK(SRSK_filename,sim.cuda_dir_path,Nt,1,num_operations_SRSK,1);

end

