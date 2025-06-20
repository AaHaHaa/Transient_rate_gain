function output = saturable_absorber_action_simple(input, saturation_power, moddepth)
%SATURABLE_ABSORBER_ACTION_SIMPLE Apply an ideal saturable absorber effect to each mode individually
%
% input - either a struct or a matrix:
%   input.fields - a (N, num_modes, m) matrix with the fields for each mode, in the time domain
%
% saturation_power - the scale power for the saturable absorber, in W
% moddepth - the modulation depth, 0-1

if isstruct(input)
    input_field = input.fields(:, :, end);
else
    error('saturable_sbsorber_action_simple:inputError',...
          '"input" must be a structure with "dt" and "fields".');
end

powers = sum(abs( input_field ).^2,2);
%trans_field = input_field.*sqrt(1 - moddepth./(1 + powers/saturation_power));
%reflect_field = input_field.*sqrt(moddepth./(1 + powers/saturation_power));
trans_fields = input_field.*sqrt(1 - moddepth*exp(-powers/saturation_power));
reflect_fields = input_field.*sqrt(moddepth*exp(-powers/saturation_power));

output = struct('fields',           trans_fields,...
                'reflect_fields', reflect_fields);
if isfield(input,'dt')
    output.dt = input.dt;
end

end