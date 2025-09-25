function E_ext = compute_external_Efield(params, x, t)
% Compute external electric field for the simulation
%
% Inputs:
%   params - simulation parameters structure
%   x      - grid points where to evaluate the external field
%   t      - time at which to evaluate the external field
%
% Output:
%   E_ext  - external electric field values at grid points

    % Check if external field is defined
    if ~isfield(params, 'E_ext') || isempty(params.E_ext)
        E_ext = 0;
        return;
    end
    
    % Evaluate external field function
    E_ext = params.E_ext(x, t);
end
