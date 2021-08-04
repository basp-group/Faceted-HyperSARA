function y = apply_direct_operator(Fx, G, varargin)

if isempty(varargin) % no DR
    y = G * Fx;
else % DR
    % Sigma = varargin{1};
    y = varargin{1} .* (G * Fx);
end 

end