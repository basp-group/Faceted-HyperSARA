function GFx = apply_direct_operator(Fx, G, varargin)

% try to avoid one operator per data block

% 

if isempty(varargin) % no DR
    for l = 1:size(Fx, 3)
        for b = 1:numel(G{l})
            y = G * Fx(W{i}{j}), G{i}{j})
            
        end
    end
else % DR
        T = varargin{1};
    end

end

    

end