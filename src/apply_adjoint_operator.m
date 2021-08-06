function Fx = apply_adjoint_operator(y, G, Sigma)
    
Fx = G' * (Sigma .* y);

end