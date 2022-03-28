function [field1, field2, tx, ty] = stationary_Gaussian_process(m, n, rho)
    % Draw a realization of a stationary random Gaussian field over an
    % :math`m \times n` grid, following :cite:p:`Kroese2015`.
    %
    % Parameters
    % ----------
    % m : int
    %     Image dimension (y axis).
    % n : int
    %     Image dimension (x axis).
    % rho : anonymous function
    %     Function handle, function of a matrix :math:`h` such that
    %     :math:`\text{cov}(X_t, Y_s) = \rho(t-s)` is the covariance of
    %     a 2-dimensional stationary Gaussian field.
    %
    % Returns
    % -------
    % field1, field2 : double[:, :]
    %     Two independent realization over the :math`m \times n` grid
    %     considered.
    % tx, ty : double[:]
    %     Two vectors such that each field can be displayed using the command
    %     ``imagesc(tx, ty, field1)`` (idem for field2).
    %
    % Example
    % -------
    % .. code-block: matlab
    %
    %    % define covariance function
    %    rho=@(h)((1-h(1)^2/50^2-h(1)*h(2)/(15*50)-h(2)^2/15^2)...
    %              *exp(-(h(1)^2/50^2+h(2)^2/15^2)));
    %    [field1, field2, tx, ty] = stationary_Gaussian_process(512,384,rho);
    %    imagesc(tx, ty, field1);
    %
    % Note
    % ----
    % - The size of covariance matrix is ``m^2*n^2``.
    % - By default, when the function is called without assigning the output,
    %   `field1` is drectly display using the ``imagesc`` MATLAB function.
    %

    %% Reference:
    % Kroese, D. P., & Botev, Z. I. (2015). Spatial Process Simulation.
    % In Stochastic Geometry, Spatial Statistics and Random Fields(pp. 369-404)
    % Springer International Publishing, DOI: 10.1007/978-3-319-10064-7_12

    % create grid for field
    tx = 0:(n - 1);
    ty = 0:(m - 1);
    Rows = zeros(m, n);
    Cols = Rows;
    for i = 1:n % sample covariance function at grid points;
        for j = 1:m
            % rows of blocks of cov matrix
            Rows(j, i) = rho([tx(i) - tx(1), ty(j) - ty(1)]);
            % columns of blocks of cov matrix
            Cols(j, i) = rho([tx(1) - tx(i), ty(j) - ty(1)]);
        end
    end
    % create the first row of the block circulant matrix with circular blocks
    % and store it as a matrix suitable for fft2;
    BlkCirc_row = [Rows, Cols(:, end:-1:2)
                   Cols(end:-1:2, :), Rows(end:-1:2, end:-1:2)];
    % compute eigen-values
    lam = real(fft2(BlkCirc_row)) / (2 * m - 1) / (2 * n - 1);
    if abs(min(lam(lam(:) < 0))) > 10^-15
        error('Could not find positive definite embedding!');
    else
        lam(lam(:) < 0) = 0;
        lam = sqrt(lam);
    end
    % generate field with covariance given by block circular matrix
    F = fft2(lam .* complex(randn(2 * m - 1, 2 * n - 1), randn(2 * m - 1, 2 * n - 1)));
    F = F(1:m, 1:n); % extract subblock with desired covariance
    field1 = real(F);
    field2 = imag(F); % two independent fields with desired covariance
    if nargout == 0
        imagesc(tx, ty, field1);
        colormap bone;
    end
