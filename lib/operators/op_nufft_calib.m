function [A, At, Gw, scale] = op_nufft_calib(p, D, N, Nn, No, Ns, S)
    % Create the nonuniform gridding matrix and fft operators
    %
    % in:
    % p[2]    - nonuniformly distributed frequency location points
    % D[2]    - DDE kernels (in the Fourier domain)
    % N[2]    - size of the reconstruction image
    % Nn[2]   - size of the kernels (number of neighbors considered on each direction)
    % No[2]   - oversampled fft from which to recover the non uniform fft via
    %           kernel convolution
    % Ns[2]   - fft shift
    % S[1]    - size of the DDE support (square support)
    %
    % out:
    % A[@]          - function handle for direct operator
    % At[@]         - function handle for adjoint operator
    % Gw[:][:]      - global convolution kernel matrix
    % scale[:][:]   - scale paremters for the oversampled FFT

    % [06/08/2019] Code modification, P.-A. Thouvenin.

    %% compute the overall gridding matrix and its associated kernels

    % st = nufft_init(p, N, Nn, No, Ns);
    % scale = st.sn;

    % [Gw, scale] = createG(p, Nn, N, No, Ns);
    [Gw, scale] = createG_calib(D, p, N, No, Nn(1), S, Ns);

    % if isa(st.p, 'fatrix2')
    %     error('fatrix2 has some very weird bugs with subindexing; force st.p to be a (sparse) matrix');
    % end

    %% operator function handles
    [A, At] = op_nu_so_fft2(N, No, scale);

    % whole G is stored in st.p
    % Gw = st.p;

end
