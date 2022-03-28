function [A, At, Gw, scale] = op_nufft_calib(p, D, N, Nn, No, Ns, S)
    % Create the nonuniform gridding matrix and fft operators.
    %
    % Parameters
    % ----------
    % p : double[:, 2]
    %     Nonuniformly distributed frequency location points.
    % D : complex[:, :]
    %     DDE kernels (in the Fourier domain).
    % N : int[2]
    %     Size of the reconstruction image.
    % Nn : int[2]
    %     Size of the kernels (number of neighbors considered on each direction).
    % No : int[2]
    %     Oversampled fft from which to recover the non uniform fft via kernel
    %     convolution.
    % Ns : int[2]
    %     FFT shift.
    % S : int[1]
    %     Size of the DDE support (square support).
    %
    % Returns
    % -------
    % A: function handle
    %     Function handle for direct operator.
    % At : function handle
    %     Function handle for adjoint operator.
    % Gw : sparse complex[:, :]
    %     Global convolution kernel matrix.
    % scale : complex[:, :]
    %     Scale paremters for the oversampled FFT.
    %

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
