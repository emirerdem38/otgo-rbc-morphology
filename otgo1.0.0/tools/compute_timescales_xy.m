function results = compute_timescales_xy(D, k_trans_xy, k_rot_xy)
% COMPUTE_TIMESCALES_XY
% Computes translational and rotational relaxation times
% using only x- and y-axis stiffness values.
%
% Inputs:
%   D           : 6x6 diffusion tensor in SI units
%   k_trans_xy  : [kx, ky] translational stiffness in pN/um
%   k_rot_xy    : [kappa_x, kappa_y] rotational stiffness in pN*um/rad
%                 use NaN for undefined modes
%
% Outputs:
%   results     : struct with friction coefficients and relaxation times

    kB = PhysConst.kB;   % J/K
    T  = 293;            % K

    % Diffusion coefficients from diagonal of D
    % Translational: rows/cols 1:3 -> x,y,z
    % Rotational:    rows/cols 4:6 -> rx,ry,rz
    Dtx = D(1,1);
    Dty = D(2,2);
    Drx = D(4,4);
    Dry = D(5,5);

    % Friction coefficients from fluctuation-dissipation relation
    gamma_tx = kB * T / Dtx;
    gamma_ty = kB * T / Dty;
    gamma_rx = kB * T / Drx;
    gamma_ry = kB * T / Dry;

    % Convert stiffness units:
    % pN/um -> N/m
    kx = k_trans_xy(1) * 1e-6;
    ky = k_trans_xy(2) * 1e-6;

    % pN*um/rad -> N*m/rad
    kappax = k_rot_xy(1);
    kappay = k_rot_xy(2);

    if ~isnan(kappax)
        kappax = kappax * 1e-18;
    end
    if ~isnan(kappay)
        kappay = kappay * 1e-18;
    end

    % Translational relaxation times
    tau_tx = gamma_tx / kx;
    tau_ty = gamma_ty / ky;

    % Rotational relaxation times
    tau_rx = NaN;
    tau_ry = NaN;

    if ~isnan(kappax) && kappax > 0
        tau_rx = gamma_rx / kappax;
    end

    if ~isnan(kappay) && kappay > 0
        tau_ry = gamma_ry / kappay;
    end

    % Store outputs
    results.Dtx = Dtx;
    results.Dty = Dty;
    results.Drx = Drx;
    results.Dry = Dry;

    results.gamma_tx = gamma_tx;
    results.gamma_ty = gamma_ty;
    results.gamma_rx = gamma_rx;
    results.gamma_ry = gamma_ry;

    results.tau_tx = tau_tx;
    results.tau_ty = tau_ty;
    results.tau_rx = tau_rx;
    results.tau_ry = tau_ry;

    % Minimum relevant relaxation time among defined modes
    relevant_taus = [tau_tx, tau_ty, tau_rx, tau_ry];
    relevant_taus = relevant_taus(~isnan(relevant_taus) & isfinite(relevant_taus));

    if isempty(relevant_taus)
        results.tau_min = NaN;
    else
        results.tau_min = min(relevant_taus);
    end
end