% Copyright (c) 2026 Elias S. G. Carotti

function state = ofdm_equalizer_init_state()
% OFDM_EQUALIZER_INIT_STATE  Initialize packet-local equalizer state.
    state = struct();
    state.residual_curve_est = [];
end
