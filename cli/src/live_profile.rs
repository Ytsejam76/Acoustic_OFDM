// Copyright (c) 2026 Elias S. G. Carotti

use clap::ValueEnum;

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
pub enum LiveProfileArg {
    Standard,
    #[value(name = "live-debug")]
    LiveDebug,
}

pub fn tx_default_log_file(profile: LiveProfileArg) -> Option<String> {
    match profile {
        LiveProfileArg::Standard => None,
        LiveProfileArg::LiveDebug => Some("/tmp/acoustic_ofdm_tx.log".to_string()),
    }
}

pub fn rx_default_log_file(profile: LiveProfileArg) -> Option<String> {
    match profile {
        LiveProfileArg::Standard => None,
        LiveProfileArg::LiveDebug => Some("/tmp/acoustic_ofdm_rx.log".to_string()),
    }
}

// vim: set ts=4 sw=4 et:
