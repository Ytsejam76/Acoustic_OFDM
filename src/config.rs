// Copyright (c) 2026 Elias S. G. Carotti

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Modulation {
    Bpsk,
    Qpsk,
}

impl Modulation {
    /// Returns the number of bits carried by one constellation symbol.
    ///
    /// Parameters:
    /// - `self`: modulation variant.
    /// Returns:
    /// - `usize`: bits per symbol (`1` for BPSK, `2` for QPSK).
    pub fn bits_per_symbol(self) -> usize {
        match self {
            Self::Bpsk => 1,
            Self::Qpsk => 2,
        }
    }

    /// Returns the packet header modulation identifier.
    ///
    /// Parameters:
    /// - `self`: modulation variant.
    /// Returns:
    /// - `u8`: modulation ID used in packet headers.
    pub fn mod_id(self) -> u8 {
        match self {
            Self::Bpsk => 1,
            Self::Qpsk => 2,
        }
    }
}

#[derive(Clone, Debug)]
pub struct OfdmConfig {
    pub fs: f32,
    pub fc: f32,
    pub nfft: usize,
    pub ncp: usize,
    pub used_bins: Vec<usize>,
    pub modulation: Modulation,
    pub wake_ms: f32,
    pub wake_freq: f32,
    pub wake_guard_ms: f32,
    pub use_chirp_sync: bool,
    pub sync_chirp_f0: f32,
    pub sync_chirp_f1: f32,
    pub sync_half_len: usize,
    pub packet_payload_bytes: usize,
    pub session_id: u16,
}

impl Default for OfdmConfig {
    /// Builds the default OFDM modem configuration.
    ///
    /// Parameters:
    /// - none.
    /// Returns:
    /// - `OfdmConfig`: default configuration values.
    fn default() -> Self {
        Self {
            fs: 48_000.0,
            fc: 17_000.0,
            nfft: 96,
            ncp: 72,
            used_bins: vec![2, 3, 4, 5],
            modulation: Modulation::Bpsk,
            wake_ms: 12.0,
            wake_freq: 16_500.0,
            wake_guard_ms: 4.0,
            use_chirp_sync: true,
            sync_chirp_f0: 4_000.0,
            sync_chirp_f1: 8_000.0,
            sync_half_len: 64,
            packet_payload_bytes: 24,
            session_id: 1234,
        }
    }
}

// vim: set ts=4 sw=4 et:
