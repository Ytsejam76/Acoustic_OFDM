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

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WakePreamble {
    Tone,
    Chirp,
    Pn,
    Gold,
}

impl WakePreamble {
    /// Parses a wake preamble mode from CLI/config text.
    ///
    /// Parameters:
    /// - `s`: mode string.
    /// Returns:
    /// - `Option<WakePreamble>`: parsed mode when recognized.
    pub fn parse(s: &str) -> Option<Self> {
        match s.to_ascii_lowercase().as_str() {
            "tone" => Some(Self::Tone),
            "chirp" => Some(Self::Chirp),
            "pn" => Some(Self::Pn),
            "gold" => Some(Self::Gold),
            _ => None,
        }
    }

    /// Returns the canonical mode name.
    ///
    /// Parameters:
    /// - `self`: wake preamble variant.
    /// Returns:
    /// - `&'static str`: printable mode name.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Tone => "tone",
            Self::Chirp => "chirp",
            Self::Pn => "pn",
            Self::Gold => "gold",
        }
    }
}

#[derive(Clone, Debug)]
pub struct OfdmConfig {
    pub fs: f32,
    pub fc: f32,
    pub nfft: usize,
    pub ncp: usize,
    pub base_freq_hz: Option<f32>,
    pub used_bins: Vec<usize>,
    pub pilot_bins: Vec<usize>,
    pub num_pilots: Option<usize>,
    pub use_pilots: Option<bool>,
    pub retrain_interval_data_symbols: Option<usize>,
    pub terminal_training_symbol: bool,
    pub modulation: Modulation,
    pub wake_ms: f32,
    pub wake_freq: f32,
    pub wake_guard_ms: f32,
    pub wake_preamble: WakePreamble,
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
            fc: 14_000.0,
            nfft: 96,
            ncp: 72,
            base_freq_hz: None,
            used_bins: vec![2, 3, 4, 5],
            pilot_bins: vec![2, 5],
            num_pilots: None,
            use_pilots: None,
            retrain_interval_data_symbols: Some(1),
            terminal_training_symbol: true,
            modulation: Modulation::Bpsk,
            wake_ms: 100.0,
            wake_freq: 12_000.0,
            wake_guard_ms: 30.0,
            wake_preamble: WakePreamble::Gold,
            sync_chirp_f0: 4_000.0,
            sync_chirp_f1: 8_000.0,
            sync_half_len: 6_000,
            packet_payload_bytes: 24,
            session_id: 1234,
        }
    }
}

// vim: set ts=4 sw=4 et:
