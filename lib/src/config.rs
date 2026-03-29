// Copyright (c) 2026 Elias S. G. Carotti

/// Constellation mapping used on OFDM data carriers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Modulation {
    /// Binary phase-shift keying, one bit per symbol.
    Bpsk,
    /// Quadrature phase-shift keying, two bits per symbol.
    Qpsk,
}

/// Bitfield describing which equalizer stages are enabled.
///
/// Rationale:
/// The equalizer is now a pipeline of optional stages rather than a single
/// monolithic mode. This bitfield lets the configuration express combinations
/// such as:
/// - training baseline + pilot phase
/// - training baseline + pilot phase + pilot amplitude
/// - pilot-only equalization
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EqualizerFeatures(u32);

impl Default for EqualizerFeatures {
    fn default() -> Self {
        Self::NONE
    }
}

impl EqualizerFeatures {
    /// No optional equalizer stages.
    pub const NONE: Self = Self(0);
    /// Use the training symbol as the baseline channel model.
    pub const TRAINING_BASELINE: Self = Self(1 << 0);
    /// Apply a pilot-derived residual phase fit on each data symbol.
    pub const PILOT_PHASE: Self = Self(1 << 1);
    /// Apply a pilot-derived residual amplitude fit on each data symbol.
    pub const PILOT_AMPLITUDE: Self = Self(1 << 2);
    /// Weight pilot observations by reliability when fitting corrections.
    pub const WEIGHTED_PILOTS: Self = Self(1 << 3);
    /// Use a noise-aware MMSE-style inverse instead of the legacy heuristic one.
    pub const NOISE_AWARE_MMSE: Self = Self(1 << 4);
    /// Use temporal least-squares tracking on pilot-derived phase-line parameters.
    pub const TEMPORAL_LS: Self = Self(1 << 5);

    /// Returns whether all requested feature bits are enabled.
    pub fn contains(self, other: Self) -> bool {
        (self.0 & other.0) == other.0
    }
}

impl std::ops::BitOr for EqualizerFeatures {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self::Output {
        Self(self.0 | rhs.0)
    }
}

impl std::ops::BitOrAssign for EqualizerFeatures {
    fn bitor_assign(&mut self, rhs: Self) {
        self.0 |= rhs.0;
    }
}

/// Equalizer subsystem configuration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EqualizerConfig {
    /// Enabled equalizer stages.
    pub features: EqualizerFeatures,
    /// Number of recent symbols used by temporal least-squares tracking.
    pub temporal_window: usize,
}

impl EqualizerConfig {
    /// Starts a builder for an equalizer configuration.
    pub fn builder() -> EqualizerBuilder {
        EqualizerBuilder::default()
    }
}

impl Default for EqualizerConfig {
    fn default() -> Self {
        EqualizerConfig::builder()
            .training_baseline()
            .pilot_phase()
            .pilot_amplitude()
            .weighted_pilots()
            .build()
    }
}

/// Builder for [`EqualizerConfig`].
///
/// Rationale:
/// The equalizer now consists of composable stages. The builder keeps call
/// sites readable while still compiling down to simple feature checks.
#[derive(Clone, Copy, Debug, Default)]
pub struct EqualizerBuilder {
    features: EqualizerFeatures,
    temporal_window: usize,
}

impl EqualizerBuilder {
    /// Enable training-symbol baseline equalization.
    pub fn training_baseline(mut self) -> Self {
        self.features |= EqualizerFeatures::TRAINING_BASELINE;
        self
    }

    /// Enable pilot-derived residual phase correction.
    pub fn pilot_phase(mut self) -> Self {
        self.features |= EqualizerFeatures::PILOT_PHASE;
        self
    }

    /// Enable pilot-derived residual amplitude correction.
    pub fn pilot_amplitude(mut self) -> Self {
        self.features |= EqualizerFeatures::PILOT_AMPLITUDE;
        self
    }

    /// Enable reliability-weighted pilot fitting.
    pub fn weighted_pilots(mut self) -> Self {
        self.features |= EqualizerFeatures::WEIGHTED_PILOTS;
        self
    }

    /// Enable noise-aware MMSE regularization in the baseline equalizer.
    pub fn noise_aware_mmse(mut self) -> Self {
        self.features |= EqualizerFeatures::NOISE_AWARE_MMSE;
        self
    }

    /// Enable temporal least-squares tracking over the last `window` symbols.
    pub fn temporal_ls(mut self, window: usize) -> Self {
        self.features |= EqualizerFeatures::TEMPORAL_LS;
        self.temporal_window = window.max(1);
        self
    }

    /// Finalize the equalizer configuration.
    pub fn build(self) -> EqualizerConfig {
        EqualizerConfig {
            features: self.features,
            temporal_window: self.temporal_window.max(1),
        }
    }
}

/// Forward-error-correction scheme applied to packet bits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FecMode {
    /// No channel coding.
    None,
    /// Hamming(7,4) block coding on the serialized packet bitstream.
    Hamming74,
}

/// Passband implementation used to reach the speaker/microphone path.
///
/// Rationale:
/// `Legacy` keeps the modem at the audio rate, while `Iq` uses a separate
/// complex-baseband rate with explicit resampling and IQ conversion.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PassbandMode {
    /// Single-rate legacy path with the modem running directly at the audio rate.
    Legacy,
    /// IQ path with separate baseband/audio rates and explicit resampling.
    Iq,
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

/// Wake-up / preamble family transmitted before the OFDM body.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WakePreamble {
    /// Single tone wake preamble.
    Tone,
    /// Linear chirp wake preamble.
    Chirp,
    /// PN-sequence wake preamble.
    Pn,
    /// Gold-like wake preamble.
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
    /// Audio-side sample rate in hertz.
    pub fs: f32,
    /// Complex-baseband sample rate in hertz before IQ up/downsampling.
    pub fs_baseband: f32,
    /// Passband carrier frequency in hertz.
    pub fc: f32,
    /// Additional payload/body gain applied before final packet shaping.
    pub payload_gain: f32,
    /// IFFT/FFT size for OFDM symbols.
    pub nfft: usize,
    /// Cyclic-prefix length in samples at the active baseband rate.
    pub ncp: usize,
    /// Optional explicit baseband frequency origin override.
    pub base_freq_hz: Option<f32>,
    /// Active FFT-bin indices used by the modem.
    pub used_bins: Vec<usize>,
    /// Active FFT-bin indices reserved for pilots.
    pub pilot_bins: Vec<usize>,
    /// Optional legacy pilot-count override kept for compatibility.
    pub num_pilots: Option<usize>,
    /// Optional legacy pilot enable flag kept for compatibility.
    pub use_pilots: Option<bool>,
    /// Optional periodic retraining interval in data symbols.
    pub retrain_interval_data_symbols: Option<usize>,
    /// Whether to append a terminal training symbol at the end of the packet.
    pub terminal_training_symbol: bool,
    /// Data modulation used on the active data carriers.
    pub modulation: Modulation,
    /// Equalizer subsystem configuration.
    pub equalizer: EqualizerConfig,
    /// Forward-error-correction mode applied to packet bits.
    pub fec_mode: FecMode,
    /// Passband conversion path: direct legacy path or IQ path.
    pub passband_mode: PassbandMode,
    /// Wake-preamble duration in milliseconds.
    pub wake_ms: f32,
    /// Tone wake frequency in hertz when tone wake-up is enabled.
    pub wake_freq: f32,
    /// Silence/guard interval inserted after the wake preamble, in milliseconds.
    pub wake_guard_ms: f32,
    /// Wake-preamble family used before the OFDM body.
    pub wake_preamble: WakePreamble,
    /// Start frequency of the chirp wake preamble in hertz.
    pub sync_chirp_f0: f32,
    /// End frequency of the chirp wake preamble in hertz.
    pub sync_chirp_f1: f32,
    /// Half-length of the repeated-half sync sequence, in baseband samples.
    pub sync_half_len: usize,
    /// Maximum application payload bytes per packet fragment before FEC.
    pub packet_payload_bytes: usize,
    /// Session identifier written into packet headers.
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
            fs: 44_100.0,
            fs_baseband: 44_100.0,
            fc: 7_500.0,
            payload_gain: 2.0,
            nfft: 2048,
            ncp: 1024,
            base_freq_hz: None,
            used_bins: vec![
                24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
            ],
            pilot_bins: vec![24, 26, 28, 31, 33, 35, 38],
            num_pilots: None,
            use_pilots: Some(true),
            retrain_interval_data_symbols: None,
            terminal_training_symbol: false,
            modulation: Modulation::Bpsk,
            equalizer: EqualizerConfig::default(),
            fec_mode: FecMode::None,
            passband_mode: PassbandMode::Legacy,
            wake_ms: 80.0,
            wake_freq: 5_500.0,
            wake_guard_ms: 20.0,
            wake_preamble: WakePreamble::Tone,
            sync_chirp_f0: 4_500.0,
            sync_chirp_f1: 6_500.0,
            sync_half_len: 2048,
            packet_payload_bytes: 24,
            session_id: 1234,
        }
    }
}

// vim: set ts=4 sw=4 et:
