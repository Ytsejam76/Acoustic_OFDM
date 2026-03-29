// Copyright (c) 2026 Elias S. G. Carotti

use rustfft::num_complex::Complex32;

use crate::config::OfdmConfig;
use crate::modem;

/// Packet-level metadata produced during oracle encoding.
///
/// Rationale:
/// - The oracle path needs exact packet boundaries and the corresponding
///   complex baseband samples in order to separate modem failures from search
///   failures.
#[derive(Clone, Debug)]
pub struct EncodedPacketMeta {
    /// Zero-based fragment index inside the encoded burst.
    pub frag_index: usize,
    /// Total number of fragments in the encoded burst.
    pub frag_count: usize,
    /// Packet start offset, in passband audio samples, inside [`EncodedBurst::audio`].
    pub packet_start: usize,
    /// Packet length, in passband audio samples.
    pub packet_len: usize,
    /// Complex baseband packet waveform before passband conversion.
    pub xbb: Vec<Complex32>,
}

/// Full encoded burst with passband audio plus oracle packet metadata.
#[derive(Clone, Debug)]
pub struct EncodedBurst {
    /// Complete burst waveform at the audio sample rate.
    pub audio: Vec<f32>,
    /// Per-packet metadata for oracle decode and debugging.
    pub packet_meta: Vec<EncodedPacketMeta>,
}

/// Summary diagnostics for one decoded passband packet window.
///
/// Rationale:
/// - This is the compact debug summary used by the CLI to decide whether a run
///   failed because of timing, CFO, channel estimation, or payload demodulation.
#[derive(Clone, Debug)]
pub struct PassbandDiagnostics {
    /// Whether the packet window contained enough samples for sync, training,
    /// and at least one OFDM data symbol.
    pub enough_samples: bool,
    /// Chosen sync offset, in baseband samples, relative to the start of the
    /// post-wake/post-guard window.
    pub sync_off: usize,
    /// Coarse carrier-frequency offset estimate in hertz.
    pub cfo_hz: f32,
    /// RMS magnitude over the repeated-half sync region.
    pub sync_rms: f32,
    /// Peak magnitude over the repeated-half sync region.
    pub sync_peak: f32,
    /// RMS magnitude over the post-sync OFDM payload region.
    pub post_rms: f32,
    /// Peak magnitude over the post-sync OFDM payload region.
    pub post_peak: f32,
    /// RMS magnitude over the training symbol after CP removal.
    pub train_rms: f32,
    /// Minimum channel-estimate magnitude observed across used bins.
    pub hest_mag_min: f32,
    /// Mean channel-estimate magnitude observed across used bins.
    pub hest_mag_mean: f32,
    /// Maximum channel-estimate magnitude observed across used bins.
    pub hest_mag_max: f32,
    /// Training-symbol reconstruction EVM after equalization.
    pub train_recon_evm: f32,
    /// Residual pilot EVM before the final pilot-based cleanup.
    pub pilot_residual_evm: f32,
    /// Pilot EVM after the final pilot-based cleanup.
    pub pilot_post_evm: f32,
    /// Decision-directed EVM over data carriers after equalization.
    pub post_eq_evm: f32,
    /// Whether the packet parser accepted the decoded packet.
    pub decoded: bool,
    /// Decoded payload length when parsing succeeded.
    pub decoded_payload_len: Option<usize>,
}

/// Equalizer input/output samples for a constellation plot.
#[derive(Clone, Debug)]
pub struct PassbandConstellationDump {
    /// Raw FFT-bin samples on data carriers before equalization.
    pub pre_eq: Vec<Complex32>,
    /// Equalized data-carrier samples after pilot/training correction.
    pub post_eq: Vec<Complex32>,
}

/// Per-symbol pilot-tracking summary for one packet.
#[derive(Clone, Debug)]
pub struct PassbandPilotTrackDump {
    /// Residual pilot phase estimate per data symbol, in radians.
    pub pilot_phase_rad: Vec<f32>,
    /// Pilot EVM before the final pilot cleanup.
    pub pilot_evm_pre: Vec<f32>,
    /// Pilot EVM after the final pilot cleanup.
    pub pilot_evm_post: Vec<f32>,
    /// Mean magnitude of the current channel estimate per tracked symbol.
    pub hest_mag_mean: Vec<f32>,
    /// Maximum magnitude of the current channel estimate per tracked symbol.
    pub hest_mag_max: Vec<f32>,
}

/// One per-bin debug row for the passband bin dump.
#[derive(Clone, Debug)]
pub struct PassbandBinDumpRow {
    /// One-based OFDM data-symbol index inside the packet payload section.
    pub data_symbol_idx: usize,
    /// FFT bin index used by this carrier.
    pub used_bin: usize,
    /// Carrier role: `"pilot"` or `"data"`.
    pub role: &'static str,
    /// Complex carrier value before equalization.
    pub pre_eq: Complex32,
    /// Complex carrier value after equalization.
    pub post_eq: Complex32,
    /// Known reference symbol when available, mainly for pilot carriers.
    pub reference: Option<Complex32>,
}

/// Full per-bin dump for a few decoded OFDM symbols.
#[derive(Clone, Debug)]
pub struct PassbandBinDump {
    /// Flat row list across symbols and used bins.
    pub rows: Vec<PassbandBinDumpRow>,
}

/// One row of the channel-comparison dump.
#[derive(Clone, Debug)]
pub struct PassbandChannelCompareRow {
    /// One-based OFDM data-symbol index inside the packet payload section.
    pub data_symbol_idx: usize,
    /// FFT bin index used by this carrier.
    pub used_bin: usize,
    /// Carrier role: `"pilot"` or `"data"`.
    pub role: &'static str,
    /// Actual per-bin channel computed from known transmitted and received bins.
    pub actual_h: Complex32,
    /// Training-only channel estimate used as the baseline equalizer model.
    pub estimated_h_train: Complex32,
    /// Pilot-adjusted channel estimate after residual per-symbol correction.
    pub estimated_h_pilot: Complex32,
}

/// Full channel-comparison dump for a few OFDM symbols.
#[derive(Clone, Debug)]
pub struct PassbandChannelCompareDump {
    /// Flat row list across symbols and used bins.
    pub rows: Vec<PassbandChannelCompareRow>,
}

/// Timing metric dump around the repeated-half synchronizer.
#[derive(Clone, Debug)]
pub struct PassbandSyncDump {
    /// Coarse Schmidl-Cox timing estimate.
    pub coarse_sync_off: usize,
    /// Refined timing estimate after the local refinement step.
    pub refined_sync_off: usize,
    /// Schmidl-Cox metric trace used for visualization/debugging.
    pub metrics: Vec<f32>,
}

/// IQ-chain dump used to inspect the RX downconversion path.
#[derive(Clone, Debug)]
pub struct PassbandIqChainDump {
    /// Complex signal after IQ downconversion and low-pass filtering at the audio rate.
    pub downconverted_audio_rate: Vec<Complex32>,
    /// Complex signal after optional resampling to the baseband rate.
    pub baseband_rate: Vec<Complex32>,
    /// Audio sample rate in hertz.
    pub fs_audio: f32,
    /// Decoder baseband sample rate in hertz.
    pub fs_baseband: f32,
}

/// Extract equalizer input/output constellation samples for a packet window.
pub fn dump_passband_constellation(
    pkt_audio: &[f32],
    cfg: &OfdmConfig,
    sync_off: f32,
) -> Option<PassbandConstellationDump> {
    modem::dump_passband_constellation_impl(pkt_audio, cfg, sync_off)
}

/// Extract per-symbol pilot tracking diagnostics for a packet window.
pub fn dump_passband_pilot_tracking(
    pkt_audio: &[f32],
    cfg: &OfdmConfig,
) -> Option<PassbandPilotTrackDump> {
    modem::dump_passband_pilot_tracking_impl(pkt_audio, cfg)
}

/// Extract pre/post-equalization carrier values for a packet window.
pub fn dump_passband_bins(pkt_audio: &[f32], cfg: &OfdmConfig) -> Option<PassbandBinDump> {
    modem::dump_passband_bins_impl(pkt_audio, cfg)
}

/// Extract pre/post-equalization carrier values using a known sync offset.
pub fn dump_passband_bins_with_sync(
    pkt_audio: &[f32],
    cfg: &OfdmConfig,
    sync_off: f32,
) -> Option<PassbandBinDump> {
    modem::dump_passband_bins_with_sync_impl(pkt_audio, cfg, sync_off)
}

/// Compare actual and estimated channel samples using a known transmitted packet.
pub fn dump_passband_channel_compare_with_sync(
    payload: &[u8],
    pkt_audio: &[f32],
    cfg: &OfdmConfig,
    sync_off: f32,
) -> Option<PassbandChannelCompareDump> {
    modem::dump_passband_channel_compare_with_sync_impl(payload, pkt_audio, cfg, sync_off)
}

/// Extract Schmidl-Cox timing metrics for one packet window.
pub fn dump_passband_sync_metric(pkt_audio: &[f32], cfg: &OfdmConfig) -> Option<PassbandSyncDump> {
    modem::dump_passband_sync_metric_impl(pkt_audio, cfg)
}

/// Dump the RX IQ chain before and after resampling.
pub fn dump_passband_iq_chain(pkt_audio: &[f32], cfg: &OfdmConfig) -> Option<PassbandIqChainDump> {
    modem::dump_passband_iq_chain_impl(pkt_audio, cfg)
}
