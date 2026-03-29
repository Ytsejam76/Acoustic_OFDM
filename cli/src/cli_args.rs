// Copyright (c) 2026 Elias S. G. Carotti

use std::error::Error;
use std::io::Read;

use acoustic_ofdm::{
    EqualizerConfig, FecMode, Modulation, OfdmConfig, PassbandMode, SpectrogramOptions,
    SpectrogramWindow, WakePreamble,
};
use clap::{Args, Parser, Subcommand, ValueEnum};

use crate::live_profile::LiveProfileArg;
use crate::logging::LogLevelArg;

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
pub(crate) enum WakePreambleArg {
    Gold,
    Pn,
    Chirp,
    Tone,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
pub(crate) enum ModulationArg {
    Bpsk,
    Qpsk,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
pub(crate) enum EqualizerModeArg {
    TrainingPilot,
    PilotOnly,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
pub(crate) enum FecModeArg {
    None,
    Hamming74,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
pub(crate) enum PassbandModeArg {
    Legacy,
    Iq,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
pub(crate) enum SpectrogramWindowArg {
    Hann,
    Hamming,
    Blackman,
    Rect,
}

impl From<SpectrogramWindowArg> for SpectrogramWindow {
    fn from(value: SpectrogramWindowArg) -> Self {
        match value {
            SpectrogramWindowArg::Hann => SpectrogramWindow::Hann,
            SpectrogramWindowArg::Hamming => SpectrogramWindow::Hamming,
            SpectrogramWindowArg::Blackman => SpectrogramWindow::Blackman,
            SpectrogramWindowArg::Rect => SpectrogramWindow::Rect,
        }
    }
}

impl From<WakePreambleArg> for WakePreamble {
    fn from(value: WakePreambleArg) -> Self {
        match value {
            WakePreambleArg::Gold => WakePreamble::Gold,
            WakePreambleArg::Pn => WakePreamble::Pn,
            WakePreambleArg::Chirp => WakePreamble::Chirp,
            WakePreambleArg::Tone => WakePreamble::Tone,
        }
    }
}

impl From<ModulationArg> for Modulation {
    fn from(value: ModulationArg) -> Self {
        match value {
            ModulationArg::Bpsk => Modulation::Bpsk,
            ModulationArg::Qpsk => Modulation::Qpsk,
        }
    }
}

impl From<EqualizerModeArg> for EqualizerConfig {
    fn from(value: EqualizerModeArg) -> Self {
        match value {
            EqualizerModeArg::TrainingPilot => EqualizerConfig::builder()
                .training_baseline()
                .pilot_phase()
                .pilot_amplitude()
                .build(),
            EqualizerModeArg::PilotOnly => EqualizerConfig::builder()
                .pilot_phase()
                .pilot_amplitude()
                .build(),
        }
    }
}

impl From<FecModeArg> for FecMode {
    fn from(value: FecModeArg) -> Self {
        match value {
            FecModeArg::None => FecMode::None,
            FecModeArg::Hamming74 => FecMode::Hamming74,
        }
    }
}

impl From<PassbandModeArg> for PassbandMode {
    fn from(value: PassbandModeArg) -> Self {
        match value {
            PassbandModeArg::Legacy => PassbandMode::Legacy,
            PassbandModeArg::Iq => PassbandMode::Iq,
        }
    }
}

#[derive(Debug, Clone, Args, Default)]
pub(crate) struct CommonCfgArgs {
    #[arg(short = 'b', long)]
    pub(crate) base_freq_hz: Option<f32>,
    #[arg(long)]
    pub(crate) fs_baseband: Option<f32>,
    #[arg(long)]
    pub(crate) nfft: Option<usize>,
    #[arg(long)]
    pub(crate) ncp: Option<usize>,
    #[arg(long)]
    pub(crate) sync_half_len: Option<usize>,
    #[arg(short = 'm', long, value_enum)]
    pub(crate) modulation: Option<ModulationArg>,
    #[arg(long, value_enum)]
    pub(crate) equalizer_mode: Option<EqualizerModeArg>,
    #[arg(long, value_enum)]
    pub(crate) fec_mode: Option<FecModeArg>,
    #[arg(long, value_enum)]
    pub(crate) passband_mode: Option<PassbandModeArg>,
    #[arg(short = 'w', long, value_enum)]
    pub(crate) wake_preamble: Option<WakePreambleArg>,
}

#[derive(Debug, Parser)]
#[command(name = "acoustic_ofdm_cli")]
pub(crate) struct Cli {
    #[command(subcommand)]
    pub(crate) command: Commands,
}

#[derive(Debug, Subcommand)]
pub(crate) enum Commands {
    Encode(EncodeCmd),
    EncodeBody(EncodeCmd),
    Decode(DecodeCmd),
    MicRoundtrip(MicRoundtripCmd),
    Roundtrip(RoundtripCmd),
    Scan(ScanCmd),
    Spectrogram(SpectrogramCmd),
    Rx(RxCmd),
    Tx(TxCmd),
}

#[derive(Debug, Args)]
pub(crate) struct EncodeCmd {
    #[command(flatten)]
    pub(crate) common: CommonCfgArgs,
    pub(crate) out_wav: String,
    pub(crate) payload_text: String,
}

#[derive(Debug, Args)]
pub(crate) struct DecodeCmd {
    #[command(flatten)]
    pub(crate) common: CommonCfgArgs,
    #[arg(short = 'p', long, value_enum, default_value_t = LiveProfileArg::LiveDebug)]
    pub(crate) profile: LiveProfileArg,
    #[arg(short = 's', long)]
    pub(crate) stdout: bool,
    #[arg(long)]
    pub(crate) start_sec: Option<f32>,
    #[arg(long)]
    pub(crate) window_sec: Option<f32>,
    #[arg(long)]
    pub(crate) sync_off: Option<f32>,
    pub(crate) in_wav: String,
}

#[derive(Debug, Args)]
pub(crate) struct RoundtripCmd {
    #[command(flatten)]
    pub(crate) common: CommonCfgArgs,
    #[arg(short = 's', long)]
    pub(crate) stdout: bool,
    pub(crate) wav_path: String,
    pub(crate) payload_text: String,
}

#[derive(Debug, Args)]
pub(crate) struct MicRoundtripCmd {
    #[command(flatten)]
    pub(crate) common: CommonCfgArgs,
    #[arg(short = 'p', long, value_enum, default_value_t = LiveProfileArg::LiveDebug)]
    pub(crate) profile: LiveProfileArg,
    #[arg(short = 'd', long)]
    pub(crate) duration_sec: Option<f32>,
    #[arg(long)]
    pub(crate) pre_delay_sec: Option<f32>,
    #[arg(short = 'r', long)]
    pub(crate) repeats: Option<usize>,
    #[arg(short = 'G', long)]
    pub(crate) gap_sec: Option<f32>,
    #[arg(short = 'g', long)]
    pub(crate) mic_gain: Option<f32>,
    #[arg(short = 'a', long)]
    pub(crate) spk_gain: Option<f32>,
    #[arg(short = 'W', long)]
    pub(crate) dump_wav: Option<String>,
    #[arg(long)]
    pub(crate) dump_tx_wav: Option<String>,
    #[arg(short = 's', long)]
    pub(crate) spectrogram: bool,
    #[arg(short = 'P', long, default_value = "/tmp/rx_spectrogram.png")]
    pub(crate) spectrogram_path: String,
    #[arg(long, default_value_t = 512)]
    pub(crate) spectrogram_nfft: usize,
    #[arg(long, default_value_t = 128)]
    pub(crate) spectrogram_hop: usize,
    #[arg(long, value_enum, default_value_t = SpectrogramWindowArg::Hann)]
    pub(crate) spectrogram_window: SpectrogramWindowArg,
    #[arg(short = 'o', long)]
    pub(crate) oracle: bool,
    #[arg(short = 'S', long)]
    pub(crate) stdout: bool,
    #[arg(short = 'v', long)]
    pub(crate) verbose: bool,
    #[arg(short = 'l', long, value_enum)]
    pub(crate) log_level: Option<LogLevelArg>,
    #[arg(short = 'L', long)]
    pub(crate) log_file: Option<String>,
    pub(crate) payload_text: Option<String>,
}

#[derive(Debug, Args)]
pub(crate) struct ScanCmd {
    #[command(flatten)]
    pub(crate) common: CommonCfgArgs,
    #[arg(short = 'p', long, value_enum, default_value_t = LiveProfileArg::LiveDebug)]
    pub(crate) profile: LiveProfileArg,
    #[arg(short = 'i', long)]
    pub(crate) in_wav: String,
    #[arg(long, default_value_t = 0.0)]
    pub(crate) start_sec: f32,
    #[arg(long)]
    pub(crate) end_sec: Option<f32>,
    #[arg(long, default_value_t = 1.40)]
    pub(crate) window_sec: f32,
    #[arg(long, default_value_t = 20.0)]
    pub(crate) step_ms: f32,
    #[arg(long, default_value_t = 0)]
    pub(crate) sync_min: i32,
    #[arg(long, default_value_t = 80)]
    pub(crate) sync_max: i32,
    #[arg(long, default_value_t = 4)]
    pub(crate) sync_step: i32,
    #[arg(long, default_value_t = 8)]
    pub(crate) top_k: usize,
}

#[derive(Debug, Args)]
pub(crate) struct SpectrogramCmd {
    #[arg(short = 'i', long)]
    pub(crate) in_wav: String,
    #[arg(short = 'o', long)]
    pub(crate) out_png: String,
    #[arg(long, default_value_t = 512)]
    pub(crate) spectrogram_nfft: usize,
    #[arg(long, default_value_t = 128)]
    pub(crate) spectrogram_hop: usize,
    #[arg(long, value_enum, default_value_t = SpectrogramWindowArg::Hann)]
    pub(crate) spectrogram_window: SpectrogramWindowArg,
}

#[derive(Debug, Args)]
pub(crate) struct TxCmd {
    #[command(flatten)]
    pub(crate) common: CommonCfgArgs,
    #[arg(short = 'p', long, value_enum, default_value_t = LiveProfileArg::LiveDebug)]
    pub(crate) profile: LiveProfileArg,
    #[arg(short = 'g', long)]
    pub(crate) spk_gain: Option<f32>,
    #[arg(short = 'd', long)]
    pub(crate) pre_delay_sec: Option<f32>,
    #[arg(short = 'r', long)]
    pub(crate) repeats: Option<usize>,
    #[arg(short = 'G', long)]
    pub(crate) gap_sec: Option<f32>,
    #[arg(short = 'W', long)]
    pub(crate) dump_wav: Option<String>,
    #[arg(short = 'o', long)]
    pub(crate) oracle: bool,
    #[arg(short = 'v', long)]
    pub(crate) verbose: bool,
    #[arg(short = 'l', long, value_enum)]
    pub(crate) log_level: Option<LogLevelArg>,
    #[arg(short = 'L', long)]
    pub(crate) log_file: Option<String>,
    pub(crate) payload_text: Option<String>,
}

#[derive(Debug, Args)]
pub(crate) struct RxCmd {
    #[command(flatten)]
    pub(crate) common: CommonCfgArgs,
    #[arg(short = 'p', long, value_enum, default_value_t = LiveProfileArg::LiveDebug)]
    pub(crate) profile: LiveProfileArg,
    #[arg(short = 'd', long)]
    pub(crate) duration_sec: Option<f32>,
    #[arg(short = 'g', long)]
    pub(crate) mic_gain: Option<f32>,
    #[arg(long, default_value_t = 12_000.0)]
    pub(crate) in_hp_hz: f32,
    #[arg(long, default_value_t = 19_000.0)]
    pub(crate) in_lp_hz: f32,
    #[arg(long)]
    pub(crate) no_input_filter: bool,
    #[arg(short = 'W', long)]
    pub(crate) dump_wav: Option<String>,
    #[arg(short = 's', long)]
    pub(crate) spectrogram: bool,
    #[arg(short = 'P', long, default_value = "/tmp/rx_spectrogram.png")]
    pub(crate) spectrogram_path: String,
    #[arg(long, default_value_t = 512)]
    pub(crate) spectrogram_nfft: usize,
    #[arg(long, default_value_t = 128)]
    pub(crate) spectrogram_hop: usize,
    #[arg(long, value_enum, default_value_t = SpectrogramWindowArg::Hann)]
    pub(crate) spectrogram_window: SpectrogramWindowArg,
    #[arg(short = 'o', long)]
    pub(crate) oracle: bool,
    #[arg(short = 'S', long)]
    pub(crate) stdout: bool,
    #[arg(short = 'v', long)]
    pub(crate) verbose: bool,
    #[arg(short = 'l', long, value_enum)]
    pub(crate) log_level: Option<LogLevelArg>,
    #[arg(short = 'L', long)]
    pub(crate) log_file: Option<String>,
}

pub(crate) fn payload_from_arg(s: &str) -> Vec<u8> {
    s.as_bytes().to_vec()
}

pub(crate) fn payload_from_arg_or_stdin(s: &str) -> Result<Vec<u8>, Box<dyn Error>> {
    if s == "-" {
        let mut buf = Vec::new();
        std::io::stdin().lock().read_to_end(&mut buf)?;
        Ok(buf)
    } else {
        Ok(payload_from_arg(s))
    }
}

#[derive(Clone, Debug)]
pub(crate) struct AudioOpts {
    pub(crate) duration_sec: f32,
    pub(crate) mic_gain: f32,
    pub(crate) spk_gain: f32,
    pub(crate) pre_delay_sec: f32,
    pub(crate) repeats: usize,
    pub(crate) gap_sec: f32,
    pub(crate) input_filter: bool,
    pub(crate) input_hp_hz: f32,
    pub(crate) input_lp_hz: f32,
    pub(crate) dump_wav: Option<String>,
    pub(crate) dump_tx_wav: Option<String>,
    pub(crate) spectrogram: bool,
    pub(crate) spectrogram_path: String,
    pub(crate) spectrogram_opts: SpectrogramOptions,
    pub(crate) oracle: bool,
    pub(crate) verbose: bool,
}

pub(crate) const ORACLE_PAYLOAD: &[u8] = b"ACOUSTIC-OFDM-ORACLE";

pub(crate) fn apply_common_cfg(
    cfg: &mut OfdmConfig,
    common: &CommonCfgArgs,
) -> Result<(), Box<dyn Error>> {
    if let Some(hz) = common.base_freq_hz {
        if !hz.is_finite() || hz <= 0.0 {
            return Err("base frequency must be a positive finite number".into());
        }
        cfg.base_freq_hz = Some(hz);
    }
    if let Some(fs_baseband) = common.fs_baseband {
        if !fs_baseband.is_finite() || fs_baseband <= 0.0 {
            return Err("baseband sample rate must be a positive finite number".into());
        }
        cfg.fs_baseband = fs_baseband;
    }
    if let Some(nfft) = common.nfft {
        if nfft == 0 {
            return Err("nfft must be positive".into());
        }
        cfg.nfft = nfft;
    }
    if let Some(ncp) = common.ncp {
        cfg.ncp = ncp;
    }
    if let Some(sync_half_len) = common.sync_half_len {
        if sync_half_len == 0 {
            return Err("sync-half-len must be positive".into());
        }
        cfg.sync_half_len = sync_half_len;
    }
    if let Some(m) = common.modulation {
        cfg.modulation = m.into();
    }
    if let Some(m) = common.equalizer_mode {
        cfg.equalizer = m.into();
    }
    if let Some(m) = common.fec_mode {
        cfg.fec_mode = m.into();
    }
    if let Some(m) = common.passband_mode {
        cfg.passband_mode = m.into();
    }
    if let Some(w) = common.wake_preamble {
        cfg.wake_preamble = w.into();
    }
    Ok(())
}

pub(crate) fn resolved_log_level(explicit: Option<LogLevelArg>, verbose: bool) -> log::LevelFilter {
    explicit.map(Into::into).unwrap_or(if verbose {
        log::LevelFilter::Debug
    } else {
        log::LevelFilter::Info
    })
}

pub(crate) fn apply_profile_cfg(
    cfg: &mut OfdmConfig,
    profile: LiveProfileArg,
    apply_default_wake: bool,
) {
    match profile {
        LiveProfileArg::Legacy4481483 => {
            cfg.fc = 17_000.0;
            cfg.used_bins = vec![2, 3, 4, 5];
            cfg.pilot_bins = vec![2, 4];
            cfg.retrain_interval_data_symbols = Some(1);
            cfg.terminal_training_symbol = true;
            cfg.wake_freq = 16_500.0;
            cfg.wake_guard_ms = 15.0;
        }
        LiveProfileArg::Standard | LiveProfileArg::LiveDebug => {}
    }
    if apply_default_wake {
        cfg.wake_preamble = match profile {
            LiveProfileArg::Standard => cfg.wake_preamble,
            LiveProfileArg::LiveDebug | LiveProfileArg::Legacy4481483 => WakePreamble::Tone,
        };
    }
}

pub(crate) fn apply_tx_profile_cfg(cfg: &mut OfdmConfig, cmd: &TxCmd) {
    apply_profile_cfg(cfg, cmd.profile, cmd.common.wake_preamble.is_none());
}

pub(crate) fn apply_rx_profile_cfg(cfg: &mut OfdmConfig, cmd: &RxCmd) {
    apply_profile_cfg(cfg, cmd.profile, cmd.common.wake_preamble.is_none());
}

pub(crate) fn apply_mic_roundtrip_profile_cfg(cfg: &mut OfdmConfig, cmd: &MicRoundtripCmd) {
    apply_profile_cfg(cfg, cmd.profile, cmd.common.wake_preamble.is_none());
}

pub(crate) fn rx_audio_opts(cmd: &RxCmd) -> AudioOpts {
    let (duration_sec, mic_gain, dump_wav, spectrogram, spectrogram_path, oracle, verbose) =
        match cmd.profile {
            LiveProfileArg::Standard => (
                10.0,
                1.0,
                None,
                false,
                "/tmp/rx_spectrogram.png".to_string(),
                false,
                false,
            ),
            LiveProfileArg::LiveDebug => (
                5.0,
                0.2,
                Some("/tmp/rx_capture.wav".to_string()),
                true,
                "/tmp/rx_spectrogram.png".to_string(),
                true,
                true,
            ),
            LiveProfileArg::Legacy4481483 => (
                5.0,
                0.2,
                Some("/tmp/rx_capture.wav".to_string()),
                true,
                "/tmp/rx_spectrogram.png".to_string(),
                true,
                true,
            ),
        };
    AudioOpts {
        duration_sec: cmd.duration_sec.unwrap_or(duration_sec),
        mic_gain: cmd.mic_gain.unwrap_or(mic_gain),
        spk_gain: 1.0,
        pre_delay_sec: 0.2,
        repeats: 3,
        gap_sec: 0.35,
        input_filter: false,
        input_hp_hz: cmd.in_hp_hz,
        input_lp_hz: cmd.in_lp_hz,
        dump_wav: cmd.dump_wav.clone().or(dump_wav),
        dump_tx_wav: None,
        spectrogram: cmd.spectrogram || spectrogram,
        spectrogram_path: if cmd.spectrogram_path != "/tmp/rx_spectrogram.png" {
            cmd.spectrogram_path.clone()
        } else {
            spectrogram_path
        },
        spectrogram_opts: SpectrogramOptions {
            nfft: cmd.spectrogram_nfft,
            hop: cmd.spectrogram_hop,
            window: cmd.spectrogram_window.into(),
        },
        oracle: cmd.oracle || oracle,
        verbose: cmd.verbose
            || verbose
            || matches!(cmd.log_level, Some(LogLevelArg::Debug | LogLevelArg::Trace)),
    }
}

pub(crate) fn tx_audio_opts(cmd: &TxCmd) -> AudioOpts {
    let (spk_gain, pre_delay_sec, repeats, gap_sec, oracle, verbose) = match cmd.profile {
        LiveProfileArg::Standard => (1.0, 0.2, 3, 0.35, false, false),
        LiveProfileArg::LiveDebug => (0.2, 0.5, 5, 0.35, true, false),
        LiveProfileArg::Legacy4481483 => (0.2, 0.5, 5, 0.35, true, false),
    };
    AudioOpts {
        duration_sec: 10.0,
        mic_gain: 1.0,
        spk_gain: cmd.spk_gain.unwrap_or(spk_gain),
        pre_delay_sec: cmd.pre_delay_sec.unwrap_or(pre_delay_sec),
        repeats: cmd.repeats.unwrap_or(repeats),
        gap_sec: cmd.gap_sec.unwrap_or(gap_sec),
        input_filter: false,
        input_hp_hz: 12_000.0,
        input_lp_hz: 19_000.0,
        dump_wav: cmd.dump_wav.clone(),
        dump_tx_wav: None,
        spectrogram: false,
        spectrogram_path: "/tmp/rx_spectrogram.png".to_string(),
        spectrogram_opts: SpectrogramOptions::default(),
        oracle: cmd.oracle || oracle,
        verbose: cmd.verbose
            || verbose
            || matches!(cmd.log_level, Some(LogLevelArg::Debug | LogLevelArg::Trace)),
    }
}

pub(crate) fn mic_roundtrip_audio_opts(cmd: &MicRoundtripCmd) -> AudioOpts {
    let (
        duration_sec,
        mic_gain,
        spk_gain,
        pre_delay_sec,
        repeats,
        gap_sec,
        dump_wav,
        spectrogram,
        spectrogram_path,
        oracle,
        verbose,
    ) = match cmd.profile {
        LiveProfileArg::Standard => (
            6.0,
            1.0,
            1.0,
            0.2,
            3,
            0.35,
            None,
            false,
            "/tmp/rx_spectrogram.png".to_string(),
            false,
            false,
        ),
        LiveProfileArg::LiveDebug => (
            10.0,
            0.2,
            0.2,
            0.5,
            5,
            0.35,
            Some("/tmp/rx_capture.wav".to_string()),
            true,
            "/tmp/rx_spectrogram.png".to_string(),
            true,
            true,
        ),
        LiveProfileArg::Legacy4481483 => (
            10.0,
            0.2,
            0.2,
            0.5,
            5,
            0.35,
            Some("/tmp/rx_capture.wav".to_string()),
            true,
            "/tmp/rx_spectrogram.png".to_string(),
            true,
            true,
        ),
    };
    AudioOpts {
        duration_sec: cmd.duration_sec.unwrap_or(duration_sec),
        mic_gain: cmd.mic_gain.unwrap_or(mic_gain),
        spk_gain: cmd.spk_gain.unwrap_or(spk_gain),
        pre_delay_sec: cmd.pre_delay_sec.unwrap_or(pre_delay_sec),
        repeats: cmd.repeats.unwrap_or(repeats),
        gap_sec: cmd.gap_sec.unwrap_or(gap_sec),
        input_filter: false,
        input_hp_hz: 12_000.0,
        input_lp_hz: 19_000.0,
        dump_wav: cmd.dump_wav.clone().or(dump_wav),
        dump_tx_wav: cmd.dump_tx_wav.clone(),
        spectrogram: cmd.spectrogram || spectrogram,
        spectrogram_path: if cmd.spectrogram_path != "/tmp/rx_spectrogram.png" {
            cmd.spectrogram_path.clone()
        } else {
            spectrogram_path
        },
        spectrogram_opts: SpectrogramOptions {
            nfft: cmd.spectrogram_nfft,
            hop: cmd.spectrogram_hop,
            window: cmd.spectrogram_window.into(),
        },
        oracle: cmd.oracle || oracle,
        verbose: cmd.verbose
            || verbose
            || matches!(cmd.log_level, Some(LogLevelArg::Debug | LogLevelArg::Trace)),
    }
}
