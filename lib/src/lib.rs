// Copyright (c) 2026 Elias S. G. Carotti

pub mod config;
pub mod constellation;
pub mod crc;
pub mod modem;
pub mod packet;
pub mod spectrogram;
pub mod wav_io;

pub use config::{Modulation, OfdmConfig, WakePreamble};
pub use constellation::save_constellation_comparison_png;
pub use modem::{
    decode_encoded_burst_oracle, decode_packet_baseband, decode_single_packet_passband,
    diagnose_passband_window, dump_passband_bins, dump_passband_constellation,
    dump_passband_pilot_tracking, dump_passband_sync_metric, encode_payload,
    encode_single_packet_passband, EncodedBurst, PassbandBinDump, PassbandBinDumpRow,
    PassbandConstellationDump, PassbandDiagnostics, PassbandPilotTrackDump, PassbandSyncDump,
};
pub use spectrogram::{
    save_spectrogram_png, save_spectrogram_png_with_options, SpectrogramOptions, SpectrogramWindow,
};
pub use wav_io::{load_wav_mono_f32, save_wav_mono_i16};

// vim: set ts=4 sw=4 et:
