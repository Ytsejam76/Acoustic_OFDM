// Copyright (c) 2026 Elias S. G. Carotti

pub mod config;
pub mod constellation;
pub mod crc;
pub mod packet;
pub mod modem;
pub mod spectrogram;
pub mod wav_io;

pub use config::{Modulation, OfdmConfig, WakePreamble};
pub use modem::{
    diagnose_passband_window,
    dump_passband_constellation,
    dump_passband_pilot_tracking,
    dump_passband_sync_metric,
    decode_encoded_burst_oracle,
    decode_packet_baseband,
    decode_single_packet_passband,
    encode_payload,
    encode_single_packet_passband,
    EncodedBurst,
    PassbandConstellationDump,
    PassbandDiagnostics,
    PassbandPilotTrackDump,
    PassbandSyncDump,
};
pub use constellation::save_constellation_comparison_png;
pub use spectrogram::save_spectrogram_png;
pub use wav_io::{load_wav_mono_f32, save_wav_mono_i16};

// vim: set ts=4 sw=4 et:
