// Copyright (c) 2026 Elias S. G. Carotti

mod baseband;
pub mod config;
pub mod constellation;
pub mod crc;
pub mod modem;
pub mod packet;
pub mod spectrogram;
pub mod wav_io;

pub use baseband::{
    decode_packet_baseband, encode_single_packet_baseband, expected_single_packet_data_symbols,
    recover_single_packet_data_symbols,
};
pub use config::{EqualizerMode, FecMode, Modulation, OfdmConfig, PassbandMode, WakePreamble};
pub use constellation::{save_channel_compare_png, save_constellation_comparison_png};
pub use modem::{
    decode_encoded_burst_oracle, decode_single_packet_passband,
    decode_single_packet_passband_with_sync, diagnose_passband_window,
    diagnose_passband_window_with_sync, dump_passband_bins, dump_passband_bins_with_sync,
    dump_passband_channel_compare_with_sync, dump_passband_constellation, dump_passband_iq_chain,
    dump_passband_pilot_tracking, dump_passband_sync_metric, encode_payload,
    encode_single_packet_passband, encode_single_packet_passband_body,
    recover_decided_packet_bytes_passband_with_sync, EncodedBurst, PassbandBinDump,
    PassbandBinDumpRow, PassbandChannelCompareDump, PassbandChannelCompareRow,
    PassbandConstellationDump, PassbandDiagnostics, PassbandIqChainDump, PassbandPilotTrackDump,
    PassbandSyncDump,
};
pub use packet::{inspect_packet_bytes, PacketParseAttempt};
pub use rustfft::num_complex::Complex32;
pub use spectrogram::{
    save_spectrogram_png, save_spectrogram_png_with_options, SpectrogramOptions, SpectrogramWindow,
};
pub use wav_io::{load_wav_mono_f32, save_wav_mono_i16};

// vim: set ts=4 sw=4 et:
