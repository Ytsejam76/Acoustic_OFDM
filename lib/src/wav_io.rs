// Copyright (c) 2026 Elias S. G. Carotti

use std::error::Error;
use std::path::Path;

/// Saves mono PCM samples to a 16-bit WAV file.
///
/// Parameters:
/// - `path`: destination WAV path.
/// - `samples`: normalized mono samples (typically in `[-1, 1]`).
/// - `sample_rate`: sample rate in Hz.
/// Returns:
/// - `Result<(), Box<dyn Error>>`: `Ok(())` on success.
pub fn save_wav_mono_i16(path: &Path, samples: &[f32], sample_rate: u32) -> Result<(), Box<dyn Error>> {
    let spec = hound::WavSpec {
        channels: 1,
        sample_rate,
        bits_per_sample: 16,
        sample_format: hound::SampleFormat::Int,
    };
    let mut wr = hound::WavWriter::create(path, spec)?;
    for &s in samples {
        let v = (s.clamp(-1.0, 1.0) * (i16::MAX as f32)).round() as i16;
        wr.write_sample(v)?;
    }
    wr.finalize()?;
    Ok(())
}

/// Loads a mono WAV file as normalized `f32` samples.
///
/// Parameters:
/// - `path`: source WAV path.
/// Returns:
/// - `Result<(Vec<f32>, u32), Box<dyn Error>>`: `(samples, sample_rate)` on success.
pub fn load_wav_mono_f32(path: &Path) -> Result<(Vec<f32>, u32), Box<dyn Error>> {
    let mut rd = hound::WavReader::open(path)?;
    let spec = rd.spec();
    if spec.channels != 1 {
        return Err("only mono WAV is supported".into());
    }
    let sr = spec.sample_rate;
    let samples = match (spec.sample_format, spec.bits_per_sample) {
        (hound::SampleFormat::Int, 16) => rd
            .samples::<i16>()
            .map(|r| r.map(|v| v as f32 / i16::MAX as f32))
            .collect::<Result<Vec<_>, _>>()?,
        _ => return Err("unsupported WAV format (expected 16-bit PCM mono)".into()),
    };
    Ok((samples, sr))
}

#[cfg(test)]
mod tests {
    use super::{load_wav_mono_f32, save_wav_mono_i16};
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    /// Builds a temporary path for WAV tests.
    ///
    /// Parameters:
    /// - none.
    /// Returns:
    /// - `PathBuf`: temporary WAV file path.
    fn tmp_wav_path() -> PathBuf {
        let ts = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock drift")
            .as_nanos();
        std::env::temp_dir().join(format!("acoustic_ofdm_wav_test_{ts}.wav"))
    }

    #[test]
    /// Verifies WAV write/read round-trip.
    fn wav_roundtrip_mono_i16() {
        let path = tmp_wav_path();
        let src = vec![0.0f32, 0.25, -0.25, 0.9, -0.9, 0.1, -0.1];
        save_wav_mono_i16(&path, &src, 48_000).expect("write wav");
        let (dst, sr) = load_wav_mono_f32(&path).expect("read wav");
        assert_eq!(sr, 48_000);
        assert_eq!(dst.len(), src.len());
        for (a, b) in src.iter().zip(dst.iter()) {
            assert!((a - b).abs() < 2.0 / i16::MAX as f32);
        }
        let _ = std::fs::remove_file(path);
    }
}

// vim: set ts=4 sw=4 et:
