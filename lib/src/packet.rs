// Copyright (c) 2026 Elias S. G. Carotti

use crate::config::{FecMode, Modulation, OfdmConfig};
use crate::crc::crc16_ccitt;

#[derive(Clone, Debug)]
pub struct PacketInfo {
    /// Packet format version from the serialized header.
    pub version: u8,
    /// On-wire modulation identifier from the packet header.
    pub mod_id: u8,
    /// Session identifier used to group fragments from the same transmission.
    pub session_id: u16,
    /// Zero-based fragment index within the session payload.
    pub frag_index: u8,
    /// Total number of fragments in the session payload.
    pub frag_count: u8,
    /// Application payload bytes carried by this fragment.
    pub payload: Vec<u8>,
}

/// Splits payload bytes into fixed-size chunks.
///
/// Parameters:
/// - `payload`: full payload bytes.
/// - `chunk_size`: max bytes per chunk.
/// Returns:
/// - `Vec<Vec<u8>>`: ordered chunks.
pub fn split_payload(payload: &[u8], chunk_size: usize) -> Vec<Vec<u8>> {
    payload
        .chunks(chunk_size)
        .map(|c| c.to_vec())
        .collect::<Vec<_>>()
}

/// Builds one packet (header + payload + CRC).
///
/// Parameters:
/// - `payload`: fragment payload bytes.
/// - `frag_index`: zero-based fragment index.
/// - `frag_count`: total fragment count.
/// - `cfg`: modem configuration (session/modulation fields used).
/// Returns:
/// - `Vec<u8>`: serialized packet bytes.
pub fn build_packet_bytes(
    payload: &[u8],
    frag_index: u8,
    frag_count: u8,
    cfg: &OfdmConfig,
) -> Vec<u8> {
    let mut body = Vec::with_capacity(11 + payload.len());
    body.extend_from_slice(&[0xA5, 0x5A]);
    body.push(1); // version
    body.push(cfg.modulation.mod_id());
    body.extend_from_slice(&cfg.session_id.to_be_bytes());
    body.push(frag_index);
    body.push(frag_count);
    body.push(payload.len() as u8);
    body.extend_from_slice(payload);
    let crc = crc16_ccitt(&body);
    body.extend_from_slice(&crc.to_be_bytes());
    body
}

/// Parses and validates one packet from a byte stream prefix.
///
/// Parameters:
/// - `rx`: received byte slice.
/// Returns:
/// - `Option<(PacketInfo, usize)>`: parsed packet and consumed byte count, or `None` if invalid.
pub fn parse_packet_bytes(rx: &[u8]) -> Option<(PacketInfo, usize)> {
    if rx.len() < 11 || rx[0] != 0xA5 || rx[1] != 0x5A {
        return None;
    }
    let plen = rx[8] as usize;
    let total = 9 + plen + 2;
    if rx.len() < total {
        return None;
    }
    let body = &rx[..9 + plen];
    let rx_crc = u16::from_be_bytes([rx[9 + plen], rx[10 + plen]]);
    if crc16_ccitt(body) != rx_crc {
        return None;
    }
    Some((
        PacketInfo {
            version: rx[2],
            mod_id: rx[3],
            session_id: u16::from_be_bytes([rx[4], rx[5]]),
            frag_index: rx[6],
            frag_count: rx[7],
            payload: rx[9..9 + plen].to_vec(),
        },
        total,
    ))
}

#[derive(Clone, Debug)]
pub struct PacketParseAttempt {
    /// Whether the leading sync/preamble bytes matched `0xA5 0x5A`.
    pub preamble_ok: bool,
    /// Whether enough bytes were available to read the fixed-size header.
    pub enough_for_header: bool,
    /// Parsed payload length from the header when available.
    pub payload_len: Option<usize>,
    /// Total packet length implied by the header when available.
    pub total_len: Option<usize>,
    /// Whether the received byte stream was long enough for the full packet.
    pub enough_for_total: bool,
    /// Whether the packet CRC matched.
    pub crc_ok: bool,
    /// Fully parsed packet when both header and CRC checks succeeded.
    pub parsed: Option<PacketInfo>,
}

pub fn inspect_packet_bytes(rx: &[u8]) -> PacketParseAttempt {
    if rx.len() < 9 {
        return PacketParseAttempt {
            preamble_ok: rx.len() >= 2 && rx[0] == 0xA5 && rx[1] == 0x5A,
            enough_for_header: false,
            payload_len: None,
            total_len: None,
            enough_for_total: false,
            crc_ok: false,
            parsed: None,
        };
    }

    let preamble_ok = rx[0] == 0xA5 && rx[1] == 0x5A;
    let plen = rx[8] as usize;
    let total = 9 + plen + 2;
    let enough_for_total = rx.len() >= total;
    let crc_ok = if preamble_ok && enough_for_total {
        let body = &rx[..9 + plen];
        let rx_crc = u16::from_be_bytes([rx[9 + plen], rx[10 + plen]]);
        crc16_ccitt(body) == rx_crc
    } else {
        false
    };
    let parsed = if preamble_ok && crc_ok {
        Some(PacketInfo {
            version: rx[2],
            mod_id: rx[3],
            session_id: u16::from_be_bytes([rx[4], rx[5]]),
            frag_index: rx[6],
            frag_count: rx[7],
            payload: rx[9..9 + plen].to_vec(),
        })
    } else {
        None
    };
    PacketParseAttempt {
        preamble_ok,
        enough_for_header: true,
        payload_len: Some(plen),
        total_len: Some(total),
        enough_for_total,
        crc_ok,
        parsed,
    }
}

/// Expands bytes into MSB-first bits.
///
/// Parameters:
/// - `bytes`: input bytes.
/// Returns:
/// - `Vec<u8>`: bit vector with values in `{0,1}`.
pub fn bytes_to_bits(bytes: &[u8]) -> Vec<u8> {
    let mut bits = Vec::with_capacity(bytes.len() * 8);
    for &b in bytes {
        for shift in (0..8).rev() {
            bits.push((b >> shift) & 1);
        }
    }
    bits
}

/// Packs MSB-first bits into bytes (truncates incomplete trailing byte).
///
/// Parameters:
/// - `bits`: input bits (`0/1` values).
/// Returns:
/// - `Vec<u8>`: packed bytes.
pub fn bits_to_bytes(bits: &[u8]) -> Vec<u8> {
    let nbytes = bits.len() / 8;
    let mut out = vec![0u8; nbytes];
    for i in 0..nbytes {
        let mut v = 0u8;
        for b in 0..8 {
            v = (v << 1) | (bits[i * 8 + b] & 1);
        }
        out[i] = v;
    }
    out
}

pub fn fec_encoded_bits_len(raw_bits_len: usize, mode: FecMode) -> usize {
    match mode {
        FecMode::None => raw_bits_len,
        FecMode::Hamming74 => raw_bits_len.div_ceil(4) * 7,
    }
}

pub fn fec_encode_bits(bits: &[u8], mode: FecMode) -> Vec<u8> {
    match mode {
        FecMode::None => bits.to_vec(),
        FecMode::Hamming74 => {
            let mut out = Vec::with_capacity(fec_encoded_bits_len(bits.len(), mode));
            let mut i = 0usize;
            while i < bits.len() {
                let d1 = bits.get(i).copied().unwrap_or(0) & 1;
                let d2 = bits.get(i + 1).copied().unwrap_or(0) & 1;
                let d3 = bits.get(i + 2).copied().unwrap_or(0) & 1;
                let d4 = bits.get(i + 3).copied().unwrap_or(0) & 1;
                let p1 = d1 ^ d2 ^ d4;
                let p2 = d1 ^ d3 ^ d4;
                let p3 = d2 ^ d3 ^ d4;
                out.extend_from_slice(&[p1, p2, d1, p3, d2, d3, d4]);
                i += 4;
            }
            out
        }
    }
}

pub fn fec_decode_bits(bits: &[u8], mode: FecMode) -> Vec<u8> {
    match mode {
        FecMode::None => bits.to_vec(),
        FecMode::Hamming74 => {
            let mut out = Vec::with_capacity((bits.len() / 7) * 4);
            for chunk in bits.chunks_exact(7) {
                let mut c = [
                    chunk[0] & 1,
                    chunk[1] & 1,
                    chunk[2] & 1,
                    chunk[3] & 1,
                    chunk[4] & 1,
                    chunk[5] & 1,
                    chunk[6] & 1,
                ];
                let s1 = c[0] ^ c[2] ^ c[4] ^ c[6];
                let s2 = c[1] ^ c[2] ^ c[5] ^ c[6];
                let s3 = c[3] ^ c[4] ^ c[5] ^ c[6];
                let syndrome = (s1 | (s2 << 1) | (s3 << 2)) as usize;
                if (1..=7).contains(&syndrome) {
                    c[syndrome - 1] ^= 1;
                }
                out.extend_from_slice(&[c[2], c[4], c[5], c[6]]);
            }
            out
        }
    }
}

/// Converts header modulation ID to enum.
///
/// Parameters:
/// - `id`: modulation ID from packet header.
/// Returns:
/// - `Option<Modulation>`: matching modulation, or `None` if unknown.
pub fn modulation_from_id(id: u8) -> Option<Modulation> {
    match id {
        1 => Some(Modulation::Bpsk),
        2 => Some(Modulation::Qpsk),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    /// Checks deterministic payload chunking.
    fn split_payload_chunks_correctly() {
        let p: Vec<u8> = (0..10).collect();
        let chunks = split_payload(&p, 4);
        assert_eq!(chunks.len(), 3);
        assert_eq!(chunks[0], vec![0, 1, 2, 3]);
        assert_eq!(chunks[1], vec![4, 5, 6, 7]);
        assert_eq!(chunks[2], vec![8, 9]);
    }

    #[test]
    /// Ensures bits<->bytes conversions are inverse.
    fn bits_bytes_roundtrip() {
        let data = vec![0x00, 0xA5, 0x5A, 0xFF, 0x13];
        let bits = bytes_to_bits(&data);
        let back = bits_to_bytes(&bits);
        assert_eq!(back, data);
    }

    #[test]
    /// Verifies packet build/parse consistency.
    fn packet_build_parse_roundtrip() {
        let cfg = OfdmConfig::default();
        let payload = vec![1, 2, 3, 4, 5];
        let pkt = build_packet_bytes(&payload, 2, 9, &cfg);
        let (info, used) = parse_packet_bytes(&pkt).expect("packet must parse");
        assert_eq!(used, pkt.len());
        assert_eq!(info.version, 1);
        assert_eq!(info.mod_id, cfg.modulation.mod_id());
        assert_eq!(info.session_id, cfg.session_id);
        assert_eq!(info.frag_index, 2);
        assert_eq!(info.frag_count, 9);
        assert_eq!(info.payload, payload);
    }

    #[test]
    /// Ensures packets with CRC corruption are rejected.
    fn packet_crc_error_is_rejected() {
        let cfg = OfdmConfig::default();
        let payload = vec![10, 20, 30];
        let mut pkt = build_packet_bytes(&payload, 0, 1, &cfg);
        let n = pkt.len();
        pkt[n - 1] ^= 0x01;
        assert!(parse_packet_bytes(&pkt).is_none());
    }

    #[test]
    fn hamming74_roundtrip_and_single_bit_correction() {
        let raw = bytes_to_bits(&[0xA5, 0x5A, 0x13, 0x7C]);
        let mut enc = fec_encode_bits(&raw, FecMode::Hamming74);
        assert_eq!(enc.len(), raw.len() / 4 * 7);
        enc[5] ^= 1;
        enc[20] ^= 1;
        let dec = fec_decode_bits(&enc, FecMode::Hamming74);
        assert_eq!(dec, raw);
    }
}

// vim: set ts=4 sw=4 et:
