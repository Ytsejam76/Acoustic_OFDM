use crate::config::{Modulation, OfdmConfig};
use crate::crc::crc16_ccitt;

#[derive(Clone, Debug)]
pub struct PacketInfo {
    pub version: u8,
    pub mod_id: u8,
    pub session_id: u16,
    pub frag_index: u8,
    pub frag_count: u8,
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
pub fn build_packet_bytes(payload: &[u8], frag_index: u8, frag_count: u8, cfg: &OfdmConfig) -> Vec<u8> {
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
}

// vim: set ts=4 sw=4 et:
