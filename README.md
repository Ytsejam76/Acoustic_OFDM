# Acoustic OFDM toy modem 

This repository is a reconstruction of an older (2014) personal experiment. 

The original code is buried deep into some dusty backup disk, but with the help of ChatGPT it is easy to tinker again with it from scratch.

The original goal was not Octave itself. The real target was to build an **audio authentication/data modem for phones**, with the production implementation intended for **Java on Android**. There was also an **Objective-C port for iOS**. Octave/Matlab was used as the exploratory environment: a place to study modulation ideas, visualize constellations, test synchronization strategies, and iterate quickly before committing to the actual mobile code.

The idea was to send a **short authentication-sized payload** over the phone audio path, using the roughly **15–19 kHz** region instead of ordinary low-frequency audio. An OFDM design with a passband carrier was explored first. In practice, the difficult part was synchronization: detecting the packet reliably, aligning the FFT window, and keeping the decoder stable in the presence of real phone audio imperfections.

Because of that, the practical implementation for the original project eventually moved to **8-FSK**, which was easier to synchronize and detect robustly with **Goertzel-based logic**. Also, since I needed it to work on the most diverse cellphones (both cheap and expensive) using higher frequency wasn't robust enough, and I reverted to audible band tones. The 8-FSK, in that respect, helped since it was less unpleasant to human ears than OFDM. 
That simpler design was the pragmatic choice at the time.

Since I never really came to terms with my failure to have a production-ready working OFDM modem, this is my new stab at it, just for hacking, study, tinkering, and pure NERD fun.
It may not work at all since I am just tinkering.

## Historical note

My original plan, in 2014 was quite ambitious, as I chose  a **16-QAM-style constellation**. 
<img src=images/16-QAM.jpg/>
So that older work  aimed at a much more complex constellation the minimal BPSK/QPSK bring-up path used here. This reconstruction deliberately starts from the simpler side again so that synchronization and packet structure can be made solid first.

## Design goals of this reconstruction

This version focuses on a short-burst packet modem rather than continuous transmission:

- short **packetized** bursts rather than a streaming link
- **wake tone** for coarse detection
- **repeated-half sync symbol** for OFDM timing and coarse CFO estimation
- **training OFDM symbol** for channel estimation
- **BPSK first** for bring-up and debugging
- **QPSK support** already present for the next step
- passband transmission around **17 kHz**

The intent is just some nerdy fun.

From Octave:

```matlab
[result, dbg, tx_audio, rx_audio, payload, p] = acoustic_ofdm_test_channel();
```

To try QPSK:

```matlab
p = struct();
p.modulation = 'QPSK';
[result, dbg, tx_audio, rx_audio, payload, p] = acoustic_ofdm_test_channel(p);
```

## What the test channel does

The simulated channel is intentionally simple but useful. It can include:

- leading timing offset
- mild passband multipath
- small carrier frequency offset
- optional amplitude ripple
- AWGN

That is enough to exercise the decoder without jumping straight to speaker/microphone tests.

## Notes

This is deliberately a **small packet modem**. That matches the original use case better than a long streaming design: short authentication-sized payloads, synchronized once per burst, with the emphasis on experimental understanding rather than production throughput.
