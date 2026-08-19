# moonbit-seismic

MoonBit library for reproducible seismic waveform processing. It provides a
small, dependency-light toolkit for SAC and miniSEED data, signal conditioning,
phase picking, event quality control, magnitude estimation, and one-dimensional
travel-time analysis.

## Positioning

The package is designed for three practical workflows:

- reading and validating waveform records before analysis;
- composing deterministic preprocessing and event-detection pipelines; and
- building teaching, research, and command-line tools without a large runtime.

The public API is organized around `Trace`, immutable configuration records,
and explicit `Result` values for malformed or physically invalid inputs.

## Core capabilities

- SAC binary parsing and serialization, including endian-aware headers;
- miniSEED fixed-section parsing and uncompressed integer/float samples;
- detrending, tapering, resampling, moving-average, IIR, FIR, and zero-phase
  filtering;
- normalized cross-correlation, spectral features, wavelet coefficients, and
  waveform comparison;
- classic STA/LTA and configurable event detection with quality reports;
- P/S phase picking, magnitude estimates, station geometry, and network event
  association helpers;
- layered velocity models, IASPEI91 travel-time calculations, interpolation,
  calibration, and instrument-response utilities;
- streaming accumulators, batch analysis, text waveform import, validation,
  and reproducible synthetic signals.

## Quick start

Install the current stable MoonBit toolchain, then run the project checks:

```bash
moon check --target all --deny-warn
moon test --target all --deny-warn
moon run cmd/main
```

The package has no third-party runtime dependency beyond
`moonbitlang/core/math` and `moonbitlang/core/strconv`.

## CLI

The demonstration program exercises synthetic waveform creation, filtering,
STA/LTA picking, magnitude estimation, travel times, interpolation, and SAC
serialization:

```bash
moon run cmd/main
```

The deterministic acceptance workload is available as a separate benchmark
command:

```bash
moon run cmd/benchmark
```

## Library example

```moonbit nocheck
let event = @seismic.SyntheticEvent::new(120, 1.0, 5.0, 80, 0.0)
let trace = @seismic.synthetic_trace(
  "DEMO",
  "BHZ",
  1024,
  100.0,
  event,
  0.02,
)
let result = @seismic.analyze_trace(
  trace,
  @seismic.default_processing_config(),
)
```

## Architecture

The implementation is intentionally split by domain so that focused APIs can
be reused independently:

```text
Trace and validation
  types, trace_ops, trace_statistics, data_validation
Signal processing
  signal_stats, signal_windows, waveform_math, filter_diagnostics
  fir_filters, adaptive_filters, spectral, spectral_features, wavelet_analysis
  streaming, resampling, preprocessing
Formats and metadata
  sac, sac_metadata, miniseed, miniseed_tools, text_formats
Analysis and workflow
  correlation, correlation_tools, phase_picking, phase_analysis
  event_detection, event_catalog, network_analysis, quality_control
  pipeline, batch_analysis, workflow_report
Earth models and calibration
  traveltime, travel_advanced, interpolation, trace_geometry, calibration
  magnitude, magnitude_advanced, instrument_response
```

## Benchmarks

`cmd/benchmark` runs the checked-in `acceptance_workload`: a deterministic
4,096-sample trace at 100 Hz with three synthetic events and 0.02 noise
amplitude. The workload uses fixed seeds and prints a checksum, processed sample
count, detected-event count, and dominant frequency so results can be compared
across toolchain or implementation changes.

The measured command output and machine/toolchain details are recorded in
[BENCHMARKS.md](BENCHMARKS.md). The numbers are local reproducibility data, not
a claim about every CPU or target.

## Testing and quality gates

The repository includes boundary tests for empty and short traces, invalid
sampling rates, malformed text and binary records, constant signals, partial
locations, clipping, duplicate samples, short windows, and degenerate geometry.

Run the complete local suite with:

```bash
moon fmt
moon info
moon check --target all --deny-warn
moon test --target all --deny-warn
```

Native coverage can be generated with `moon test --target native --enable-coverage`
followed by `moon coverage analyze`.

## Continuous integration

GitHub Actions runs on Linux, macOS, and Windows using the stable MoonBit
installer. The workflow checks all targets, denies warnings, runs the complete
test suite, verifies formatting and generated interfaces, collects native
coverage on Linux, and enforces a minimum of 8,000 production `.mbt` source
lines (tests and build artifacts excluded).

## License

Released under the [Apache License 2.0](LICENSE).
