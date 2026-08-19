///|
/// Configuration for a reproducible analysis pipeline.
pub struct ProcessingConfig {
  remove_mean : Bool
  remove_trend : Bool
  taper_fraction : Double
  lowpass_frequency : Double?
  highpass_frequency : Double?
  fir_cutoff : Double?
  fir_taps : Int
  resample_rate : Double?
  zero_phase : Bool
  quality_clip_limit : Double
  quality_flat_tolerance : Double
}

///|
/// Output of the deterministic processing pipeline.
pub struct ProcessingResult {
  original : TraceSummary
  processed : Trace
  quality : QualityReport
  dominant_frequency : Double
  events : Array[EventCandidate]
}

///|
/// Return a conservative pipeline configuration.
pub fn default_processing_config() -> ProcessingConfig {
  {
    remove_mean: true,
    remove_trend: true,
    taper_fraction: 0.05,
    lowpass_frequency: None,
    highpass_frequency: None,
    fir_cutoff: None,
    fir_taps: 31,
    resample_rate: None,
    zero_phase: true,
    quality_clip_limit: 1.0e6,
    quality_flat_tolerance: 0.0000001,
  }
}

///|
/// Apply a configured sequence of preprocessing operations.
pub fn process_trace(
  trace : Trace,
  config : ProcessingConfig,
) -> Result[Trace, String] {
  let validation = validate_trace(trace)
  match validation {
    Err(message) => Err(message)
    Ok(_) => {
      let mut result = preprocess_trace(
        trace,
        config.remove_mean,
        config.remove_trend,
        config.taper_fraction,
      )
      match config.highpass_frequency {
        Some(frequency) => highpass(result, frequency, config.zero_phase)
        None => ()
      }
      match config.lowpass_frequency {
        Some(frequency) => lowpass(result, frequency, config.zero_phase)
        None => ()
      }
      match config.fir_cutoff {
        Some(cutoff) =>
          match
            design_fir_lowpass(
              result.sampling_rate,
              cutoff,
              config.fir_taps,
              fir_hann(),
            ) {
            Ok(filter) =>
              result = filter_trace_fir(result, filter, config.zero_phase)
            Err(message) => return Err(message)
          }
        None => ()
      }
      match config.resample_rate {
        Some(rate) => result = resample_trace(result, rate)
        None => ()
      }
      Ok(result)
    }
  }
}

///|
/// Run preprocessing, QC, spectral analysis, and event detection together.
pub fn analyze_trace(
  trace : Trace,
  config : ProcessingConfig,
) -> Result[ProcessingResult, String] {
  let original = summarize_trace(trace)
  let processed = match process_trace(trace, config) {
    Ok(value) => value
    Err(message) => return Err(message)
  }
  let quality = assess_quality(
    processed,
    config.quality_clip_limit,
    config.quality_flat_tolerance,
  )
  let spectrum = hann_periodogram(processed.data, processed.sampling_rate)
  let detector = default_event_config(processed.delta(), processed.start_time)
  Ok({
    original,
    processed,
    quality,
    dominant_frequency: dominant_frequency(spectrum),
    events: detect_events(processed.data, detector),
  })
}

///|
/// Run a low-pass processing shortcut with documented defaults.
pub fn clean_seismic_trace(
  trace : Trace,
  cutoff : Double,
) -> Result[Trace, String] {
  let config = {
    remove_mean: true,
    remove_trend: true,
    taper_fraction: 0.05,
    lowpass_frequency: Some(cutoff),
    highpass_frequency: None,
    fir_cutoff: None,
    fir_taps: 31,
    resample_rate: None,
    zero_phase: true,
    quality_clip_limit: 1.0e6,
    quality_flat_tolerance: 0.0000001,
  }
  process_trace(trace, config)
}

///|
/// Return a quality score after processing a trace.
pub fn processed_quality_score(
  trace : Trace,
  config : ProcessingConfig,
) -> Result[Double, String] {
  match analyze_trace(trace, config) {
    Ok(result) => Ok(result.quality.score)
    Err(message) => Err(message)
  }
}

///|
/// Return the strongest event from an analysis result.
pub fn strongest_processed_event(result : ProcessingResult) -> EventCandidate? {
  strongest_event(result.events)
}

///|
/// Produce a compact pipeline status line for the CLI.
pub fn processing_status(result : ProcessingResult) -> String {
  "samples=" +
  result.processed.npts().to_string() +
  " quality=" +
  result.quality.score.to_string() +
  " dominant_hz=" +
  result.dominant_frequency.to_string() +
  " events=" +
  result.events.length().to_string()
}
