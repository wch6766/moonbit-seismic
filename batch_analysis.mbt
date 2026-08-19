///|
/// Per-window feature record for batch triage.
pub struct WindowFeatureRow {
  start_time : Double
  end_time : Double
  rms : Double
  peak : Double
  dominant_frequency : Double
  spectral_entropy : Double
  quality_score : Double
  event_count : Int
}

///|
/// Summary of a batch analysis.
pub struct BatchAnalysis {
  rows : Array[WindowFeatureRow]
  total_windows : Int
  accepted_windows : Int
  mean_quality : Double
  strongest_peak : Double
}

///|
/// Analyze each trace window using the same deterministic feature pipeline.
pub fn analyze_trace_windows(
  trace : Trace,
  window_size : Int,
  step : Int,
  minimum_quality : Double,
) -> BatchAnalysis {
  let rows = Array::new()
  let windows = trace_windows(trace, window_size, step, false)
  let mut quality_sum = 0.0
  let mut accepted = 0
  let mut strongest = 0.0
  for window in windows {
    let values = window.data
    let report = assess_quality(
      Trace::new(
        trace.station,
        trace.channel,
        trace.sampling_rate,
        window.start_time,
        values,
      ),
      1.0e6,
      0.0000001,
    )
    let spectrum = hann_periodogram(window.data, trace.sampling_rate)
    let features = extract_spectral_features(spectrum)
    let detector = default_event_config(trace.delta(), window.start_time)
    let events = detect_events(window.data, detector)
    let peak = maximum_sample(window.data)
      .abs()
      .max(minimum_sample(window.data).abs())
    rows.push({
      start_time: window.start_time,
      end_time: window.end_time,
      rms: rms_samples(window.data),
      peak,
      dominant_frequency: features.dominant_frequency,
      spectral_entropy: features.entropy,
      quality_score: report.score,
      event_count: events.length(),
    })
    quality_sum += report.score
    if report.score >= minimum_quality {
      accepted += 1
    }
    if peak > strongest {
      strongest = peak
    }
  }
  let total = rows.length()
  {
    rows,
    total_windows: total,
    accepted_windows: accepted,
    mean_quality: if total == 0 {
      0.0
    } else {
      quality_sum / total.to_double()
    },
    strongest_peak: strongest,
  }
}

///|
/// Return rows containing at least one detected event.
pub fn event_windows(batch : BatchAnalysis) -> Array[WindowFeatureRow] {
  batch.rows.filter(row => row.event_count > 0)
}

///|
/// Return rows whose RMS exceeds a physical amplitude threshold.
pub fn energetic_windows(
  batch : BatchAnalysis,
  threshold : Double,
) -> Array[WindowFeatureRow] {
  batch.rows.filter(row => row.rms >= threshold)
}

///|
/// Return a quality-weighted peak score for ranking windows.
pub fn window_priority(row : WindowFeatureRow) -> Double {
  row.peak * row.quality_score * (1.0 + row.event_count.to_double())
}

///|
/// Return the highest-priority row, if any.
pub fn strongest_window(batch : BatchAnalysis) -> WindowFeatureRow? {
  if batch.rows.is_empty() {
    return None
  }
  let mut best = batch.rows[0]
  for index = 1; index < batch.rows.length(); index = index + 1 {
    if window_priority(batch.rows[index]) > window_priority(best) {
      best = batch.rows[index]
    }
  }
  Some(best)
}

///|
/// Aggregate several traces into one stack and analyze the stack.
pub fn analyze_trace_collection(
  traces : Array[Trace],
  window_size : Int,
  step : Int,
) -> Result[BatchAnalysis, String] {
  match stack_traces(traces) {
    Ok(stacked) => Ok(analyze_trace_windows(stacked, window_size, step, 0.5))
    Err(message) => Err(message)
  }
}

///|
/// Return a stable textual table for human inspection.
pub fn batch_analysis_text(batch : BatchAnalysis) -> String {
  let mut output = "start,end,rms,peak,dominant_hz,entropy,quality,events\n"
  for row in batch.rows {
    output = output +
      row.start_time.to_string() +
      "," +
      row.end_time.to_string() +
      "," +
      row.rms.to_string() +
      "," +
      row.peak.to_string() +
      "," +
      row.dominant_frequency.to_string() +
      "," +
      row.spectral_entropy.to_string() +
      "," +
      row.quality_score.to_string() +
      "," +
      row.event_count.to_string() +
      "\n"
  }
  output
}
