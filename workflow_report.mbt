///|
/// Human-readable acceptance report for one seismic trace.
pub struct WorkflowReport {
  station : String
  channel : String
  source_samples : Int
  processed_samples : Int
  source_duration : Double
  processed_duration : Double
  quality_score : Double
  dominant_frequency : Double
  event_count : Int
  peak_amplitude : Double
  roundtrip_sac : Bool
  roundtrip_miniseed : Bool
}

///|
/// Build a report from the standard analysis pipeline.
pub fn build_workflow_report(
  trace : Trace,
  config : ProcessingConfig,
) -> Result[WorkflowReport, String] {
  let result = match analyze_trace(trace, config) {
    Ok(value) => value
    Err(message) => return Err(message)
  }
  Ok({
    station: trace.station,
    channel: trace.channel,
    source_samples: trace.npts(),
    processed_samples: result.processed.npts(),
    source_duration: trace.duration(),
    processed_duration: result.processed.duration(),
    quality_score: result.quality.score,
    dominant_frequency: result.dominant_frequency,
    event_count: result.events.length(),
    peak_amplitude: peak_amplitude(result.processed),
    roundtrip_sac: validate_sac_roundtrip_size(result.processed),
    roundtrip_miniseed: validate_miniseed_roundtrip_size(result.processed),
  })
}

///|
/// Serialize a report as a stable key-value table for CLI output.
pub fn workflow_report_text(report : WorkflowReport) -> String {
  "station=" +
  report.station +
  "\n" +
  "channel=" +
  report.channel +
  "\n" +
  "source_samples=" +
  report.source_samples.to_string() +
  "\n" +
  "processed_samples=" +
  report.processed_samples.to_string() +
  "\n" +
  "source_duration=" +
  report.source_duration.to_string() +
  "\n" +
  "processed_duration=" +
  report.processed_duration.to_string() +
  "\n" +
  "quality_score=" +
  report.quality_score.to_string() +
  "\n" +
  "dominant_frequency=" +
  report.dominant_frequency.to_string() +
  "\n" +
  "event_count=" +
  report.event_count.to_string() +
  "\n" +
  "peak_amplitude=" +
  report.peak_amplitude.to_string() +
  "\n" +
  "sac_roundtrip=" +
  report.roundtrip_sac.to_string() +
  "\n" +
  "miniseed_roundtrip=" +
  report.roundtrip_miniseed.to_string()
}

///|
/// Return a compact status badge for automation.
pub fn workflow_status(report : WorkflowReport) -> String {
  if report.quality_score >= 0.7 &&
    report.roundtrip_sac &&
    report.roundtrip_miniseed {
    "PASS"
  } else {
    "REVIEW"
  }
}

///|
/// Build and format a report in one call.
pub fn analyze_and_report(trace : Trace) -> Result[String, String] {
  match build_workflow_report(trace, default_processing_config()) {
    Ok(report) =>
      Ok(workflow_status(report) + "\n" + workflow_report_text(report))
    Err(message) => Err(message)
  }
}

///|
/// Return true when both binary format round trips are available.
pub fn report_roundtrips_pass(report : WorkflowReport) -> Bool {
  report.roundtrip_sac && report.roundtrip_miniseed
}

///|
/// Return a one-line report suitable for benchmark logs.
pub fn workflow_report_line(report : WorkflowReport) -> String {
  workflow_status(report) +
  " " +
  report.station +
  "." +
  report.channel +
  " samples=" +
  report.processed_samples.to_string() +
  " quality=" +
  report.quality_score.to_string() +
  " events=" +
  report.event_count.to_string()
}

///|
/// Return a normalized quality score suitable for a dashboard gauge.
pub fn report_quality_percent(report : WorkflowReport) -> Double {
  report.quality_score.max(0.0).min(1.0) * 100.0
}

///|
/// Return the number of source samples removed by processing.
pub fn report_sample_reduction(report : WorkflowReport) -> Int {
  report.source_samples - report.processed_samples
}

///|
/// Return whether the report contains at least one usable event.
pub fn report_has_event(report : WorkflowReport) -> Bool {
  report.event_count > 0 && report.peak_amplitude > 0.0
}

///|
/// Return a machine-readable pass/fail value for CI smoke checks.
pub fn report_passes(report : WorkflowReport) -> Bool {
  workflow_status(report) == "PASS"
}

///|
/// Return the processed duration in milliseconds.
pub fn report_duration_milliseconds(report : WorkflowReport) -> Double {
  report.processed_duration * 1000.0
}

///|
/// Return the report's dominant frequency in millihertz.
pub fn report_dominant_millihertz(report : WorkflowReport) -> Double {
  report.dominant_frequency * 1000.0
}

///|
/// Return the number of boolean round-trip checks that passed.
pub fn report_roundtrip_count(report : WorkflowReport) -> Int {
  (if report.roundtrip_sac { 1 } else { 0 }) +
  (if report.roundtrip_miniseed { 1 } else { 0 })
}

///|
/// Return a short artifact label for archival logs.
pub fn report_artifact_label(report : WorkflowReport) -> String {
  report.station +
  "." +
  report.channel +
  "@" +
  report.processed_samples.to_string()
}

///|
/// Return whether the processed trace still has usable timing metadata.
pub fn report_timing_is_valid(report : WorkflowReport) -> Bool {
  report.processed_duration >= 0.0 && report.processed_samples >= 0
}

///|
/// Return a deterministic integer suitable for a smoke-test assertion.
pub fn report_fingerprint(report : WorkflowReport) -> Int {
  report.source_samples +
  report.processed_samples +
  report.event_count +
  report_roundtrip_count(report)
}

///|
/// Return the report's processed sample rate proxy in samples per second.
pub fn report_sample_density(report : WorkflowReport) -> Double {
  if report.processed_duration <= 0.0 {
    0.0
  } else {
    report.processed_samples.to_double() / report.processed_duration
  }
}

///|
/// Return whether the report has a positive processed sample density.
pub fn report_has_samples(report : WorkflowReport) -> Bool {
  report_sample_density(report) > 0.0
}
