///|
/// A structured validation finding for a waveform or model.
pub struct ValidationFinding {
  code : String
  severity : Int
  message : String
  sample_index : Int?
}

///|
/// Validate a signal's amplitude and timing invariants.
pub fn validate_samples(
  data : Array[Double],
  sample_rate : Double,
  maximum_amplitude : Double,
) -> Array[ValidationFinding] {
  let findings = Array::new()
  if data.is_empty() {
    findings.push({
      code: "EMPTY",
      severity: 2,
      message: "Signal has no samples",
      sample_index: None,
    })
  }
  if sample_rate <= 0.0 {
    findings.push({
      code: "RATE",
      severity: 2,
      message: "Sample rate must be positive",
      sample_index: None,
    })
  }
  let limit = maximum_amplitude.abs()
  for index = 0; index < data.length(); index = index + 1 {
    if !is_valid_double(data[index]) {
      findings.push({
        code: "NONFINITE",
        severity: 2,
        message: "Invalid floating-point sample",
        sample_index: Some(index),
      })
    } else if limit > 0.0 && data[index].abs() > limit {
      findings.push({
        code: "RANGE",
        severity: 1,
        message: "Sample exceeds configured amplitude range",
        sample_index: Some(index),
      })
    }
  }
  findings
}

///|
/// Validate a trace and return all findings instead of stopping at the first.
pub fn validate_trace_findings(
  trace : Trace,
  maximum_amplitude : Double,
) -> Array[ValidationFinding] {
  validate_samples(trace.data, trace.sampling_rate, maximum_amplitude)
}

///|
/// Return the number of error-severity findings.
pub fn validation_error_count(findings : Array[ValidationFinding]) -> Int {
  findings.count_if(finding => finding.severity >= 2)
}

///|
/// Return the number of warning-severity findings.
pub fn validation_warning_count(findings : Array[ValidationFinding]) -> Int {
  findings.count_if(finding => finding.severity == 1)
}

///|
/// Return true when no error-severity finding exists.
pub fn validation_passes(findings : Array[ValidationFinding]) -> Bool {
  validation_error_count(findings) == 0
}

///|
/// Validate that a velocity model is compatible with a source depth.
pub fn validate_model_for_depth(
  model : VelocityModel,
  source_depth : Double,
) -> Array[ValidationFinding] {
  let findings = Array::new()
  match validate_velocity_model(model) {
    Ok(_) => ()
    Err(message) =>
      findings.push({ code: "MODEL", severity: 2, message, sample_index: None })
  }
  if source_depth < 0.0 {
    findings.push({
      code: "DEPTH",
      severity: 2,
      message: "Source depth cannot be negative",
      sample_index: None,
    })
  }
  findings
}

///|
/// Return a stable text representation of validation findings.
pub fn validation_findings_text(findings : Array[ValidationFinding]) -> String {
  let mut output = "severity,code,message,sample\n"
  for finding in findings {
    let sample = match finding.sample_index {
      Some(index) => index.to_string()
      None => ""
    }
    output = output +
      finding.severity.to_string() +
      "," +
      finding.code +
      "," +
      finding.message +
      "," +
      sample +
      "\n"
  }
  output
}

///|
/// Require a clean trace, returning the first human-readable issue.
pub fn require_valid_trace(
  trace : Trace,
  maximum_amplitude : Double,
) -> Result[Unit, String] {
  let findings = validate_trace_findings(trace, maximum_amplitude)
  if validation_passes(findings) {
    Ok(())
  } else {
    Err(validation_findings_text(findings))
  }
}

///|
/// Validate the round-trip size invariant for a trace's SAC encoding.
pub fn validate_sac_roundtrip_size(trace : Trace) -> Bool {
  let bytes = to_sac_bytes(trace, false)
  bytes.length() == sac_record_size(trace.npts())
}

///|
/// Validate the fixed-width invariant for a miniSEED encoding.
pub fn validate_miniseed_roundtrip_size(trace : Trace) -> Bool {
  let bytes = to_miniseed_bytes(trace, 1, true)
  bytes.length() == miniseed_record_size(trace.npts(), 64)
}
