///|
/// A sensitivity curve sampled at calibration frequencies.
pub struct CalibrationCurve {
  frequencies : Array[Double]
  sensitivities : Array[Double]
  unit : WaveQuantity
}

///|
/// Construct a calibration curve.
pub fn CalibrationCurve::new(
  frequencies : Array[Double],
  sensitivities : Array[Double],
  unit : WaveQuantity,
) -> Result[CalibrationCurve, String] {
  if frequencies.length() != sensitivities.length() || frequencies.is_empty() {
    return Err("Calibration frequency and sensitivity arrays must match")
  }
  for index = 0; index < frequencies.length(); index = index + 1 {
    if frequencies[index] <= 0.0 || sensitivities[index] <= 0.0 {
      return Err("Calibration frequencies and sensitivities must be positive")
    }
    if index > 0 && frequencies[index] <= frequencies[index - 1] {
      return Err("Calibration frequencies must be strictly increasing")
    }
  }
  Ok({ frequencies, sensitivities, unit })
}

///|
/// Interpolate the sensitivity of a calibration curve in log-frequency space.
pub fn CalibrationCurve::sensitivity(
  self : CalibrationCurve,
  frequency : Double,
) -> Double {
  if frequency <= self.frequencies[0] {
    return self.sensitivities[0]
  }
  let last = self.frequencies.length() - 1
  if frequency >= self.frequencies[last] {
    return self.sensitivities[last]
  }
  let mut left = 0
  let mut right = last
  while left < right - 1 {
    let middle = (left + right) / 2
    if self.frequencies[middle] <= frequency {
      left = middle
    } else {
      right = middle
    }
  }
  let x0 = @math.ln(self.frequencies[left])
  let x1 = @math.ln(self.frequencies[right])
  let t = (@math.ln(frequency) - x0) / (x1 - x0)
  self.sensitivities[left] +
  t * (self.sensitivities[right] - self.sensitivities[left])
}

///|
/// Apply a frequency-dependent calibration to a trace using a selected band.
pub fn calibrate_with_curve(
  trace : Trace,
  curve : CalibrationCurve,
  representative_frequency : Double,
) -> Trace {
  let sensitivity = curve.sensitivity(representative_frequency)
  if sensitivity == 0.0 {
    trace.clone()
  } else {
    let data = trace.data.copy()
    for index = 0; index < data.length(); index = index + 1 {
      data[index] = data[index] / sensitivity
    }
    trace_with_data(trace, data)
  }
}

///|
/// Create a constant sensitivity curve for synthetic examples.
pub fn constant_calibration(
  sensitivity : Double,
  unit : WaveQuantity,
) -> CalibrationCurve {
  match
    CalibrationCurve::new(
      [0.1, 1.0, 10.0],
      [sensitivity, sensitivity, sensitivity],
      unit,
    ) {
    Ok(curve) => curve
    Err(_) => { frequencies: [1.0], sensitivities: [1.0], unit }
  }
}
