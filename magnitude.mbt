///|
/// Find the peak absolute amplitude in a seismic trace.
pub fn peak_amplitude(trace : Trace) -> Double {
  let n = trace.npts()
  if n == 0 {
    return 0.0
  }
  let mut max_val = 0.0
  for v in trace.data {
    let abs_v = v.abs()
    if abs_v > max_val {
      max_val = abs_v
    }
  }
  max_val
}

///|
/// Calculate Richter Local Magnitude (M_L).
/// * `peak_amp_mm` : The peak amplitude measured on a simulated Wood-Anderson seismograph in mm.
/// * `distance_km` : The epicentral distance in km.
pub fn local_magnitude(
  peak_amp_mm : Double,
  distance_km : Double,
) -> Result[Double, String] {
  if distance_km <= 0.0 {
    return Err("Distance must be positive")
  }
  if peak_amp_mm <= 0.0 {
    return Err("Amplitude must be positive")
  }

  let log_a = @math.log10(peak_amp_mm)
  // Standard Richter calibration function: -log10(A_0)
  // Approximated for distances under 200 km
  let calib = 1.11 * @math.log10(distance_km) + 0.00189 * distance_km - 2.09
  Ok(log_a + calib)
}

///|
/// Calculate Body Wave Magnitude (M_b).
/// * `amplitude_micrometer` : Peak displacement amplitude in micrometers.
/// * `period_seconds` : Period of the wave in seconds.
/// * `distance_degrees` : Epicentral distance in degrees.
pub fn body_wave_magnitude(
  amplitude_micrometer : Double,
  period_seconds : Double,
  distance_degrees : Double,
) -> Result[Double, String] {
  if amplitude_micrometer <= 0.0 ||
    period_seconds <= 0.0 ||
    distance_degrees <= 0.0 {
    return Err("All parameters must be positive")
  }

  let ratio = amplitude_micrometer / period_seconds
  let log_at = @math.log10(ratio)
  let calib = 3.75 + 1.2 * @math.log10(distance_degrees)
  Ok(log_at + calib)
}

///|
/// Calculate Surface Wave Magnitude (M_s).
/// * `amplitude_micrometer` : Peak displacement amplitude in micrometers.
/// * `period_seconds` : Period of the wave in seconds.
/// * `distance_degrees` : Epicentral distance in degrees.
pub fn surface_wave_magnitude(
  amplitude_micrometer : Double,
  period_seconds : Double,
  distance_degrees : Double,
) -> Result[Double, String] {
  if amplitude_micrometer <= 0.0 ||
    period_seconds <= 0.0 ||
    distance_degrees <= 0.0 {
    return Err("All parameters must be positive")
  }

  let ratio = amplitude_micrometer / period_seconds
  let log_at = @math.log10(ratio)
  let calib = 1.66 * @math.log10(distance_degrees) + 3.3
  Ok(log_at + calib)
}
