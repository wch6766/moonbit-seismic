///|
/// Standard magnitude scales supported by the analysis helpers.
pub enum MagnitudeScale {
  LocalMagnitude
  BodyWaveMagnitude
  SurfaceWaveMagnitude
  DurationMagnitude
  MomentMagnitude
}

///|
/// Return scale constructors for configuration-driven callers.
pub fn scale_local() -> MagnitudeScale {
  LocalMagnitude
}

///|
pub fn scale_body() -> MagnitudeScale {
  BodyWaveMagnitude
}

///|
pub fn scale_surface() -> MagnitudeScale {
  SurfaceWaveMagnitude
}

///|
pub fn scale_duration() -> MagnitudeScale {
  DurationMagnitude
}

///|
pub fn scale_moment() -> MagnitudeScale {
  MomentMagnitude
}

///|
/// A magnitude estimate with provenance and uncertainty metadata.
pub struct MagnitudeEstimate {
  scale : MagnitudeScale
  value : Double
  amplitude : Double
  distance : Double
  uncertainty : Double
}

///|
/// Estimate duration magnitude from event duration in seconds and distance.
pub fn duration_magnitude(
  duration_seconds : Double,
  distance_km : Double,
) -> Result[Double, String] {
  if duration_seconds <= 0.0 || distance_km <= 0.0 {
    return Err("Duration and distance must be positive")
  }
  Ok(2.0 * @math.log10(duration_seconds) + 0.003 * distance_km - 0.87)
}

///|
/// Estimate scalar moment magnitude from seismic moment in N-m.
pub fn moment_magnitude(moment_newton_meter : Double) -> Result[Double, String] {
  if moment_newton_meter <= 0.0 {
    return Err("Seismic moment must be positive")
  }
  Ok(2.0 / 3.0 * (@math.log10(moment_newton_meter) - 9.1))
}

///|
/// Convert peak millimetres on an analog record to micrometres.
pub fn millimeters_to_micrometers(value : Double) -> Double {
  value * 1000.0
}

///|
/// Apply a simple geometric spreading correction.
pub fn geometric_spreading(amplitude : Double, distance_km : Double) -> Double {
  if distance_km <= 0.0 {
    0.0
  } else {
    amplitude * distance_km.sqrt()
  }
}

///|
/// Estimate local magnitude and retain the measurement inputs.
pub fn estimate_local(
  peak_amplitude_mm : Double,
  distance_km : Double,
  uncertainty : Double,
) -> Result[MagnitudeEstimate, String] {
  match local_magnitude(peak_amplitude_mm, distance_km) {
    Ok(value) =>
      Ok({
        scale: LocalMagnitude,
        value,
        amplitude: peak_amplitude_mm,
        distance: distance_km,
        uncertainty: if uncertainty < 0.0 {
          0.0
        } else {
          uncertainty
        },
      })
    Err(message) => Err(message)
  }
}

///|
/// Estimate a local magnitude directly from a processed trace.
pub fn magnitude_from_trace(
  trace : Trace,
  distance_km : Double,
) -> Result[MagnitudeEstimate, String] {
  estimate_local(wood_anderson_peak(trace), distance_km, 0.2)
}

///|
/// Return a magnitude estimate as a compact text record.
pub fn magnitude_summary(estimate : MagnitudeEstimate) -> String {
  let label = match estimate.scale {
    LocalMagnitude => "ML"
    BodyWaveMagnitude => "mb"
    SurfaceWaveMagnitude => "Ms"
    DurationMagnitude => "Md"
    MomentMagnitude => "Mw"
  }
  label +
  "=" +
  estimate.value.to_string() +
  " +/- " +
  estimate.uncertainty.to_string()
}

///|
/// Build a duration magnitude estimate.
pub fn estimate_duration(
  duration_seconds : Double,
  distance_km : Double,
) -> Result[MagnitudeEstimate, String] {
  match duration_magnitude(duration_seconds, distance_km) {
    Ok(value) =>
      Ok({
        scale: DurationMagnitude,
        value,
        amplitude: duration_seconds,
        distance: distance_km,
        uncertainty: 0.3,
      })
    Err(message) => Err(message)
  }
}

///|
/// Build a moment magnitude estimate.
pub fn estimate_moment(
  moment_newton_meter : Double,
) -> Result[MagnitudeEstimate, String] {
  match moment_magnitude(moment_newton_meter) {
    Ok(value) =>
      Ok({
        scale: MomentMagnitude,
        value,
        amplitude: moment_newton_meter,
        distance: 0.0,
        uncertainty: 0.2,
      })
    Err(message) => Err(message)
  }
}
