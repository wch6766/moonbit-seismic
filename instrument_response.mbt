///|
/// Physical quantity represented by a calibrated waveform.
pub enum WaveQuantity {
  Counts
  Displacement
  Velocity
  Acceleration
}

///|
/// Construct common quantity tags.
pub fn quantity_counts() -> WaveQuantity {
  Counts
}

///|
/// Construct the displacement quantity tag.
pub fn quantity_displacement() -> WaveQuantity {
  Displacement
}

///|
/// Construct the velocity quantity tag.
pub fn quantity_velocity() -> WaveQuantity {
  Velocity
}

///|
/// Construct the acceleration quantity tag.
pub fn quantity_acceleration() -> WaveQuantity {
  Acceleration
}

///|
/// One real-valued pole/zero response stage.
pub struct ResponseStage {
  gain : Double
  poles : Array[Double]
  zeros : Array[Double]
  corner_frequency : Double
}

///|
/// A portable response description for calibration and unit conversion.
pub struct InstrumentResponse {
  name : String
  input_quantity : WaveQuantity
  output_quantity : WaveQuantity
  sensitivity : Double
  reference_frequency : Double
  stages : Array[ResponseStage]
}

///|
/// Build a response stage.
pub fn ResponseStage::new(
  gain : Double,
  poles : Array[Double],
  zeros : Array[Double],
  corner_frequency : Double,
) -> ResponseStage {
  { gain, poles, zeros, corner_frequency }
}

///|
/// Build an instrument response.
pub fn InstrumentResponse::new(
  name : String,
  input_quantity : WaveQuantity,
  output_quantity : WaveQuantity,
  sensitivity : Double,
  reference_frequency : Double,
  stages : Array[ResponseStage],
) -> InstrumentResponse {
  {
    name,
    input_quantity,
    output_quantity,
    sensitivity,
    reference_frequency,
    stages,
  }
}

///|
/// Validate response metadata before applying calibration.
pub fn validate_response(response : InstrumentResponse) -> Result[Unit, String] {
  if response.sensitivity <= 0.0 || response.reference_frequency <= 0.0 {
    return Err(
      "Instrument sensitivity and reference frequency must be positive",
    )
  }
  for stage in response.stages {
    if stage.gain == 0.0 || stage.corner_frequency <= 0.0 {
      return Err(
        "Response stages require non-zero gain and positive corner frequency",
      )
    }
  }
  Ok(())
}

///|
/// Return the scalar gain of one response stage at a frequency.
pub fn stage_gain(stage : ResponseStage, frequency : Double) -> Double {
  if frequency <= 0.0 {
    return stage.gain
  }
  let ratio = frequency / stage.corner_frequency
  let mut gain = stage.gain
  for pole in stage.poles {
    gain = gain /
      (1.0 +
      @math.pow(frequency / (stage.corner_frequency * (pole.abs() + 1.0)), 2.0)).sqrt()
  }
  for zero in stage.zeros {
    gain = gain *
      (1.0 +
      @math.pow(frequency / (stage.corner_frequency * (zero.abs() + 1.0)), 2.0)).sqrt()
  }
  gain * (1.0 + ratio * ratio).sqrt() / (1.0 + ratio * ratio).sqrt()
}

///|
/// Evaluate the scalar response gain at a frequency.
pub fn response_gain(
  response : InstrumentResponse,
  frequency : Double,
) -> Double {
  let mut gain = response.sensitivity
  for stage in response.stages {
    gain = gain * stage_gain(stage, frequency)
  }
  gain
}

///|
/// Convert sample units by the specified quantity relationship.
pub fn quantity_conversion(
  source : WaveQuantity,
  target : WaveQuantity,
  sample_rate : Double,
) -> Array[Double] {
  let result = Array::make(3, 0.0)
  let dt = if sample_rate > 0.0 { 1.0 / sample_rate } else { 0.0 }
  match (source, target) {
    (Counts, Counts) => result[0] = 1.0
    (Displacement, Displacement) => result[0] = 1.0
    (Velocity, Velocity) => result[0] = 1.0
    (Acceleration, Acceleration) => result[0] = 1.0
    (Displacement, Velocity) => result[1] = 1.0 / dt
    (Velocity, Displacement) => result[2] = dt
    (Velocity, Acceleration) => result[1] = 1.0 / dt
    (Acceleration, Velocity) => result[2] = dt
    _ => result[0] = 1.0
  }
  result
}

///|
/// Apply a calibration sensitivity to a trace copy.
pub fn calibrate_trace(
  trace : Trace,
  response : InstrumentResponse,
  frequency : Double,
) -> Trace {
  let gain = response_gain(response, frequency)
  if gain == 0.0 {
    return trace.clone()
  }
  let data = trace.data.copy()
  for index = 0; index < data.length(); index = index + 1 {
    data[index] = data[index] / gain
  }
  trace_with_data(trace, data)
}

///|
/// Apply a scalar response gain to a trace copy.
pub fn apply_instrument_gain(
  trace : Trace,
  response : InstrumentResponse,
  frequency : Double,
) -> Trace {
  let gain = response_gain(response, frequency)
  let data = trace.data.copy()
  for index = 0; index < data.length(); index = index + 1 {
    data[index] = data[index] * gain
  }
  trace_with_data(trace, data)
}

///|
/// Convert a displacement-like signal to velocity with finite differences.
pub fn displacement_to_velocity(trace : Trace) -> Trace {
  differentiate_trace(trace)
}

///|
/// Convert a velocity-like signal to displacement by integration.
pub fn velocity_to_displacement(trace : Trace) -> Trace {
  integrate_trace(trace)
}

///|
/// Convert acceleration to velocity and then displacement.
pub fn acceleration_to_displacement(trace : Trace) -> Trace {
  integrate_trace(integrate_trace(trace))
}

///|
/// Apply a simple Wood-Anderson-like band response for local magnitude work.
pub fn wood_anderson_response(trace : Trace) -> Trace {
  let filtered = trace.clone()
  if filtered.sampling_rate > 2.0 {
    highpass(filtered, 0.8, false)
    lowpass(filtered, 20.0.min(filtered.sampling_rate / 2.1), false)
  }
  filtered
}

///|
/// Estimate Wood-Anderson peak amplitude from a trace.
pub fn wood_anderson_peak(trace : Trace) -> Double {
  peak_amplitude(wood_anderson_response(trace))
}

///|
/// Return a standard response for a broadband teaching/example channel.
pub fn generic_broadband_response() -> InstrumentResponse {
  InstrumentResponse::new("generic-broadband", Counts, Displacement, 1.0, 1.0, [
    ResponseStage::new(1.0, [0.0, 0.0], [], 1.0),
  ])
}
