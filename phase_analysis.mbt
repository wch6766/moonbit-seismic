///|
/// A phase pick with confidence and local signal measurements.
pub struct PhaseMeasurement {
  sample_index : Int
  time : Double
  amplitude : Double
  phase : Double
  frequency : Double
  confidence : Double
}

///|
/// Return a wrapped phase in `[-PI, PI]`.
pub fn wrap_phase(phase : Double) -> Double {
  let mut value = phase
  while value > @math.PI {
    value = value - 2.0 * @math.PI
  }
  while value < -@math.PI {
    value = value + 2.0 * @math.PI
  }
  value
}

///|
/// Unwrap a phase sequence by adding whole turns at discontinuities.
pub fn unwrap_phase(phases : Array[Double]) -> Array[Double] {
  let result = phases.copy()
  if result.length() < 2 {
    return result
  }
  let mut offset = 0.0
  for index = 1; index < result.length(); index = index + 1 {
    let difference = result[index] - result[index - 1]
    if difference > @math.PI {
      offset = offset - 2.0 * @math.PI
    } else if difference < -@math.PI {
      offset = offset + 2.0 * @math.PI
    }
    result[index] = result[index] + offset
  }
  result
}

///|
/// Estimate instantaneous phase from a centered quadrature approximation.
pub fn instantaneous_phase(data : Array[Double]) -> Array[Double] {
  let result = Array::make(data.length(), 0.0)
  if data.is_empty() {
    return result
  }
  for index = 0; index < data.length(); index = index + 1 {
    let previous = if index == 0 { data[index] } else { data[index - 1] }
    let next = if index + 1 >= data.length() {
      data[index]
    } else {
      data[index + 1]
    }
    result[index] = @math.atan2(0.5 * (next - previous), data[index])
  }
  unwrap_phase(result)
}

///|
/// Estimate instantaneous frequency from an unwrapped phase sequence.
pub fn instantaneous_frequency(
  phase : Array[Double],
  sample_rate : Double,
) -> Array[Double] {
  let result = Array::make(phase.length(), 0.0)
  if phase.length() < 2 || sample_rate <= 0.0 {
    return result
  }
  result[0] = (phase[1] - phase[0]) * sample_rate / (2.0 * @math.PI)
  for index = 1; index < phase.length() - 1; index = index + 1 {
    result[index] = (phase[index + 1] - phase[index - 1]) *
      sample_rate /
      (4.0 * @math.PI)
  }
  result[phase.length() - 1] = (
      phase[phase.length() - 1] - phase[phase.length() - 2]
    ) *
    sample_rate /
    (2.0 * @math.PI)
  result
}

///|
/// Estimate local phase and frequency at a selected sample.
pub fn phase_at(
  data : Array[Double],
  sample_rate : Double,
  index : Int,
) -> PhaseMeasurement {
  let phases = instantaneous_phase(data)
  let frequencies = instantaneous_frequency(phases, sample_rate)
  let safe = if index < 0 {
    0
  } else if index >= data.length() {
    data.length() - 1
  } else {
    index
  }
  if data.is_empty() {
    {
      sample_index: 0,
      time: 0.0,
      amplitude: 0.0,
      phase: 0.0,
      frequency: 0.0,
      confidence: 0.0,
    }
  } else {
    let amplitude = data[safe].abs()
    let confidence = amplitude / (rms_samples(data) + 0.0000000001)
    {
      sample_index: safe,
      time: if sample_rate > 0.0 {
        safe.to_double() / sample_rate
      } else {
        0.0
      },
      amplitude,
      phase: phases[safe],
      frequency: frequencies[safe],
      confidence: confidence.min(1.0),
    }
  }
}

///|
/// Estimate a dominant period from zero crossings.
pub fn estimate_period(data : Array[Double], sample_rate : Double) -> Double {
  let crossings = zero_crossing_indices(data)
  if crossings.length() < 2 || sample_rate <= 0.0 {
    return 0.0
  }
  let mut total = 0.0
  for index = 1; index < crossings.length(); index = index + 1 {
    total += (crossings[index] - crossings[index - 1]).to_double()
  }
  let average_half_period = total / (crossings.length() - 1).to_double()
  2.0 * average_half_period / sample_rate
}

///|
/// Estimate phase coherence between two same-length signals.
pub fn phase_coherence(left : Array[Double], right : Array[Double]) -> Double {
  let count = if left.length() < right.length() {
    left.length()
  } else {
    right.length()
  }
  if count == 0 {
    return 0.0
  }
  let left_phase = instantaneous_phase(slice_samples(left, 0, count))
  let right_phase = instantaneous_phase(slice_samples(right, 0, count))
  let mut real = 0.0
  let mut imaginary = 0.0
  for index = 0; index < count; index = index + 1 {
    let difference = left_phase[index] - right_phase[index]
    real += @math.cos(difference)
    imaginary += @math.sin(difference)
  }
  (real * real + imaginary * imaginary).sqrt() / count.to_double()
}

///|
/// Pick the first high-confidence phase measurement in a window.
pub fn pick_phase(
  trace : Trace,
  minimum_amplitude : Double,
  minimum_confidence : Double,
) -> PhaseMeasurement? {
  let phase = instantaneous_phase(trace.data)
  let frequency = instantaneous_frequency(phase, trace.sampling_rate)
  for index = 0; index < trace.npts(); index = index + 1 {
    let amplitude = trace.data[index].abs()
    let confidence = amplitude / (rms_samples(trace.data) + 0.0000000001)
    if amplitude >= minimum_amplitude && confidence >= minimum_confidence {
      return Some({
        sample_index: index,
        time: sample_to_time(trace, index),
        amplitude,
        phase: phase[index],
        frequency: frequency[index],
        confidence: confidence.min(1.0),
      })
    }
  }
  None
}
