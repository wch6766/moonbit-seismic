///|
/// Deterministic noise models for reproducible synthetic waveforms.
pub enum NoiseModel {
  NoNoise
  UniformNoise
  SinusoidalNoise
  ColoredNoise
}

///|
/// Return the no-noise model.
pub fn noise_none() -> NoiseModel {
  NoNoise
}

///|
/// Return the deterministic uniform-noise model.
pub fn noise_uniform() -> NoiseModel {
  UniformNoise
}

///|
/// Return the sinusoidal-noise model.
pub fn noise_sinusoidal() -> NoiseModel {
  SinusoidalNoise
}

///|
/// Return the low-frequency colored-noise model.
pub fn noise_colored() -> NoiseModel {
  ColoredNoise
}

///|
/// Description of a synthetic seismic arrival.
pub struct SyntheticEvent {
  onset_sample : Int
  amplitude : Double
  dominant_frequency : Double
  duration_samples : Int
  phase : Double
}

///|
/// Construct a synthetic event descriptor.
pub fn SyntheticEvent::new(
  onset_sample : Int,
  amplitude : Double,
  dominant_frequency : Double,
  duration_samples : Int,
  phase : Double,
) -> SyntheticEvent {
  { onset_sample, amplitude, dominant_frequency, duration_samples, phase }
}

///|
/// Generate a sinusoid with an optional phase offset.
pub fn sine_wave(
  length : Int,
  frequency : Double,
  sample_rate : Double,
  amplitude : Double,
  phase : Double,
) -> Array[Double] {
  let result = Array::make(if length < 0 { 0 } else { length }, 0.0)
  if length <= 0 || sample_rate <= 0.0 {
    return result
  }
  for index = 0; index < length; index = index + 1 {
    let time = index.to_double() / sample_rate
    result[index] = amplitude *
      @math.sin(2.0 * @math.PI * frequency * time + phase)
  }
  result
}

///|
/// Generate a linear frequency sweep.
pub fn chirp_wave(
  length : Int,
  start_frequency : Double,
  end_frequency : Double,
  sample_rate : Double,
  amplitude : Double,
) -> Array[Double] {
  let result = Array::make(if length < 0 { 0 } else { length }, 0.0)
  if length <= 0 || sample_rate <= 0.0 {
    return result
  }
  let duration = length.to_double() / sample_rate
  let slope = if duration == 0.0 {
    0.0
  } else {
    (end_frequency - start_frequency) / duration
  }
  for index = 0; index < length; index = index + 1 {
    let time = index.to_double() / sample_rate
    let phase = 2.0 *
      @math.PI *
      (start_frequency * time + 0.5 * slope * time * time)
    result[index] = amplitude * @math.sin(phase)
  }
  result
}

///|
/// Generate a deterministic pseudo-random value in `[-1, 1]`.
fn lcg_noise(state : Int) -> (Int, Double) {
  let next = (state * 1103515245 + 12345) & 0x7fffffff
  let value = next.to_double() / 1073741824.0 - 1.0
  (next, value)
}

///|
/// Generate deterministic uniform noise for tests and examples.
pub fn deterministic_noise(
  length : Int,
  amplitude : Double,
  seed : Int,
) -> Array[Double] {
  let result = Array::make(if length < 0 { 0 } else { length }, 0.0)
  let mut state = if seed == 0 { 1 } else { seed.abs() }
  for index = 0; index < result.length(); index = index + 1 {
    let (next, value) = lcg_noise(state)
    state = next
    result[index] = amplitude * value
  }
  result
}

///|
/// Add a noise model to a signal in place.
pub fn add_noise(
  data : Array[Double],
  sample_rate : Double,
  model : NoiseModel,
  amplitude : Double,
  frequency : Double,
  seed : Int,
) -> Unit {
  match model {
    NoNoise => ()
    UniformNoise => {
      let noise = deterministic_noise(data.length(), amplitude, seed)
      for index = 0; index < data.length(); index = index + 1 {
        data[index] = data[index] + noise[index]
      }
    }
    SinusoidalNoise => {
      let noise = sine_wave(
        data.length(),
        frequency,
        sample_rate,
        amplitude,
        0.0,
      )
      for index = 0; index < data.length(); index = index + 1 {
        data[index] = data[index] + noise[index]
      }
    }
    ColoredNoise => {
      let noise = deterministic_noise(data.length(), amplitude, seed)
      let mut previous = 0.0
      for index = 0; index < data.length(); index = index + 1 {
        previous = 0.85 * previous + 0.15 * noise[index]
        data[index] = data[index] + previous
      }
    }
  }
}

///|
/// Add a tapered Ricker event at a sample offset.
pub fn add_synthetic_event(
  data : Array[Double],
  sample_rate : Double,
  event : SyntheticEvent,
) -> Unit {
  if sample_rate <= 0.0 || event.duration_samples <= 0 {
    return
  }
  let wavelet = ricker_wavelet(
    event.duration_samples,
    event.dominant_frequency,
    sample_rate,
  )
  for index = 0; index < wavelet.length(); index = index + 1 {
    let target = event.onset_sample + index
    if target >= 0 && target < data.length() {
      let envelope = if event.duration_samples <= 1 {
        1.0
      } else {
        let phase = index.to_double() / (event.duration_samples - 1).to_double()
        @math.sin(@math.PI * phase)
      }
      let carrier = @math.cos(
        event.phase +
        2.0 *
        @math.PI *
        event.dominant_frequency *
        index.to_double() /
        sample_rate,
      )
      data[target] = data[target] +
        event.amplitude * wavelet[index] * envelope * carrier
    }
  }
}

///|
/// Generate an impulse with an optional triangular spread.
pub fn impulse_signal(
  length : Int,
  index : Int,
  amplitude : Double,
  spread : Int,
) -> Array[Double] {
  let result = Array::make(if length < 0 { 0 } else { length }, 0.0)
  let radius = if spread < 0 { 0 } else { spread }
  if index < 0 || index >= result.length() {
    return result
  }
  for offset = -radius; offset <= radius; offset = offset + 1 {
    let target = index + offset
    if target >= 0 && target < result.length() {
      let weight = if radius == 0 {
        1.0
      } else {
        1.0 - offset.abs().to_double() / (radius + 1).to_double()
      }
      result[target] = amplitude * weight
    }
  }
  result
}

///|
/// Generate a step and optional exponential relaxation.
pub fn step_signal(
  length : Int,
  onset : Int,
  amplitude : Double,
  relaxation : Double,
) -> Array[Double] {
  let result = Array::make(if length < 0 { 0 } else { length }, 0.0)
  for index = 0; index < result.length(); index = index + 1 {
    if index >= onset {
      let age = (index - onset).to_double()
      result[index] = amplitude *
        (if relaxation > 0.0 {
          1.0 - @math.exp(-age / relaxation)
        } else {
          1.0
        })
    } else {
      ()
    }
  }
  result
}

///|
/// Add two signals, treating missing tail samples as zero.
pub fn mix_signals(
  left : Array[Double],
  right : Array[Double],
  left_weight : Double,
  right_weight : Double,
) -> Array[Double] {
  let length = if left.length() > right.length() {
    left.length()
  } else {
    right.length()
  }
  let result = Array::make(length, 0.0)
  for index = 0; index < length; index = index + 1 {
    let left_value = if index < left.length() { left[index] } else { 0.0 }
    let right_value = if index < right.length() { right[index] } else { 0.0 }
    result[index] = left_weight * left_value + right_weight * right_value
  }
  result
}

///|
/// Generate a complete reproducible trace for documentation and benchmarks.
pub fn synthetic_trace(
  station : String,
  channel : String,
  length : Int,
  sample_rate : Double,
  event : SyntheticEvent,
  noise_amplitude : Double,
) -> Trace {
  let data = deterministic_noise(length, noise_amplitude, 20260819)
  add_synthetic_event(data, sample_rate, event)
  Trace::new(station, channel, sample_rate, 0.0, data)
}

///|
/// Generate a cosine-tapered sinusoid for filter response checks.
pub fn tapered_sine(
  length : Int,
  frequency : Double,
  sample_rate : Double,
  amplitude : Double,
  taper_fraction : Double,
) -> Array[Double] {
  let result = sine_wave(length, frequency, sample_rate, amplitude, 0.0)
  cosine_taper(result, taper_fraction)
  result
}

///|
/// Return an amplitude envelope for a synthetic event shape.
pub fn synthetic_event_envelope(
  length : Int,
  onset : Int,
  decay : Double,
) -> Array[Double] {
  let result = Array::make(if length < 0 { 0 } else { length }, 0.0)
  for index = 0; index < result.length(); index = index + 1 {
    if index >= onset {
      let age = (index - onset).to_double()
      result[index] = if decay <= 0.0 { 1.0 } else { @math.exp(-age / decay) }
    }
  }
  result
}
