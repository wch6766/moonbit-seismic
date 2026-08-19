///|
/// A contiguous flat or saturated portion of a waveform.
pub struct SampleSegment {
  start_index : Int
  end_index : Int
  duration : Double
}

///|
/// Machine-readable quality findings for a trace.
pub enum QualityFlag {
  EmptyTrace
  InvalidSampleRate
  TooShort
  Flatline
  Clipped
  BaselineDrift
  LowSignalToNoise
  ExcessiveGaps
}

///|
/// Summary of deterministic waveform quality checks.
pub struct QualityReport {
  score : Double
  flags : Array[QualityFlag]
  sample_count : Int
  peak_amplitude : Double
  rms_amplitude : Double
  noise_rms : Double
  signal_to_noise : Double
  clipped_fraction : Double
  baseline_slope : Double
  flat_fraction : Double
}

///|
/// Return a stable finite check for the values representable in portable
/// MoonBit arithmetic. Infinity and NaN are rejected by the self-inequality.
fn is_valid_double(value : Double) -> Bool {
  value == value && value.abs() < 1.0e300
}

///|
/// Count invalid floating-point samples.
pub fn invalid_sample_count(data : Array[Double]) -> Int {
  let mut count = 0
  for value in data {
    if !is_valid_double(value) {
      count += 1
    }
  }
  count
}

///|
/// Detect samples at or beyond a symmetric amplitude limit.
pub fn clipped_fraction(data : Array[Double], limit : Double) -> Double {
  if data.is_empty() || limit <= 0.0 {
    return 0.0
  }
  let mut clipped = 0
  for value in data {
    if value.abs() >= limit {
      clipped += 1
    }
  }
  clipped.to_double() / data.length().to_double()
}

///|
/// Find contiguous runs with absolute amplitude no greater than a tolerance.
pub fn flat_segments(
  data : Array[Double],
  tolerance : Double,
  sample_rate : Double,
) -> Array[SampleSegment] {
  let result = Array::new()
  if data.is_empty() {
    return result
  }
  let threshold = if tolerance < 0.0 { 0.0 } else { tolerance }
  let mut active = false
  let mut start = 0
  for index = 0; index < data.length(); index = index + 1 {
    let flat = data[index].abs() <= threshold
    if flat && !active {
      active = true
      start = index
    } else if !flat && active {
      let end = index - 1
      result.push({
        start_index: start,
        end_index: end,
        duration: if sample_rate > 0.0 {
          (end - start + 1).to_double() / sample_rate
        } else {
          0.0
        },
      })
      active = false
    }
  }
  if active {
    let end = data.length() - 1
    result.push({
      start_index: start,
      end_index: end,
      duration: if sample_rate > 0.0 {
        (end - start + 1).to_double() / sample_rate
      } else {
        0.0
      },
    })
  }
  result
}

///|
/// Detect saturated runs using a symmetric clipping threshold.
pub fn clipping_segments(
  data : Array[Double],
  limit : Double,
  sample_rate : Double,
) -> Array[SampleSegment] {
  let result = Array::new()
  if data.is_empty() || limit <= 0.0 {
    return result
  }
  let mut active = false
  let mut start = 0
  for index = 0; index < data.length(); index = index + 1 {
    let clipped = data[index].abs() >= limit
    if clipped && !active {
      active = true
      start = index
    } else if !clipped && active {
      let end = index - 1
      result.push({
        start_index: start,
        end_index: end,
        duration: if sample_rate > 0.0 {
          (end - start + 1).to_double() / sample_rate
        } else {
          0.0
        },
      })
      active = false
    }
  }
  if active {
    let end = data.length() - 1
    result.push({
      start_index: start,
      end_index: end,
      duration: if sample_rate > 0.0 {
        (end - start + 1).to_double() / sample_rate
      } else {
        0.0
      },
    })
  }
  result
}

///|
/// Estimate a baseline slope in amplitude units per sample.
pub fn baseline_slope(data : Array[Double]) -> Double {
  let n = data.length()
  if n <= 1 {
    return 0.0
  }
  let mean_x = (n - 1).to_double() / 2.0
  let mean_y = mean_samples(data)
  let mut covariance = 0.0
  let mut variance = 0.0
  for index = 0; index < n; index = index + 1 {
    let dx = index.to_double() - mean_x
    covariance += dx * (data[index] - mean_y)
    variance += dx * dx
  }
  if variance == 0.0 {
    0.0
  } else {
    covariance / variance
  }
}

///|
/// Estimate noise RMS from first differences.
pub fn difference_noise_rms(data : Array[Double]) -> Double {
  if data.length() < 2 {
    return 0.0
  }
  let differences = Array::new(capacity=data.length() - 1)
  for index = 1; index < data.length(); index = index + 1 {
    differences.push(data[index] - data[index - 1])
  }
  standard_deviation(differences) / 2.0.sqrt()
}

///|
/// Estimate signal-to-noise ratio in decibels using difference noise.
pub fn signal_to_noise_db(data : Array[Double]) -> Double {
  let signal = rms_samples(data)
  let noise = difference_noise_rms(data)
  if signal <= 0.0 || noise <= 0.0 {
    0.0
  } else {
    20.0 * @math.log10(signal / noise)
  }
}

///|
/// Replace isolated spikes with the local median.
pub fn repair_spikes(
  data : Array[Double],
  window_size : Int,
  threshold : Double,
) -> Array[Double] {
  let result = data.copy()
  let median = median_filter(data, window_size)
  let limit = if threshold < 0.0 { 0.0 } else { threshold }
  for index = 0; index < result.length(); index = index + 1 {
    if (result[index] - median[index]).abs() > limit {
      result[index] = median[index]
    }
  }
  result
}

///|
/// Linearly interpolate samples in a short clipped run.
pub fn repair_clipped_samples(
  data : Array[Double],
  limit : Double,
) -> Array[Double] {
  let result = data.copy()
  if limit <= 0.0 {
    return result
  }
  for index = 0; index < result.length(); index = index + 1 {
    if result[index].abs() >= limit {
      let left = if index == 0 { result[index] } else { result[index - 1] }
      let right = if index + 1 >= result.length() {
        result[index]
      } else {
        result[index + 1]
      }
      result[index] = 0.5 * (left + right)
    }
  }
  result
}

///|
/// Build a quality report using conservative, documented defaults.
pub fn assess_quality(
  trace : Trace,
  clip_limit : Double,
  flat_tolerance : Double,
) -> QualityReport {
  let data = trace.data
  let count = data.length()
  let peak = maximum_sample(data).abs()
  let rms = rms_samples(data)
  let noise = difference_noise_rms(data)
  let snr = signal_to_noise_db(data)
  let clipped = clipped_fraction(data, clip_limit)
  let flats = flat_segments(data, flat_tolerance, trace.sampling_rate)
  let flat_samples = flats.fold(init=0, (total, segment) => {
    total + segment.end_index - segment.start_index + 1
  })
  let flat_ratio = if count == 0 {
    0.0
  } else {
    flat_samples.to_double() / count.to_double()
  }
  let slope = baseline_slope(data)
  let flags = Array::new()
  if count == 0 {
    flags.push(EmptyTrace)
  }
  if trace.sampling_rate <= 0.0 {
    flags.push(InvalidSampleRate)
  }
  if count > 0 && count < 16 {
    flags.push(TooShort)
  }
  if flat_ratio >= 0.5 && count > 0 {
    flags.push(Flatline)
  }
  if clipped >= 0.01 {
    flags.push(Clipped)
  }
  if slope.abs() > 0.01 * (peak + 1.0) && count > 2 {
    flags.push(BaselineDrift)
  }
  if snr < 3.0 && count > 1 {
    flags.push(LowSignalToNoise)
  }
  if invalid_sample_count(data) > 0 {
    flags.push(ExcessiveGaps)
  }
  let penalty = flags.length().to_double() * 0.15
  let score = if 1.0 - penalty < 0.0 { 0.0 } else { 1.0 - penalty }
  {
    score,
    flags,
    sample_count: count,
    peak_amplitude: peak,
    rms_amplitude: rms,
    noise_rms: noise,
    signal_to_noise: snr,
    clipped_fraction: clipped,
    baseline_slope: slope,
    flat_fraction: flat_ratio,
  }
}

///|
/// Validate metadata and sample values before entering a processing pipeline.
pub fn validate_trace(trace : Trace) -> Result[Unit, String] {
  if trace.sampling_rate <= 0.0 {
    return Err("Trace sampling rate must be positive")
  }
  if invalid_sample_count(trace.data) > 0 {
    return Err("Trace contains an invalid floating-point sample")
  }
  Ok(())
}

///|
/// Return true when no quality flag is present.
pub fn is_good_quality(report : QualityReport) -> Bool {
  report.flags.is_empty()
}
