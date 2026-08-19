///|
/// Configuration for an adaptive least-mean-squares noise canceller.
pub struct LmsConfig {
  step_size : Double
  filter_length : Int
  leakage : Double
}

///|
/// Return a stable default LMS configuration.
pub fn default_lms_config() -> LmsConfig {
  { step_size: 0.01, filter_length: 8, leakage: 0.0001 }
}

///|
/// Adaptive LMS noise cancellation using a reference signal.
pub fn lms_denoise(
  signal : Array[Double],
  reference : Array[Double],
  config : LmsConfig,
) -> Array[Double] {
  let result = Array::make(signal.length(), 0.0)
  let length = if config.filter_length < 1 { 1 } else { config.filter_length }
  let step = if config.step_size <= 0.0 { 0.001 } else { config.step_size }
  let leakage = if config.leakage < 0.0 { 0.0 } else { config.leakage }
  let weights = Array::make(length, 0.0)
  let input = Array::make(length, 0.0)
  for index = 0; index < signal.length(); index = index + 1 {
    for tap = length - 1; tap > 0; tap = tap - 1 {
      input[tap] = input[tap - 1]
    }
    input[0] = if index < reference.length() { reference[index] } else { 0.0 }
    let mut estimate = 0.0
    for tap = 0; tap < length; tap = tap + 1 {
      estimate += weights[tap] * input[tap]
    }
    let error = signal[index] - estimate
    result[index] = error
    for tap = 0; tap < length; tap = tap + 1 {
      weights[tap] = (1.0 - leakage) * weights[tap] + step * error * input[tap]
    }
  }
  result
}

///|
/// Smooth a signal with a scalar alpha-beta tracker.
pub fn alpha_beta_smooth(
  data : Array[Double],
  alpha : Double,
  beta : Double,
  sample_interval : Double,
) -> Array[Double] {
  let result = Array::make(data.length(), 0.0)
  if data.is_empty() {
    return result
  }
  let a = alpha.min(1.0).max(0.0)
  let b = beta.min(1.0).max(0.0)
  let dt = if sample_interval <= 0.0 { 1.0 } else { sample_interval }
  let mut position = data[0]
  let mut velocity = 0.0
  result[0] = position
  for index = 1; index < data.length(); index = index + 1 {
    let prediction = position + velocity * dt
    let residual = data[index] - prediction
    position = prediction + a * residual
    velocity = velocity + b * residual / dt
    result[index] = position
  }
  result
}

///|
/// Apply a robust Huber-style smoother to suppress isolated impulses.
pub fn huber_smooth(
  data : Array[Double],
  window_size : Int,
  threshold : Double,
  passes : Int,
) -> Array[Double] {
  let result = data.copy()
  let count = if passes < 1 { 1 } else { passes }
  let limit = if threshold <= 0.0 { 1.0 } else { threshold }
  for _ in 0..<count {
    let mean = moving_mean_samples(result, window_size)
    for index = 0; index < result.length(); index = index + 1 {
      let residual = result[index] - mean[index]
      let weight = if residual.abs() <= limit {
        1.0
      } else {
        limit / residual.abs()
      }
      result[index] = mean[index] + residual * weight
    }
  }
  result
}

///|
/// Apply a three-point second-difference smoother.
pub fn laplacian_smooth(
  data : Array[Double],
  strength : Double,
) -> Array[Double] {
  let result = data.copy()
  let amount = strength.max(0.0).min(0.5)
  if data.length() < 3 {
    return result
  }
  for index = 1; index < data.length() - 1; index = index + 1 {
    let curvature = data[index - 1] - 2.0 * data[index] + data[index + 1]
    result[index] = data[index] + amount * curvature
  }
  result
}

///|
/// Estimate the noise floor from a low percentile of absolute amplitudes.
pub fn noise_floor(data : Array[Double], percentile : Double) -> Double {
  let magnitudes = Array::new(capacity=data.length())
  for value in data {
    magnitudes.push(value.abs())
  }
  quantile_samples(magnitudes, percentile)
}

///|
/// Return a signal-to-noise weighted adaptive threshold.
pub fn adaptive_threshold(data : Array[Double], multiplier : Double) -> Double {
  noise_floor(data, 0.5) + multiplier.abs() * median_absolute_deviation(data)
}

///|
/// Apply an adaptive threshold while keeping the original sign.
pub fn threshold_denoise(
  data : Array[Double],
  threshold : Double,
) -> Array[Double] {
  let result = data.copy()
  let limit = threshold.abs()
  for index = 0; index < result.length(); index = index + 1 {
    if result[index].abs() < limit {
      result[index] = 0.0
    }
  }
  result
}
