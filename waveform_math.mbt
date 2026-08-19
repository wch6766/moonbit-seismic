///|
/// Return a Hann-weighted dot product for two same-length windows.
pub fn weighted_dot(left : Array[Double], right : Array[Double]) -> Double {
  let count = if left.length() < right.length() {
    left.length()
  } else {
    right.length()
  }
  if count == 0 {
    return 0.0
  }
  let mut total = 0.0
  for index = 0; index < count; index = index + 1 {
    let weight = 0.5 -
      0.5 * @math.cos(2.0 * @math.PI * index.to_double() / count.to_double())
    total += weight * left[index] * right[index]
  }
  total
}

///|
/// Return the normalized weighted mean-square error.
pub fn weighted_mse(left : Array[Double], right : Array[Double]) -> Double {
  let count = if left.length() < right.length() {
    left.length()
  } else {
    right.length()
  }
  if count == 0 {
    return 0.0
  }
  let mut total = 0.0
  let mut weights = 0.0
  for index = 0; index < count; index = index + 1 {
    let weight = 0.5 -
      0.5 * @math.cos(2.0 * @math.PI * index.to_double() / count.to_double())
    let difference = left[index] - right[index]
    total += weight * difference * difference
    weights += weight
  }
  if weights == 0.0 {
    0.0
  } else {
    total / weights
  }
}

///|
/// Return a normalized root-mean-square error.
pub fn normalized_rmse(
  reference : Array[Double],
  candidate : Array[Double],
) -> Double {
  let denominator = rms_samples(reference)
  if denominator == 0.0 {
    rms_samples(candidate)
  } else {
    weighted_mse(reference, candidate).sqrt() / denominator
  }
}

///|
/// Apply a polynomial baseline of degree one to a copy.
pub fn remove_polynomial_baseline(
  data : Array[Double],
  sample_rate : Double,
) -> Array[Double] {
  let result = data.copy()
  if data.length() < 2 || sample_rate <= 0.0 {
    return result
  }
  let slope = baseline_slope(data)
  let intercept = mean_samples(data) -
    slope * (data.length() - 1).to_double() / 2.0
  for index = 0; index < result.length(); index = index + 1 {
    result[index] = result[index] - (intercept + slope * index.to_double())
  }
  result
}

///|
/// Return a stable sample index for a fractional position.
pub fn fractional_index(length : Int, fraction : Double) -> Int {
  if length <= 0 {
    return 0
  }
  let p = fraction.max(0.0).min(1.0)
  (p * (length - 1).to_double()).round().to_int()
}

///|
/// Interpolate a sample at a fractional index.
pub fn sample_at_fraction(data : Array[Double], fraction : Double) -> Double {
  if data.is_empty() {
    return 0.0
  }
  if data.length() == 1 {
    return data[0]
  }
  let p = fraction.max(0.0).min(1.0) * (data.length() - 1).to_double()
  let left = p.to_int()
  if left >= data.length() - 1 {
    data[data.length() - 1]
  } else {
    let t = p - left.to_double()
    data[left] + t * (data[left + 1] - data[left])
  }
}

///|
/// Compute a normalized crossfade between two waveforms.
pub fn crossfade(
  left : Array[Double],
  right : Array[Double],
  fraction : Double,
) -> Array[Double] {
  let count = if left.length() < right.length() {
    left.length()
  } else {
    right.length()
  }
  let result = Array::make(count, 0.0)
  let t = fraction.max(0.0).min(1.0)
  for index = 0; index < count; index = index + 1 {
    result[index] = (1.0 - t) * left[index] + t * right[index]
  }
  result
}

///|
/// Return a normalized triangular weighting vector.
pub fn triangular_weights(length : Int) -> Array[Double] {
  let result = Array::make(if length < 0 { 0 } else { length }, 0.0)
  if length <= 0 {
    return result
  }
  let center = (length - 1).to_double() / 2.0
  for index = 0; index < length; index = index + 1 {
    result[index] = if center == 0.0 {
      1.0
    } else {
      1.0 - (index.to_double() - center).abs() / (center + 1.0)
    }
  }
  result
}

///|
/// Apply arbitrary weights and return a weighted average.
pub fn weighted_average(
  data : Array[Double],
  weights : Array[Double],
) -> Double {
  let count = if data.length() < weights.length() {
    data.length()
  } else {
    weights.length()
  }
  let mut numerator = 0.0
  let mut denominator = 0.0
  for index = 0; index < count; index = index + 1 {
    numerator += data[index] * weights[index]
    denominator += weights[index]
  }
  if denominator == 0.0 {
    0.0
  } else {
    numerator / denominator
  }
}
