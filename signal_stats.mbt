///|
/// Descriptive statistics for a seismic sample vector.
///
/// The structure deliberately keeps the aggregates together so a quality
/// report can compute them once and reuse the values in a CLI or a pipeline.
pub struct DescriptiveStats {
  count : Int
  minimum : Double
  maximum : Double
  mean : Double
  variance : Double
  standard_deviation : Double
  rms : Double
  energy : Double
  median : Double
  lower_quartile : Double
  upper_quartile : Double
  peak_to_peak : Double
  zero_crossings : Int
}

///|
/// Return a sorted copy of a sample array.
fn sorted_copy(data : Array[Double]) -> Array[Double] {
  let result = data.copy()
  result.sort()
  result
}

///|
/// Sum a signal without changing its order or contents.
pub fn sum_samples(data : Array[Double]) -> Double {
  let mut total = 0.0
  for value in data {
    total += value
  }
  total
}

///|
/// Return the arithmetic mean, or zero for an empty signal.
pub fn mean_samples(data : Array[Double]) -> Double {
  if data.is_empty() {
    0.0
  } else {
    sum_samples(data) / data.length().to_double()
  }
}

///|
/// Return the population variance of a sample vector.
pub fn variance_samples(data : Array[Double]) -> Double {
  let n = data.length()
  if n == 0 {
    return 0.0
  }
  let mean = mean_samples(data)
  let mut total = 0.0
  for value in data {
    let distance = value - mean
    total += distance * distance
  }
  total / n.to_double()
}

///|
/// Return the population standard deviation of a sample vector.
pub fn standard_deviation(data : Array[Double]) -> Double {
  variance_samples(data).sqrt()
}

///|
/// Return the root-mean-square amplitude.
pub fn rms_samples(data : Array[Double]) -> Double {
  if data.is_empty() {
    return 0.0
  }
  let mut total = 0.0
  for value in data {
    total += value * value
  }
  (total / data.length().to_double()).sqrt()
}

///|
/// Return the sum of squared amplitudes, useful as a discrete energy proxy.
pub fn signal_energy(data : Array[Double]) -> Double {
  let mut total = 0.0
  for value in data {
    total += value * value
  }
  total
}

///|
/// Return the smallest sample, or zero for an empty array.
pub fn minimum_sample(data : Array[Double]) -> Double {
  if data.is_empty() {
    return 0.0
  }
  let mut value = data[0]
  for current in data {
    if current < value {
      value = current
    }
  }
  value
}

///|
/// Return the largest sample, or zero for an empty array.
pub fn maximum_sample(data : Array[Double]) -> Double {
  if data.is_empty() {
    return 0.0
  }
  let mut value = data[0]
  for current in data {
    if current > value {
      value = current
    }
  }
  value
}

///|
/// Return the peak-to-peak range.
pub fn peak_to_peak_samples(data : Array[Double]) -> Double {
  maximum_sample(data) - minimum_sample(data)
}

///|
/// Return a linear-interpolated quantile of a sample vector.
///
/// `probability` is clamped to the closed interval `[0, 1]`. This makes the
/// helper useful for telemetry and user supplied quality thresholds without
/// turning a harmless configuration typo into an indexing panic.
pub fn quantile_samples(data : Array[Double], probability : Double) -> Double {
  if data.is_empty() {
    return 0.0
  }
  let values = sorted_copy(data)
  let p = if probability < 0.0 {
    0.0
  } else if probability > 1.0 {
    1.0
  } else {
    probability
  }
  let position = p * (values.length() - 1).to_double()
  let lower = position.to_int()
  let upper = if lower + 1 >= values.length() { lower } else { lower + 1 }
  let fraction = position - lower.to_double()
  values[lower] + fraction * (values[upper] - values[lower])
}

///|
/// Return the median of a sample vector.
pub fn median_samples(data : Array[Double]) -> Double {
  quantile_samples(data, 0.5)
}

///|
/// Count sign changes between adjacent samples.
pub fn zero_crossing_count(data : Array[Double]) -> Int {
  if data.length() < 2 {
    return 0
  }
  let mut count = 0
  let mut previous = data[0]
  for index = 1; index < data.length(); index = index + 1 {
    let current = data[index]
    if (previous < 0.0 && current >= 0.0) || (previous >= 0.0 && current < 0.0) {
      count += 1
    }
    if current != 0.0 {
      previous = current
    }
  }
  count
}

///|
/// Return the fraction of samples with an absolute value above `threshold`.
pub fn exceedance_ratio(data : Array[Double], threshold : Double) -> Double {
  if data.is_empty() {
    return 0.0
  }
  let mut count = 0
  for value in data {
    if value.abs() > threshold {
      count += 1
    }
  }
  count.to_double() / data.length().to_double()
}

///|
/// Return a complete descriptive summary.
pub fn describe_samples(data : Array[Double]) -> DescriptiveStats {
  let count = data.length()
  let mean = mean_samples(data)
  let variance = variance_samples(data)
  {
    count,
    minimum: minimum_sample(data),
    maximum: maximum_sample(data),
    mean,
    variance,
    standard_deviation: variance.sqrt(),
    rms: rms_samples(data),
    energy: signal_energy(data),
    median: median_samples(data),
    lower_quartile: quantile_samples(data, 0.25),
    upper_quartile: quantile_samples(data, 0.75),
    peak_to_peak: peak_to_peak_samples(data),
    zero_crossings: zero_crossing_count(data),
  }
}

///|
/// Return a copy normalized by its maximum absolute amplitude.
pub fn normalize_peak(data : Array[Double]) -> Array[Double] {
  let result = data.copy()
  let peak = maximum_sample(data).abs()
  let minimum = minimum_sample(data).abs()
  let scale = if peak > minimum { peak } else { minimum }
  if scale == 0.0 {
    return result
  }
  for index = 0; index < result.length(); index = index + 1 {
    result[index] = result[index] / scale
  }
  result
}

///|
/// Return a z-score normalized copy. Constant signals become all zeros.
pub fn normalize_zscore(data : Array[Double]) -> Array[Double] {
  let result = data.copy()
  let mean = mean_samples(data)
  let deviation = standard_deviation(data)
  if deviation == 0.0 {
    result.fill(0.0)
    return result
  }
  for index = 0; index < result.length(); index = index + 1 {
    result[index] = (result[index] - mean) / deviation
  }
  result
}

///|
/// Scale a signal into `[low, high]`, preserving constant signals at `low`.
pub fn scale_samples(
  data : Array[Double],
  low : Double,
  high : Double,
) -> Array[Double] {
  let result = data.copy()
  let source_low = minimum_sample(data)
  let source_high = maximum_sample(data)
  let width = source_high - source_low
  if width == 0.0 {
    result.fill(low)
    return result
  }
  for index = 0; index < result.length(); index = index + 1 {
    let t = (result[index] - source_low) / width
    result[index] = low + t * (high - low)
  }
  result
}

///|
/// Apply an absolute-value floor and ceiling in place.
pub fn clamp_samples(
  data : Array[Double],
  minimum : Double,
  maximum : Double,
) -> Unit {
  for index = 0; index < data.length(); index = index + 1 {
    if data[index] < minimum {
      data[index] = minimum
    } else if data[index] > maximum {
      data[index] = maximum
    }
  }
}

///|
/// Compute a centered moving mean in O(N) using a prefix sum.
pub fn moving_mean_samples(
  data : Array[Double],
  window_size : Int,
) -> Array[Double] {
  let n = data.length()
  let result = Array::make(n, 0.0)
  if n == 0 || window_size <= 0 {
    return result
  }
  let prefix = Array::make(n + 1, 0.0)
  for index = 0; index < n; index = index + 1 {
    prefix[index + 1] = prefix[index] + data[index]
  }
  let half = window_size / 2
  for index = 0; index < n; index = index + 1 {
    let start = if index < half { 0 } else { index - half }
    let end = if index + half + 1 > n { n } else { index + half + 1 }
    result[index] = (prefix[end] - prefix[start]) / (end - start).to_double()
  }
  result
}

///|
/// Compute a centered moving RMS in O(N) using a prefix sum of squares.
pub fn moving_rms_samples(
  data : Array[Double],
  window_size : Int,
) -> Array[Double] {
  let squares = Array::new(capacity=data.length())
  for value in data {
    squares.push(value * value)
  }
  let means = moving_mean_samples(squares, window_size)
  for index = 0; index < means.length(); index = index + 1 {
    means[index] = means[index].sqrt()
  }
  means
}
