///|
/// Peak information extracted from a correlation curve.
pub struct CorrelationPeak {
  lag : Int
  coefficient : Double
  time_offset : Double
  sharpness : Double
}

///|
/// Return an unbiased autocorrelation curve for non-negative lags.
pub fn autocorrelation(data : Array[Double], max_lag : Int) -> Array[Double] {
  let result = Array::new()
  if data.is_empty() || max_lag < 0 {
    return result
  }
  let mean = mean_samples(data)
  let mut total_variance = 0.0
  for value in data {
    let centered = value - mean
    total_variance += centered * centered
  }
  let limit = if max_lag >= data.length() { data.length() - 1 } else { max_lag }
  for lag = 0; lag <= limit; lag = lag + 1 {
    let mut covariance = 0.0
    let count = data.length() - lag
    for index = 0; index < count; index = index + 1 {
      covariance += (data[index] - mean) * (data[index + lag] - mean)
    }
    result.push(
      if total_variance == 0.0 {
        0.0
      } else {
        covariance / total_variance
      },
    )
  }
  result
}

///|
/// Compare two equal-length windows using Pearson correlation.
pub fn pearson_correlation(
  left : Array[Double],
  right : Array[Double],
) -> Double {
  let count = if left.length() < right.length() {
    left.length()
  } else {
    right.length()
  }
  if count == 0 {
    return 0.0
  }
  let left_mean = mean_samples(slice_samples(left, 0, count))
  let right_mean = mean_samples(slice_samples(right, 0, count))
  let mut covariance = 0.0
  let mut left_energy = 0.0
  let mut right_energy = 0.0
  for index = 0; index < count; index = index + 1 {
    let a = left[index] - left_mean
    let b = right[index] - right_mean
    covariance += a * b
    left_energy += a * a
    right_energy += b * b
  }
  if left_energy == 0.0 || right_energy == 0.0 {
    0.0
  } else {
    covariance / (left_energy * right_energy).sqrt()
  }
}

///|
/// Locate a correlation peak and estimate its local sharpness.
pub fn correlation_peak(
  result : CorrelationResult,
  sample_rate : Double,
) -> CorrelationPeak {
  if result.coefficients.is_empty() {
    return { lag: 0, coefficient: 0.0, time_offset: 0.0, sharpness: 0.0 }
  }
  let mut index = 0
  for candidate = 1
      candidate < result.coefficients.length()
      candidate = candidate + 1 {
    if result.coefficients[candidate] > result.coefficients[index] {
      index = candidate
    }
  }
  let left = if index == 0 {
    result.coefficients[index]
  } else {
    result.coefficients[index - 1]
  }
  let right = if index + 1 >= result.coefficients.length() {
    result.coefficients[index]
  } else {
    result.coefficients[index + 1]
  }
  let sharpness = result.coefficients[index] - 0.5 * (left + right)
  let dt = if sample_rate > 0.0 { 1.0 / sample_rate } else { 0.0 }
  {
    lag: result.lags[index],
    coefficient: result.coefficients[index],
    time_offset: result.lags[index].to_double() * dt,
    sharpness,
  }
}

///|
/// Correlate every valid template position in a signal.
pub fn sliding_template_match(
  signal : Array[Double],
  template : Array[Double],
  step : Int,
) -> Array[CorrelationPeak] {
  let result = Array::new()
  if signal.is_empty() ||
    template.is_empty() ||
    template.length() > signal.length() ||
    step <= 0 {
    return result
  }
  let mut start = 0
  while start + template.length() <= signal.length() {
    let window = slice_samples(signal, start, start + template.length())
    result.push({
      lag: start,
      coefficient: pearson_correlation(window, template),
      time_offset: start.to_double(),
      sharpness: 0.0,
    })
    start = start + step
  }
  result
}

///|
/// Convert a correlation lag into seconds relative to a trace.
pub fn lag_seconds(lag : Int, sample_rate : Double) -> Double {
  if sample_rate <= 0.0 {
    0.0
  } else {
    lag.to_double() / sample_rate
  }
}

///|
/// Estimate a stable lag from several pairwise correlations.
pub fn median_lag(results : Array[CorrelationResult]) -> Int {
  if results.is_empty() {
    return 0
  }
  let lags = Array::new(capacity=results.length())
  for result in results {
    lags.push(result.max_lag)
  }
  lags.sort()
  lags[lags.length() / 2]
}

///|
/// Stack aligned traces by arithmetic mean.
pub fn stack_traces(traces : Array[Trace]) -> Result[Trace, String] {
  if traces.is_empty() {
    return Err("Cannot stack an empty trace collection")
  }
  let reference = traces[0]
  let data = Array::make(reference.npts(), 0.0)
  let mut used = 0
  for trace in traces {
    if trace.npts() == reference.npts() &&
      trace.sampling_rate == reference.sampling_rate {
      for index = 0; index < data.length(); index = index + 1 {
        data[index] += trace.data[index]
      }
      used += 1
    }
  }
  if used == 0 {
    return Err("No traces share the reference sample grid")
  }
  for index = 0; index < data.length(); index = index + 1 {
    data[index] = data[index] / used.to_double()
  }
  Ok(
    Trace::new(
      reference.station,
      reference.channel,
      reference.sampling_rate,
      reference.start_time,
      data,
    ),
  )
}

///|
/// Return the average correlation of a trace collection to a reference.
pub fn stack_coherence(reference : Trace, traces : Array[Trace]) -> Double {
  if traces.is_empty() {
    return 0.0
  }
  let mut total = 0.0
  let mut used = 0
  for trace in traces {
    if trace.npts() == reference.npts() &&
      trace.sampling_rate == reference.sampling_rate {
      total += pearson_correlation(reference.data, trace.data)
      used += 1
    }
  }
  if used == 0 {
    0.0
  } else {
    total / used.to_double()
  }
}

///|
/// Whiten a signal by subtracting its mean and dividing by standard deviation.
pub fn whiten_for_correlation(data : Array[Double]) -> Array[Double] {
  normalize_zscore(data)
}

///|
/// Calculate a normalized correlation over a bounded lag range.
pub fn bounded_correlation(
  signal : Array[Double],
  template : Array[Double],
  maximum_lag : Int,
  sample_rate : Double,
) -> CorrelationPeak {
  correlation_peak(cross_correlate(signal, template, maximum_lag), sample_rate)
}

///|
/// Return a similarity score in `[0, 1]` for a candidate template match.
pub fn similarity_score(coefficient : Double) -> Double {
  let value = coefficient.abs()
  if value < 0.0 {
    0.0
  } else if value > 1.0 {
    1.0
  } else {
    value
  }
}
