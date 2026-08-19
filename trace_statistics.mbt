///|
/// A fixed-bin amplitude histogram for waveform quality dashboards.
pub struct AmplitudeHistogram {
  lower : Double
  upper : Double
  counts : Array[Int]
  total : Int
}

///|
/// Build a symmetric amplitude histogram.
pub fn amplitude_histogram(
  data : Array[Double],
  lower : Double,
  upper : Double,
  bins : Int,
) -> AmplitudeHistogram {
  let count = if bins < 1 { 1 } else { bins }
  let counts = Array::make(count, 0)
  let low = if lower < upper { lower } else { upper }
  let high = if upper > lower { upper } else { lower + 1.0 }
  let width = (high - low) / count.to_double()
  let mut total = 0
  for value in data {
    if value >= low && value <= high {
      let raw = ((value - low) / width).to_int()
      let index = if raw < 0 {
        0
      } else if raw >= count {
        count - 1
      } else {
        raw
      }
      counts[index] = counts[index] + 1
      total += 1
    }
  }
  { lower: low, upper: high, counts, total }
}

///|
/// Return the probability mass of a histogram bin.
pub fn histogram_probability(
  histogram : AmplitudeHistogram,
  index : Int,
) -> Double {
  if index < 0 || index >= histogram.counts.length() || histogram.total == 0 {
    0.0
  } else {
    histogram.counts[index].to_double() / histogram.total.to_double()
  }
}

///|
/// Compute covariance for two aligned sample arrays.
pub fn covariance_samples(
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
  let mut total = 0.0
  for index = 0; index < count; index = index + 1 {
    total += (left[index] - left_mean) * (right[index] - right_mean)
  }
  total / count.to_double()
}

///|
/// Compute a symmetric covariance matrix for a small trace collection.
pub fn covariance_matrix(
  signals : Array[Array[Double]],
) -> Array[Array[Double]] {
  let matrix = Array::new(capacity=signals.length())
  for left in signals {
    let row = Array::new(capacity=signals.length())
    for right in signals {
      row.push(covariance_samples(left, right))
    }
    matrix.push(row)
  }
  matrix
}

///|
/// Return a robust median center and scale pair.
pub fn robust_center_scale(data : Array[Double]) -> (Double, Double) {
  (median_samples(data), 1.4826 * median_absolute_deviation(data))
}

///|
/// Return a clipped, centered signal using a robust scale.
pub fn robust_standardize(
  data : Array[Double],
  multiplier : Double,
) -> Array[Double] {
  let (center, scale) = robust_center_scale(data)
  let result = data.copy()
  let safe_scale = if scale <= 0.0 { 1.0 } else { scale }
  let limit = multiplier.abs()
  for index = 0; index < result.length(); index = index + 1 {
    let value = (result[index] - center) / safe_scale
    result[index] = if value > limit {
      limit
    } else if value < -limit {
      -limit
    } else {
      value
    }
  }
  result
}

///|
/// Compute the Ljung-Box-like sum of squared autocorrelations.
pub fn whiteness_statistic(data : Array[Double], max_lag : Int) -> Double {
  let curve = autocorrelation(data, max_lag)
  let mut total = 0.0
  for index = 1; index < curve.length(); index = index + 1 {
    total += curve[index] * curve[index]
  }
  total
}

///|
/// Return an empirical threshold at a requested false-positive quantile.
pub fn empirical_threshold(
  data : Array[Double],
  probability : Double,
) -> Double {
  let absolute = data.map(value => value.abs())
  quantile_samples(absolute, probability)
}

///|
/// Compute a rolling median baseline for drift monitoring.
pub fn rolling_median(data : Array[Double], window_size : Int) -> Array[Double] {
  let result = Array::make(data.length(), 0.0)
  let width = if window_size < 1 { 1 } else { window_size }
  let half = width / 2
  for index = 0; index < data.length(); index = index + 1 {
    let start = if index < half { 0 } else { index - half }
    let end = if index + half + 1 > data.length() {
      data.length()
    } else {
      index + half + 1
    }
    result[index] = median_samples(slice_samples(data, start, end))
  }
  result
}

///|
/// Return a normalized baseline residual.
pub fn baseline_residual(
  data : Array[Double],
  window_size : Int,
) -> Array[Double] {
  let baseline = rolling_median(data, window_size)
  let result = data.copy()
  for index = 0; index < result.length(); index = index + 1 {
    result[index] = result[index] - baseline[index]
  }
  result
}
