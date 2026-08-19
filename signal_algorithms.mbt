///|
/// A detected local peak and its neighboring baseline.
pub struct SignalPeak {
  index : Int
  amplitude : Double
  prominence : Double
  width_samples : Int
}

///|
/// Return the median absolute deviation, a robust noise scale estimator.
pub fn median_absolute_deviation(data : Array[Double]) -> Double {
  if data.is_empty() {
    return 0.0
  }
  let center = median_samples(data)
  let distances = Array::new(capacity=data.length())
  for value in data {
    distances.push((value - center).abs())
  }
  median_samples(distances)
}

///|
/// Remove outliers with a median/MAD rule while preserving the signal length.
pub fn robust_clip(
  data : Array[Double],
  mad_multiplier : Double,
) -> Array[Double] {
  let result = data.copy()
  if data.is_empty() {
    return result
  }
  let center = median_samples(data)
  let mad = median_absolute_deviation(data)
  let scale = if mad == 0.0 { 0.0 } else { mad_multiplier.abs() * mad }
  if scale == 0.0 {
    return result
  }
  for index = 0; index < result.length(); index = index + 1 {
    if result[index] > center + scale {
      result[index] = center + scale
    } else if result[index] < center - scale {
      result[index] = center - scale
    }
  }
  result
}

///|
/// Find local maxima with a minimum height and separation.
pub fn find_peaks(
  data : Array[Double],
  minimum_height : Double,
  minimum_distance : Int,
) -> Array[SignalPeak] {
  let result = Array::new()
  if data.length() < 3 {
    return result
  }
  let distance = if minimum_distance < 1 { 1 } else { minimum_distance }
  let mut last_peak = -distance
  for index = 1; index < data.length() - 1; index = index + 1 {
    let value = data[index]
    if value >= minimum_height &&
      value >= data[index - 1] &&
      value >= data[index + 1] &&
      index - last_peak >= distance {
      let mut left = index
      while left > 0 && data[left - 1] <= value {
        left = left - 1
        if index - left >= distance {
          break
        }
      }
      let mut right = index
      while right + 1 < data.length() && data[right + 1] <= value {
        right = right + 1
        if right - index >= distance {
          break
        }
      }
      let baseline = 0.5 * (data[left] + data[right])
      result.push({
        index,
        amplitude: value,
        prominence: value - baseline,
        width_samples: right - left + 1,
      })
      last_peak = index
    }
  }
  result
}

///|
/// Return zero-crossing sample indices.
pub fn zero_crossing_indices(data : Array[Double]) -> Array[Int] {
  let result = Array::new()
  if data.length() < 2 {
    return result
  }
  for index = 1; index < data.length(); index = index + 1 {
    if (data[index - 1] < 0.0 && data[index] >= 0.0) ||
      (data[index - 1] >= 0.0 && data[index] < 0.0) {
      result.push(index)
    }
  }
  result
}

///|
/// Measure the duration for which a signal exceeds an absolute threshold.
pub fn threshold_duration(
  data : Array[Double],
  threshold : Double,
  sample_rate : Double,
) -> Double {
  if data.is_empty() || sample_rate <= 0.0 {
    return 0.0
  }
  let mut count = 0
  for value in data {
    if value.abs() >= threshold.abs() {
      count += 1
    }
  }
  count.to_double() / sample_rate
}

///|
/// Integrate only the positive part of a signal.
pub fn positive_area(data : Array[Double], sample_rate : Double) -> Double {
  if sample_rate <= 0.0 {
    return 0.0
  }
  let mut area = 0.0
  for value in data {
    if value > 0.0 {
      area += value / sample_rate
    }
  }
  area
}

///|
/// Integrate the absolute amplitude of a signal.
pub fn absolute_area(data : Array[Double], sample_rate : Double) -> Double {
  if sample_rate <= 0.0 {
    0.0
  } else {
    sum_samples(data.map(value => value.abs())) / sample_rate
  }
}

///|
/// Compute a normalized first-order difference of successive samples.
pub fn normalized_difference(data : Array[Double]) -> Array[Double] {
  let result = difference_characteristic(data)
  let peak = maximum_sample(result)
  if peak == 0.0 {
    return result
  }
  for index = 0; index < result.length(); index = index + 1 {
    result[index] = result[index] / peak
  }
  result
}

///|
/// Calculate the lag-one autocorrelation coefficient.
pub fn lag_one_autocorrelation(data : Array[Double]) -> Double {
  let values = autocorrelation(data, 1)
  if values.length() < 2 {
    0.0
  } else {
    values[1]
  }
}

///|
/// Return the fraction of energy in a selected index interval.
pub fn energy_fraction(data : Array[Double], start : Int, end : Int) -> Double {
  let total = signal_energy(data)
  if total == 0.0 {
    0.0
  } else {
    signal_energy(slice_samples(data, start, end)) / total
  }
}

///|
/// Smooth a signal repeatedly with a centered moving average.
pub fn repeated_smooth(
  data : Array[Double],
  window_size : Int,
  passes : Int,
) -> Array[Double] {
  let mut result = data.copy()
  let count = if passes < 0 { 0 } else { passes }
  for _ in 0..<count {
    result = moving_mean_samples(result, window_size)
  }
  result
}

///|
/// Return a robust threshold based on the median and MAD.
pub fn robust_threshold(data : Array[Double], multiplier : Double) -> Double {
  median_samples(data) +
  multiplier.abs() * 1.4826 * median_absolute_deviation(data)
}

///|
/// Return the nearest sample index to a physical time.
pub fn time_to_sample(trace : Trace, time : Double) -> Int {
  if trace.sampling_rate <= 0.0 {
    0
  } else {
    ((time - trace.start_time) * trace.sampling_rate).round().to_int()
  }
}

///|
/// Return the physical time represented by a sample index.
pub fn sample_to_time(trace : Trace, index : Int) -> Double {
  trace.start_time + index.to_double() * trace.delta()
}
