///|
/// A time-indexed view of a signal window.
pub struct SampleWindow {
  start_index : Int
  end_index : Int
  start_time : Double
  end_time : Double
  data : Array[Double]
}

///|
/// Return the number of windows produced by a sliding configuration.
pub fn window_count(
  sample_count : Int,
  window_size : Int,
  step : Int,
  include_partial : Bool,
) -> Int {
  if sample_count <= 0 || window_size <= 0 || step <= 0 {
    return 0
  }
  if include_partial {
    (sample_count + step - 1) / step
  } else if sample_count < window_size {
    0
  } else {
    (sample_count - window_size) / step + 1
  }
}

///|
/// Copy a half-open range. Invalid bounds are clipped to the signal.
pub fn slice_samples(
  data : Array[Double],
  start : Int,
  end : Int,
) -> Array[Double] {
  let first = if start < 0 {
    0
  } else if start > data.length() {
    data.length()
  } else {
    start
  }
  let last = if end < first {
    first
  } else if end > data.length() {
    data.length()
  } else {
    end
  }
  let result = Array::new(capacity=last - first)
  for index = first; index < last; index = index + 1 {
    result.push(data[index])
  }
  result
}

///|
/// Extract a window and optionally pad a short tail with `fill`.
pub fn extract_window(
  data : Array[Double],
  start : Int,
  window_size : Int,
  fill : Double,
) -> Array[Double] {
  if window_size <= 0 {
    return []
  }
  let result = Array::make(window_size, fill)
  for offset = 0; offset < window_size; offset = offset + 1 {
    let index = start + offset
    if index >= 0 && index < data.length() {
      result[offset] = data[index]
    }
  }
  result
}

///|
/// Split a signal into overlapping or non-overlapping windows.
pub fn split_windows(
  data : Array[Double],
  window_size : Int,
  step : Int,
  include_partial : Bool,
  fill : Double,
) -> Array[Array[Double]] {
  let result = Array::new()
  if data.is_empty() || window_size <= 0 || step <= 0 {
    return result
  }
  let mut start = 0
  while start < data.length() {
    if !include_partial && start + window_size > data.length() {
      break
    }
    result.push(extract_window(data, start, window_size, fill))
    start = start + step
  }
  result
}

///|
/// Produce windows carrying sample and physical time bounds.
pub fn trace_windows(
  trace : Trace,
  window_size : Int,
  step : Int,
  include_partial : Bool,
) -> Array[SampleWindow] {
  let result = Array::new()
  if trace.npts() == 0 || window_size <= 0 || step <= 0 {
    return result
  }
  let mut start = 0
  while start < trace.npts() {
    if !include_partial && start + window_size > trace.npts() {
      break
    }
    let end = if start + window_size > trace.npts() {
      trace.npts()
    } else {
      start + window_size
    }
    let values = slice_samples(trace.data, start, end)
    result.push({
      start_index: start,
      end_index: end,
      start_time: trace.start_time + start.to_double() * trace.delta(),
      end_time: trace.start_time + (end - 1).to_double() * trace.delta(),
      data: values,
    })
    start = start + step
  }
  result
}

///|
/// Return a copy of a trace restricted to a sample range.
pub fn trace_slice(trace : Trace, start : Int, end : Int) -> Trace {
  let first = if start < 0 {
    0
  } else if start > trace.npts() {
    trace.npts()
  } else {
    start
  }
  let values = slice_samples(trace.data, first, end)
  Trace::new(
    trace.station,
    trace.channel,
    trace.sampling_rate,
    trace.start_time + first.to_double() * trace.delta(),
    values,
  )
}

///|
/// Add constant padding around a signal.
pub fn pad_samples(
  data : Array[Double],
  left : Int,
  right : Int,
  fill : Double,
) -> Array[Double] {
  let left_count = if left < 0 { 0 } else { left }
  let right_count = if right < 0 { 0 } else { right }
  let result = Array::make(left_count, fill)
  for value in data {
    result.push(value)
  }
  for _ in 0..<right_count {
    result.push(fill)
  }
  result
}

///|
/// Generate a Hann analysis window.
pub fn hann_window(length : Int) -> Array[Double] {
  let result = Array::make(if length < 0 { 0 } else { length }, 0.0)
  if length <= 0 {
    return result
  }
  if length == 1 {
    result[0] = 1.0
    return result
  }
  for index = 0; index < length; index = index + 1 {
    result[index] = 0.5 -
      0.5 *
      @math.cos(2.0 * @math.PI * index.to_double() / (length - 1).to_double())
  }
  result
}

///|
/// Generate a Hamming analysis window.
pub fn hamming_window(length : Int) -> Array[Double] {
  let result = Array::make(if length < 0 { 0 } else { length }, 0.0)
  if length <= 0 {
    return result
  }
  if length == 1 {
    result[0] = 1.0
    return result
  }
  for index = 0; index < length; index = index + 1 {
    result[index] = 0.54 -
      0.46 *
      @math.cos(2.0 * @math.PI * index.to_double() / (length - 1).to_double())
  }
  result
}

///|
/// Generate a Blackman analysis window.
pub fn blackman_window(length : Int) -> Array[Double] {
  let result = Array::make(if length < 0 { 0 } else { length }, 0.0)
  if length <= 0 {
    return result
  }
  if length == 1 {
    result[0] = 1.0
    return result
  }
  for index = 0; index < length; index = index + 1 {
    let phase = 2.0 * @math.PI * index.to_double() / (length - 1).to_double()
    result[index] = 0.42 -
      0.5 * @math.cos(phase) +
      0.08 * @math.cos(2.0 * phase)
  }
  result
}

///|
/// Apply a supplied taper/window coefficient array in place.
pub fn apply_window(data : Array[Double], window : Array[Double]) -> Unit {
  let count = if data.length() < window.length() {
    data.length()
  } else {
    window.length()
  }
  for index = 0; index < count; index = index + 1 {
    data[index] = data[index] * window[index]
  }
}

///|
/// Apply a cosine taper to the first and last fraction of a signal.
pub fn cosine_taper(data : Array[Double], fraction : Double) -> Unit {
  let n = data.length()
  if n < 2 {
    return
  }
  let f = if fraction < 0.0 {
    0.0
  } else if fraction > 0.5 {
    0.5
  } else {
    fraction
  }
  let edge = (f * n.to_double()).to_int()
  if edge <= 0 {
    return
  }
  for index = 0; index < edge; index = index + 1 {
    let weight = 0.5 -
      0.5 * @math.cos(@math.PI * index.to_double() / edge.to_double())
    data[index] = data[index] * weight
    let mirror = n - 1 - index
    data[mirror] = data[mirror] * weight
  }
}

///|
/// Integrate a uniformly sampled signal with the trapezoidal rule.
pub fn cumulative_integral(
  data : Array[Double],
  sampling_rate : Double,
) -> Array[Double] {
  let result = Array::make(data.length(), 0.0)
  if data.is_empty() || sampling_rate <= 0.0 {
    return result
  }
  let dt = 1.0 / sampling_rate
  for index = 1; index < data.length(); index = index + 1 {
    result[index] = result[index - 1] +
      0.5 * (data[index - 1] + data[index]) * dt
  }
  result
}

///|
/// Differentiate a uniformly sampled signal with a centered finite difference.
pub fn differentiate_samples(
  data : Array[Double],
  sampling_rate : Double,
) -> Array[Double] {
  let result = Array::make(data.length(), 0.0)
  if data.length() < 2 || sampling_rate <= 0.0 {
    return result
  }
  let dt = 1.0 / sampling_rate
  result[0] = (data[1] - data[0]) / dt
  for index = 1; index < data.length() - 1; index = index + 1 {
    result[index] = (data[index + 1] - data[index - 1]) / (2.0 * dt)
  }
  result[data.length() - 1] = (
      data[data.length() - 1] - data[data.length() - 2]
    ) /
    dt
  result
}

///|
/// Downsample by taking every `factor`-th value.
pub fn decimate_samples(data : Array[Double], factor : Int) -> Array[Double] {
  let result = Array::new()
  if factor <= 0 {
    return result
  }
  for index = 0; index < data.length(); index = index + factor {
    result.push(data[index])
  }
  result
}

///|
/// Resample a signal with linear interpolation and a new sample rate.
pub fn resample_linear(
  data : Array[Double],
  old_rate : Double,
  new_rate : Double,
) -> Array[Double] {
  let result = Array::new()
  if data.is_empty() || old_rate <= 0.0 || new_rate <= 0.0 {
    return result
  }
  if data.length() == 1 {
    result.push(data[0])
    return result
  }
  let duration = (data.length() - 1).to_double() / old_rate
  let count = (duration * new_rate).to_int() + 1
  for index = 0; index < count; index = index + 1 {
    let position = index.to_double() * old_rate / new_rate
    let left = position.to_int()
    if left >= data.length() - 1 {
      result.push(data[data.length() - 1])
    } else {
      let fraction = position - left.to_double()
      result.push(data[left] + fraction * (data[left + 1] - data[left]))
    }
  }
  result
}

///|
/// Resample a trace while preserving its station metadata and start time.
pub fn resample_trace(trace : Trace, new_rate : Double) -> Trace {
  let values = resample_linear(trace.data, trace.sampling_rate, new_rate)
  Trace::new(
    trace.station,
    trace.channel,
    if new_rate > 0.0 {
      new_rate
    } else {
      trace.sampling_rate
    },
    trace.start_time,
    values,
  )
}

///|
/// Convolve a signal with a finite impulse response kernel in `same` mode.
pub fn convolve_same(
  data : Array[Double],
  kernel : Array[Double],
) -> Array[Double] {
  let result = Array::make(data.length(), 0.0)
  if data.is_empty() || kernel.is_empty() {
    return result
  }
  let center = kernel.length() / 2
  for output = 0; output < data.length(); output = output + 1 {
    let mut total = 0.0
    for tap = 0; tap < kernel.length(); tap = tap + 1 {
      let source = output + tap - center
      if source >= 0 && source < data.length() {
        total += data[source] * kernel[tap]
      }
    }
    result[output] = total
  }
  result
}
