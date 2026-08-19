///|
/// Window choices for finite impulse response design.
pub enum FirWindow {
  Rectangular
  Hann
  Hamming
  Blackman
}

///|
/// Return the rectangular window choice.
pub fn fir_rectangular() -> FirWindow {
  Rectangular
}

///|
/// Return the Hann window choice.
pub fn fir_hann() -> FirWindow {
  Hann
}

///|
/// Return the Hamming window choice.
pub fn fir_hamming() -> FirWindow {
  Hamming
}

///|
/// Return the Blackman window choice.
pub fn fir_blackman() -> FirWindow {
  Blackman
}

///|
/// A linear-phase finite impulse response filter.
pub struct FirFilter {
  coefficients : Array[Double]
  sample_rate : Double
  low_cut : Double
  high_cut : Double
  group_delay : Double
  window : FirWindow
}

///|
/// Return a window coefficient for a normalized tap position.
fn fir_window_value(window : FirWindow, index : Int, length : Int) -> Double {
  if length <= 1 {
    return 1.0
  }
  let phase = 2.0 * @math.PI * index.to_double() / (length - 1).to_double()
  match window {
    Rectangular => 1.0
    Hann => 0.5 - 0.5 * @math.cos(phase)
    Hamming => 0.54 - 0.46 * @math.cos(phase)
    Blackman => 0.42 - 0.5 * @math.cos(phase) + 0.08 * @math.cos(2.0 * phase)
  }
}

///|
/// Return the normalized sinc function.
fn normalized_sinc(value : Double) -> Double {
  if value.abs() < 0.0000000001 {
    1.0
  } else {
    @math.sin(@math.PI * value) / (@math.PI * value)
  }
}

///|
/// Validate and normalize an FIR tap count.
fn normalized_tap_count(taps : Int) -> Int {
  let minimum = if taps < 3 { 3 } else { taps }
  if minimum % 2 == 0 {
    minimum + 1
  } else {
    minimum
  }
}

///|
/// Construct a low-pass FIR filter with a windowed-sinc kernel.
pub fn design_fir_lowpass(
  sample_rate : Double,
  cutoff : Double,
  taps : Int,
  window : FirWindow,
) -> Result[FirFilter, String] {
  if sample_rate <= 0.0 || cutoff <= 0.0 || cutoff >= sample_rate / 2.0 {
    return Err("FIR low-pass requires 0 < cutoff < Nyquist")
  }
  let count = normalized_tap_count(taps)
  let coefficients = Array::make(count, 0.0)
  let center = count / 2
  let normalized_cutoff = cutoff / sample_rate
  for index = 0; index < count; index = index + 1 {
    let distance = (index - center).to_double()
    let ideal = 2.0 *
      normalized_cutoff *
      normalized_sinc(2.0 * normalized_cutoff * distance)
    coefficients[index] = ideal * fir_window_value(window, index, count)
  }
  let total = sum_samples(coefficients)
  if total != 0.0 {
    for index = 0; index < count; index = index + 1 {
      coefficients[index] = coefficients[index] / total
    }
  }
  Ok({
    coefficients,
    sample_rate,
    low_cut: 0.0,
    high_cut: cutoff,
    group_delay: center.to_double() / sample_rate,
    window,
  })
}

///|
/// Construct a high-pass FIR filter by spectral inversion.
pub fn design_fir_highpass(
  sample_rate : Double,
  cutoff : Double,
  taps : Int,
  window : FirWindow,
) -> Result[FirFilter, String] {
  let lowpass = design_fir_lowpass(sample_rate, cutoff, taps, window)
  match lowpass {
    Err(message) => Err(message)
    Ok(filter) => {
      let coefficients = filter.coefficients.copy()
      let center = coefficients.length() / 2
      for index = 0; index < coefficients.length(); index = index + 1 {
        coefficients[index] = -coefficients[index]
      }
      coefficients[center] = coefficients[center] + 1.0
      Ok({
        coefficients,
        sample_rate,
        low_cut: cutoff,
        high_cut: sample_rate / 2.0,
        group_delay: filter.group_delay,
        window,
      })
    }
  }
}

///|
/// Construct a band-pass FIR filter as the difference of two low-pass kernels.
pub fn design_fir_bandpass(
  sample_rate : Double,
  low_cut : Double,
  high_cut : Double,
  taps : Int,
  window : FirWindow,
) -> Result[FirFilter, String] {
  if sample_rate <= 0.0 ||
    low_cut <= 0.0 ||
    high_cut <= low_cut ||
    high_cut >= sample_rate / 2.0 {
    return Err("FIR band-pass requires 0 < low < high < Nyquist")
  }
  let count = normalized_tap_count(taps)
  let coefficients = Array::make(count, 0.0)
  let center = count / 2
  let low = low_cut / sample_rate
  let high = high_cut / sample_rate
  for index = 0; index < count; index = index + 1 {
    let distance = (index - center).to_double()
    let upper = 2.0 * high * normalized_sinc(2.0 * high * distance)
    let lower = 2.0 * low * normalized_sinc(2.0 * low * distance)
    coefficients[index] = (upper - lower) *
      fir_window_value(window, index, count)
  }
  Ok({
    coefficients,
    sample_rate,
    low_cut,
    high_cut,
    group_delay: center.to_double() / sample_rate,
    window,
  })
}

///|
/// Construct a band-stop FIR filter by spectral inversion of a band-pass.
pub fn design_fir_bandstop(
  sample_rate : Double,
  low_cut : Double,
  high_cut : Double,
  taps : Int,
  window : FirWindow,
) -> Result[FirFilter, String] {
  let bandpass = design_fir_bandpass(
    sample_rate, low_cut, high_cut, taps, window,
  )
  match bandpass {
    Err(message) => Err(message)
    Ok(filter) => {
      let coefficients = filter.coefficients.copy()
      let center = coefficients.length() / 2
      for index = 0; index < coefficients.length(); index = index + 1 {
        coefficients[index] = -coefficients[index]
      }
      coefficients[center] = coefficients[center] + 1.0
      Ok({
        coefficients,
        sample_rate,
        low_cut,
        high_cut,
        group_delay: filter.group_delay,
        window,
      })
    }
  }
}

///|
/// Apply an FIR filter with edge samples reflected at the boundary.
pub fn apply_fir(data : Array[Double], filter : FirFilter) -> Array[Double] {
  let result = Array::make(data.length(), 0.0)
  if data.is_empty() || filter.coefficients.is_empty() {
    return result
  }
  let center = filter.coefficients.length() / 2
  for output = 0; output < data.length(); output = output + 1 {
    let mut total = 0.0
    for tap = 0; tap < filter.coefficients.length(); tap = tap + 1 {
      let raw = output + tap - center
      let source = if raw < 0 {
        -raw
      } else if raw >= data.length() {
        2 * data.length() - raw - 2
      } else {
        raw
      }
      if source >= 0 && source < data.length() {
        total += filter.coefficients[tap] * data[source]
      }
    }
    result[output] = total
  }
  result
}

///|
/// Apply an FIR filter forward and backward for zero-phase output.
pub fn apply_fir_zero_phase(
  data : Array[Double],
  filter : FirFilter,
) -> Array[Double] {
  let forward = apply_fir(data, filter)
  forward.rev_in_place()
  let backward = apply_fir(forward, filter)
  backward.rev_in_place()
  backward
}

///|
/// Filter a trace with an FIR filter and return a new trace.
pub fn filter_trace_fir(
  trace : Trace,
  filter : FirFilter,
  zero_phase : Bool,
) -> Trace {
  let values = if zero_phase {
    apply_fir_zero_phase(trace.data, filter)
  } else {
    apply_fir(trace.data, filter)
  }
  Trace::new(
    trace.station,
    trace.channel,
    trace.sampling_rate,
    trace.start_time,
    values,
  )
}

///|
/// A robust median filter with clipped edge windows.
pub fn median_filter(data : Array[Double], window_size : Int) -> Array[Double] {
  let result = Array::make(data.length(), 0.0)
  if data.is_empty() || window_size <= 0 {
    return result
  }
  let width = if window_size % 2 == 0 { window_size + 1 } else { window_size }
  let half = width / 2
  for index = 0; index < data.length(); index = index + 1 {
    let values = Array::new(capacity=width)
    let start = if index < half { 0 } else { index - half }
    let end = if index + half + 1 > data.length() {
      data.length()
    } else {
      index + half + 1
    }
    for source = start; source < end; source = source + 1 {
      values.push(data[source])
    }
    values.sort()
    result[index] = values[values.length() / 2]
  }
  result
}

///|
/// Construct a normalized discrete Gaussian kernel.
pub fn gaussian_kernel(radius : Int, sigma : Double) -> Array[Double] {
  let safe_radius = if radius < 0 { 0 } else { radius }
  let safe_sigma = if sigma <= 0.0 { 1.0 } else { sigma }
  let length = 2 * safe_radius + 1
  let kernel = Array::make(length, 0.0)
  let mut total = 0.0
  for index = 0; index < length; index = index + 1 {
    let distance = (index - safe_radius).to_double()
    let value = @math.exp(
      -0.5 * distance * distance / (safe_sigma * safe_sigma),
    )
    kernel[index] = value
    total += value
  }
  if total != 0.0 {
    for index = 0; index < length; index = index + 1 {
      kernel[index] = kernel[index] / total
    }
  }
  kernel
}

///|
/// Smooth a signal using a Gaussian kernel.
pub fn gaussian_smooth(
  data : Array[Double],
  radius : Int,
  sigma : Double,
) -> Array[Double] {
  convolve_same(data, gaussian_kernel(radius, sigma))
}

///|
/// Build a normalized Ricker wavelet for synthetic event tests.
pub fn ricker_wavelet(
  length : Int,
  dominant_frequency : Double,
  sample_rate : Double,
) -> Array[Double] {
  let result = Array::make(if length < 0 { 0 } else { length }, 0.0)
  if length <= 0 || dominant_frequency <= 0.0 || sample_rate <= 0.0 {
    return result
  }
  let center = (length - 1).to_double() / 2.0
  for index = 0; index < length; index = index + 1 {
    let time = (index.to_double() - center) / sample_rate
    let value = @math.PI * dominant_frequency * time
    let squared = value * value
    result[index] = (1.0 - 2.0 * squared) * @math.exp(-squared)
  }
  result
}
