///| Constant value of PI for frequency calculations.
const PI : Double = 3.14159265358979323846

///| Helper to reverse array elements in place.
fn reverse_data(data : Array[Double]) -> Unit {
  let n = data.length()
  for i = 0; i < n / 2; i = i + 1 {
    let temp = data[i]
    data[i] = data[n - 1 - i]
    data[n - 1 - i] = temp
  }
}

///| Apply a 2nd-order Biquad filter forward in-place.
pub fn apply_biquad(
  trace : Trace,
  b0 : Double,
  b1 : Double,
  b2 : Double,
  a1 : Double,
  a2 : Double
) -> Unit {
  let n = trace.npts()
  if n == 0 {
    return
  }
  
  let mut x1 = 0.0
  let mut x2 = 0.0
  let mut y1 = 0.0
  let mut y2 = 0.0
  
  for i = 0; i < n; i = i + 1 {
    let xi = trace.data[i]
    let yi = b0 * xi + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2
    
    x2 = x1
    x1 = xi
    y2 = y1
    y1 = yi
    
    trace.data[i] = yi
  }
}

///| Apply a 2nd-order Biquad filter forward and backward to achieve zero phase shift.
pub fn apply_biquad_zero_phase(
  trace : Trace,
  b0 : Double,
  b1 : Double,
  b2 : Double,
  a1 : Double,
  a2 : Double
) -> Unit {
  apply_biquad(trace, b0, b1, b2, a1, a2)
  reverse_data(trace.data)
  apply_biquad(trace, b0, b1, b2, a1, a2)
  reverse_data(trace.data)
}

///| Apply a 2nd-order Butterworth low-pass filter to the trace.
pub fn lowpass(trace : Trace, freq : Double, zero_phase : Bool) -> Unit {
  let fs = trace.sampling_rate
  if freq >= fs / 2.0 || freq <= 0.0 {
    return
  }
  
  let omega0 = 2.0 * PI * freq / fs
  let k = @math.tan(omega0 / 2.0)
  let k_sq = k * k
  let sqrt2 = 1.4142135623730951
  
  let d = 1.0 + sqrt2 * k + k_sq
  let b0 = k_sq / d
  let b1 = 2.0 * b0
  let b2 = b0
  let a1 = 2.0 * (k_sq - 1.0) / d
  let a2 = (1.0 - sqrt2 * k + k_sq) / d
  
  if zero_phase {
    apply_biquad_zero_phase(trace, b0, b1, b2, a1, a2)
  } else {
    apply_biquad(trace, b0, b1, b2, a1, a2)
  }
}

///| Apply a 2nd-order Butterworth high-pass filter to the trace.
pub fn highpass(trace : Trace, freq : Double, zero_phase : Bool) -> Unit {
  let fs = trace.sampling_rate
  if freq >= fs / 2.0 || freq <= 0.0 {
    return
  }
  
  let omega0 = 2.0 * PI * freq / fs
  let k = @math.tan(omega0 / 2.0)
  let k_sq = k * k
  let sqrt2 = 1.4142135623730951
  
  let d = 1.0 + sqrt2 * k + k_sq
  let b0 = 1.0 / d
  let b1 = -2.0 * b0
  let b2 = b0
  let a1 = 2.0 * (k_sq - 1.0) / d
  let a2 = (1.0 - sqrt2 * k + k_sq) / d
  
  if zero_phase {
    apply_biquad_zero_phase(trace, b0, b1, b2, a1, a2)
  } else {
    apply_biquad(trace, b0, b1, b2, a1, a2)
  }
}

///| Apply a 2nd-order Butterworth band-pass filter to the trace.
pub fn bandpass(trace : Trace, f_min : Double, f_max : Double, zero_phase : Bool) -> Unit {
  let fs = trace.sampling_rate
  if f_min <= 0.0 || f_max >= fs / 2.0 || f_min >= f_max {
    return
  }
  
  let f0 = (f_min * f_max).sqrt()
  let bw = f_max - f_min
  let q = f0 / bw
  let omega0 = 2.0 * PI * f0 / fs
  let alpha = @math.sin(omega0) / (2.0 * q)
  
  let d = 1.0 + alpha
  let b0 = alpha / d
  let b1 = 0.0
  let b2 = -alpha / d
  let a1 = -2.0 * @math.cos(omega0) / d
  let a2 = (1.0 - alpha) / d
  
  if zero_phase {
    apply_biquad_zero_phase(trace, b0, b1, b2, a1, a2)
  } else {
    apply_biquad(trace, b0, b1, b2, a1, a2)
  }
}

///| Apply a 2nd-order Butterworth band-stop (notch) filter to the trace.
pub fn bandstop(trace : Trace, f_min : Double, f_max : Double, zero_phase : Bool) -> Unit {
  let fs = trace.sampling_rate
  if f_min <= 0.0 || f_max >= fs / 2.0 || f_min >= f_max {
    return
  }
  
  let f0 = (f_min * f_max).sqrt()
  let bw = f_max - f_min
  let q = f0 / bw
  let omega0 = 2.0 * PI * f0 / fs
  let alpha = @math.sin(omega0) / (2.0 * q)
  
  let d = 1.0 + alpha
  let b0 = 1.0 / d
  let b1 = -2.0 * @math.cos(omega0) / d
  let b2 = 1.0 / d
  let a1 = -2.0 * @math.cos(omega0) / d
  let a2 = (1.0 - alpha) / d
  
  if zero_phase {
    apply_biquad_zero_phase(trace, b0, b1, b2, a1, a2)
  } else {
    apply_biquad(trace, b0, b1, b2, a1, a2)
  }
}

///| Apply a simple Moving Average smoothing filter to the trace.
pub fn moving_average(trace : Trace, window_size : Int) -> Unit {
  let n = trace.npts()
  if n == 0 || window_size <= 1 {
    return
  }
  
  let half = window_size / 2
  let temp = Array::new(capacity=n)
  for v in trace.data {
    temp.push(v)
  }
  
  for i = 0; i < n; i = i + 1 {
    let mut sum = 0.0
    let mut count = 0
    let start = if i - half < 0 { 0 } else { i - half }
    let end = if i + half >= n { n - 1 } else { i + half }
    for j = start; j <= end; j = j + 1 {
      sum = sum + temp[j]
      count = count + 1
    }
    trace.data[i] = sum / count.to_double()
  }
}
