///|
/// Core structure representing a single-channel seismic trace.
pub struct Trace {
  station : String
  channel : String
  sampling_rate : Double
  start_time : Double
  data : Array[Double]
}

///|
/// Create a new Trace.
pub fn Trace::new(
  station : String,
  channel : String,
  sampling_rate : Double,
  start_time : Double,
  data : Array[Double],
) -> Trace {
  { station, channel, sampling_rate, start_time, data }
}

///|
/// Get the sampling interval in seconds (delta = 1.0 / sampling_rate).
pub fn Trace::delta(self : Trace) -> Double {
  if self.sampling_rate <= 0.0 {
    0.0
  } else {
    1.0 / self.sampling_rate
  }
}

///|
/// Get the number of data points in the trace.
pub fn Trace::npts(self : Trace) -> Int {
  self.data.length()
}

///|
/// Get the duration of the trace in seconds.
pub fn Trace::duration(self : Trace) -> Double {
  let n = self.npts()
  if n <= 1 {
    0.0
  } else {
    (n - 1).to_double() * self.delta()
  }
}

///|
/// Return the time array corresponding to each sample in the trace.
pub fn Trace::times(self : Trace) -> Array[Double] {
  let n = self.npts()
  let dt = self.delta()
  let t = Array::new(capacity=n)
  for i = 0; i < n; i = i + 1 {
    t.push(self.start_time + i.to_double() * dt)
  }
  t
}

///|
/// Deep clone the Trace.
pub fn Trace::clone(self : Trace) -> Trace {
  let cloned_data = Array::new(capacity=self.data.length())
  for v in self.data {
    cloned_data.push(v)
  }
  {
    station: self.station,
    channel: self.channel,
    sampling_rate: self.sampling_rate,
    start_time: self.start_time,
    data: cloned_data,
  }
}

///|
/// Calculate the mean value of the trace data.
pub fn Trace::mean(self : Trace) -> Double {
  let n = self.npts()
  if n == 0 {
    return 0.0
  }
  let mut sum = 0.0
  for v in self.data {
    sum = sum + v
  }
  sum / n.to_double()
}

///|
/// Calculate the standard deviation of the trace data.
pub fn Trace::std(self : Trace) -> Double {
  let n = self.npts()
  if n <= 1 {
    return 0.0
  }
  let mu = self.mean()
  let mut sum_sq = 0.0
  for v in self.data {
    let diff = v - mu
    sum_sq = sum_sq + diff * diff
  }
  (sum_sq / (n - 1).to_double()).sqrt()
}

///|
/// Remove the mean value from the trace data in-place.
pub fn Trace::demean(self : Trace) -> Unit {
  let mu = self.mean()
  let n = self.npts()
  for i = 0; i < n; i = i + 1 {
    self.data[i] = self.data[i] - mu
  }
}

///|
/// Remove linear trend from the trace data in-place using least squares regression.
pub fn Trace::detrend(self : Trace) -> Unit {
  let n = self.npts()
  if n <= 1 {
    return
  }

  // Calculate means of x (sample indices) and y (sample values)
  let mut sum_y = 0.0
  for v in self.data {
    sum_y = sum_y + v
  }
  let mean_y = sum_y / n.to_double()
  let mean_x = (n - 1).to_double() / 2.0

  // Calculate covariance S_xy and variance S_xx
  let mut s_xy = 0.0
  let mut s_xx = 0.0
  for i = 0; i < n; i = i + 1 {
    let dx = i.to_double() - mean_x
    let dy = self.data[i] - mean_y
    s_xy = s_xy + dx * dy
    s_xx = s_xx + dx * dx
  }

  if s_xx == 0.0 {
    return
  }

  let slope = s_xy / s_xx
  let intercept = mean_y - slope * mean_x

  // Subtract the trend: y_i = y_i - (slope * i + intercept)
  for i = 0; i < n; i = i + 1 {
    self.data[i] = self.data[i] - (slope * i.to_double() + intercept)
  }
}
