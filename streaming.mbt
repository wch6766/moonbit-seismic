///|
/// Online mean and variance accumulator using Welford's stable update.
pub struct RollingAccumulator {
  mut count : Int
  mut mean : Double
  mut m2 : Double
}

///|
/// Construct an empty accumulator.
pub fn RollingAccumulator::new() -> RollingAccumulator {
  { count: 0, mean: 0.0, m2: 0.0 }
}

///|
/// Add one sample to an accumulator.
pub fn RollingAccumulator::push(
  self : RollingAccumulator,
  value : Double,
) -> Unit {
  self.count = self.count + 1
  let difference = value - self.mean
  self.mean = self.mean + difference / self.count.to_double()
  let second_difference = value - self.mean
  self.m2 = self.m2 + difference * second_difference
}

///|
/// Return the number of samples seen.
pub fn RollingAccumulator::count(self : RollingAccumulator) -> Int {
  self.count
}

///|
/// Return the current mean.
pub fn RollingAccumulator::mean(self : RollingAccumulator) -> Double {
  self.mean
}

///|
/// Return the population variance.
pub fn RollingAccumulator::variance(self : RollingAccumulator) -> Double {
  if self.count == 0 {
    0.0
  } else {
    self.m2 / self.count.to_double()
  }
}

///|
/// Return the current standard deviation.
pub fn RollingAccumulator::standard_deviation(
  self : RollingAccumulator,
) -> Double {
  self.variance().sqrt()
}

///|
/// Reset an accumulator without allocating a new object.
pub fn RollingAccumulator::reset(self : RollingAccumulator) -> Unit {
  self.count = 0
  self.mean = 0.0
  self.m2 = 0.0
}

///|
/// Streaming STA/LTA detector state.
pub struct StreamingStaLta {
  sta_samples : Int
  lta_samples : Int
  sta_buffer : Array[Double]
  lta_buffer : Array[Double]
  mut sta_sum : Double
  mut lta_sum : Double
}

///|
/// Construct a streaming STA/LTA detector.
pub fn StreamingStaLta::new(
  sta_samples : Int,
  lta_samples : Int,
) -> StreamingStaLta {
  {
    sta_samples: if sta_samples < 1 {
      1
    } else {
      sta_samples
    },
    lta_samples: if lta_samples < 1 {
      1
    } else {
      lta_samples
    },
    sta_buffer: [],
    lta_buffer: [],
    sta_sum: 0.0,
    lta_sum: 0.0,
  }
}

///|
/// Push an amplitude and return the current recursive ratio.
pub fn StreamingStaLta::push(self : StreamingStaLta, value : Double) -> Double {
  let amplitude = value.abs()
  self.sta_buffer.push(amplitude)
  self.sta_sum = self.sta_sum + amplitude
  if self.sta_buffer.length() > self.sta_samples {
    self.sta_sum = self.sta_sum - self.sta_buffer.remove(0)
  }
  self.lta_buffer.push(amplitude)
  self.lta_sum = self.lta_sum + amplitude
  if self.lta_buffer.length() > self.lta_samples {
    self.lta_sum = self.lta_sum - self.lta_buffer.remove(0)
  }
  let sta = self.sta_sum / self.sta_buffer.length().to_double()
  let lta = self.lta_sum / self.lta_buffer.length().to_double()
  if lta == 0.0 {
    1.0
  } else {
    sta / lta
  }
}

///|
/// Clear streaming detector state.
pub fn StreamingStaLta::reset(self : StreamingStaLta) -> Unit {
  self.sta_buffer.clear()
  self.lta_buffer.clear()
  self.sta_sum = 0.0
  self.lta_sum = 0.0
}

///|
/// A fixed-size ring buffer for streaming waveform samples.
pub struct SampleRing {
  capacity : Int
  values : Array[Double]
}

///|
/// Construct a ring buffer.
pub fn SampleRing::new(capacity : Int) -> SampleRing {
  { capacity: if capacity < 1 { 1 } else { capacity }, values: [] }
}

///|
/// Add a sample, evicting the oldest value when full.
pub fn SampleRing::push(self : SampleRing, value : Double) -> Unit {
  self.values.push(value)
  if self.values.length() > self.capacity {
    ignore(self.values.remove(0))
  }
}

///|
/// Return a snapshot of the current ring content.
pub fn SampleRing::snapshot(self : SampleRing) -> Array[Double] {
  self.values.copy()
}

///|
/// Return whether a ring has reached its configured capacity.
pub fn SampleRing::is_full(self : SampleRing) -> Bool {
  self.values.length() >= self.capacity
}

///|
/// A streaming quality monitor with configurable batch size.
pub struct StreamingQualityMonitor {
  batch_size : Int
  mut batch : Array[Double]
  reports : Array[QualityReport]
}

///|
/// Construct a streaming quality monitor.
pub fn StreamingQualityMonitor::new(
  batch_size : Int,
) -> StreamingQualityMonitor {
  {
    batch_size: if batch_size < 1 {
      1
    } else {
      batch_size
    },
    batch: [],
    reports: [],
  }
}

///|
/// Add a sample and optionally emit a report when a batch is complete.
pub fn StreamingQualityMonitor::push(
  self : StreamingQualityMonitor,
  value : Double,
  sample_rate : Double,
) -> QualityReport? {
  self.batch.push(value)
  if self.batch.length() < self.batch_size {
    return None
  }
  let trace = Trace::new("stream", "data", sample_rate, 0.0, self.batch)
  let report = assess_quality(trace, 1.0e6, 0.0000001)
  self.reports.push(report)
  self.batch = []
  Some(report)
}

///|
/// Return all completed quality reports.
pub fn StreamingQualityMonitor::reports(
  self : StreamingQualityMonitor,
) -> Array[QualityReport] {
  self.reports.copy()
}
