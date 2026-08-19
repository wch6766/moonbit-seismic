///|
/// Summary of a trace's timing and sample layout.
pub struct TraceSummary {
  station : String
  channel : String
  sample_count : Int
  sample_rate : Double
  start_time : Double
  end_time : Double
  duration : Double
  mean : Double
  standard_deviation : Double
  peak : Double
}

///|
/// Create a compact summary suitable for logs and QC tables.
pub fn summarize_trace(trace : Trace) -> TraceSummary {
  {
    station: trace.station,
    channel: trace.channel,
    sample_count: trace.npts(),
    sample_rate: trace.sampling_rate,
    start_time: trace.start_time,
    end_time: if trace.npts() == 0 {
      trace.start_time
    } else {
      trace.start_time + trace.duration()
    },
    duration: trace.duration(),
    mean: mean_samples(trace.data),
    standard_deviation: standard_deviation(trace.data),
    peak: maximum_sample(trace.data).abs().max(minimum_sample(trace.data).abs()),
  }
}

///|
/// Return a trace with a replacement data array and preserved metadata.
pub fn trace_with_data(trace : Trace, data : Array[Double]) -> Trace {
  Trace::new(
    trace.station,
    trace.channel,
    trace.sampling_rate,
    trace.start_time,
    data,
  )
}

///|
/// Add two aligned traces sample by sample.
pub fn add_traces(left : Trace, right : Trace) -> Result[Trace, String] {
  if left.sampling_rate != right.sampling_rate {
    return Err("Cannot add traces with different sample rates")
  }
  if left.start_time != right.start_time {
    return Err("Cannot add traces with different start times")
  }
  if left.npts() != right.npts() {
    return Err("Cannot add traces with different lengths")
  }
  let data = Array::new(capacity=left.npts())
  for index = 0; index < left.npts(); index = index + 1 {
    data.push(left.data[index] + right.data[index])
  }
  Ok(
    Trace::new(
      left.station,
      left.channel,
      left.sampling_rate,
      left.start_time,
      data,
    ),
  )
}

///|
/// Subtract two aligned traces sample by sample.
pub fn subtract_traces(left : Trace, right : Trace) -> Result[Trace, String] {
  if left.sampling_rate != right.sampling_rate ||
    left.start_time != right.start_time {
    return Err("Traces are not time aligned")
  }
  if left.npts() != right.npts() {
    return Err("Cannot subtract traces with different lengths")
  }
  let data = Array::new(capacity=left.npts())
  for index = 0; index < left.npts(); index = index + 1 {
    data.push(left.data[index] - right.data[index])
  }
  Ok(
    Trace::new(
      left.station,
      left.channel,
      left.sampling_rate,
      left.start_time,
      data,
    ),
  )
}

///|
/// Concatenate traces when their sampling grids meet exactly.
pub fn concatenate_traces(
  first : Trace,
  second : Trace,
) -> Result[Trace, String] {
  if first.sampling_rate != second.sampling_rate {
    return Err("Cannot concatenate traces with different sample rates")
  }
  let expected = first.start_time + first.npts().to_double() * first.delta()
  if (second.start_time - expected).abs() > first.delta() * 0.001 {
    return Err("Trace timing contains a gap or overlap")
  }
  let data = first.data.copy()
  for value in second.data {
    data.push(value)
  }
  Ok(
    Trace::new(
      first.station,
      first.channel,
      first.sampling_rate,
      first.start_time,
      data,
    ),
  )
}

///|
/// Return a trace aligned to an integer sample offset.
pub fn shift_trace(trace : Trace, sample_offset : Int) -> Trace {
  let data = Array::make(trace.npts(), 0.0)
  for index = 0; index < trace.npts(); index = index + 1 {
    let source = index - sample_offset
    if source >= 0 && source < trace.npts() {
      data[index] = trace.data[source]
    }
  }
  Trace::new(
    trace.station,
    trace.channel,
    trace.sampling_rate,
    trace.start_time + sample_offset.to_double() * trace.delta(),
    data,
  )
}

///|
/// Align a trace to a reference using normalized cross-correlation.
pub fn align_to_reference(
  trace : Trace,
  reference : Trace,
  max_lag : Int,
) -> Result[Trace, String] {
  if trace.sampling_rate != reference.sampling_rate {
    return Err("Cannot align traces with different sample rates")
  }
  let result = cross_correlate(trace.data, reference.data, max_lag)
  Ok(shift_trace(trace, -result.max_lag))
}

///|
/// Return the overlapping sample interval for two traces.
pub fn common_time_range(left : Trace, right : Trace) -> (Double, Double)? {
  if left.sampling_rate != right.sampling_rate {
    return None
  }
  let left_end = left.start_time + left.duration()
  let right_end = right.start_time + right.duration()
  let start = if left.start_time > right.start_time {
    left.start_time
  } else {
    right.start_time
  }
  let end = if left_end < right_end { left_end } else { right_end }
  if end < start {
    None
  } else {
    Some((start, end))
  }
}

///|
/// Return an amplitude envelope using a quadrature approximation from a
/// centered first difference. It is stable for short traces and sufficient
/// for event-duration and amplitude measurements.
pub fn amplitude_envelope(data : Array[Double]) -> Array[Double] {
  let result = Array::make(data.length(), 0.0)
  if data.is_empty() {
    return result
  }
  for index = 0; index < data.length(); index = index + 1 {
    let previous = if index == 0 { data[index] } else { data[index - 1] }
    let next = if index + 1 >= data.length() {
      data[index]
    } else {
      data[index + 1]
    }
    let quadrature = 0.5 * (next - previous)
    result[index] = (data[index] * data[index] + quadrature * quadrature).sqrt()
  }
  result
}

///|
/// Return a trace containing an amplitude envelope.
pub fn envelope_trace(trace : Trace) -> Trace {
  trace_with_data(trace, amplitude_envelope(trace.data))
}

///|
/// Integrate a trace and preserve its physical sampling metadata.
pub fn integrate_trace(trace : Trace) -> Trace {
  trace_with_data(trace, cumulative_integral(trace.data, trace.sampling_rate))
}

///|
/// Differentiate a trace with a finite-difference operator.
pub fn differentiate_trace(trace : Trace) -> Trace {
  trace_with_data(trace, differentiate_samples(trace.data, trace.sampling_rate))
}

///|
/// Apply common preprocessing operations to a trace copy.
pub fn preprocess_trace(
  trace : Trace,
  remove_mean : Bool,
  remove_trend : Bool,
  taper_fraction : Double,
) -> Trace {
  let result = trace.clone()
  if remove_mean {
    result.demean()
  }
  if remove_trend {
    result.detrend()
  }
  cosine_taper(result.data, taper_fraction)
  result
}

///|
/// Return a trace normalized to unit peak amplitude.
pub fn normalize_trace_peak(trace : Trace) -> Trace {
  trace_with_data(trace, normalize_peak(trace.data))
}

///|
/// Return a trace normalized to zero mean and unit standard deviation.
pub fn normalize_trace_zscore(trace : Trace) -> Trace {
  trace_with_data(trace, normalize_zscore(trace.data))
}

///|
/// Return a trace with all samples clipped to a physical range.
pub fn clamp_trace(trace : Trace, minimum : Double, maximum : Double) -> Trace {
  let data = trace.data.copy()
  clamp_samples(data, minimum, maximum)
  trace_with_data(trace, data)
}

///|
/// Compute a signal-to-noise ratio using two time windows.
pub fn trace_snr_db(
  trace : Trace,
  noise_start : Double,
  noise_end : Double,
  signal_start : Double,
  signal_end : Double,
) -> Double {
  let noise = slice_samples(
    trace.data,
    ((noise_start - trace.start_time) * trace.sampling_rate).to_int(),
    ((noise_end - trace.start_time) * trace.sampling_rate).to_int(),
  )
  let signal = slice_samples(
    trace.data,
    ((signal_start - trace.start_time) * trace.sampling_rate).to_int(),
    ((signal_end - trace.start_time) * trace.sampling_rate).to_int(),
  )
  let noise_rms = rms_samples(noise)
  let signal_rms = rms_samples(signal)
  if noise_rms <= 0.0 {
    0.0
  } else {
    20.0 * @math.log10(signal_rms / noise_rms)
  }
}

///|
/// Return true when trace metadata and samples satisfy basic invariants.
pub fn trace_is_well_formed(trace : Trace) -> Bool {
  validate_trace(trace) is Ok(_) && trace.npts() >= 0
}
