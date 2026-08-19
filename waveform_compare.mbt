///|
/// A comparable error report for two traces on a shared grid.
pub struct TraceComparison {
  sample_count : Int
  mean_error : Double
  rms_error : Double
  maximum_error : Double
  correlation : Double
  lag : Int
}

///|
/// Compare two traces after checking their sample rates.
pub fn compare_traces(
  left : Trace,
  right : Trace,
  maximum_lag : Int,
) -> Result[TraceComparison, String] {
  if left.sampling_rate != right.sampling_rate {
    return Err("Trace comparison requires equal sample rates")
  }
  let count = if left.npts() < right.npts() {
    left.npts()
  } else {
    right.npts()
  }
  if count == 0 {
    return Ok({
      sample_count: 0,
      mean_error: 0.0,
      rms_error: 0.0,
      maximum_error: 0.0,
      correlation: 0.0,
      lag: 0,
    })
  }
  let errors = Array::new(capacity=count)
  let mut maximum = 0.0
  for index = 0; index < count; index = index + 1 {
    let error = left.data[index] - right.data[index]
    errors.push(error)
    if error.abs() > maximum {
      maximum = error.abs()
    }
  }
  let correlation = pearson_correlation(
    slice_samples(left.data, 0, count),
    slice_samples(right.data, 0, count),
  )
  let correlation_result = cross_correlate(
    slice_samples(left.data, 0, count),
    slice_samples(right.data, 0, count),
    maximum_lag,
  )
  Ok({
    sample_count: count,
    mean_error: mean_samples(errors),
    rms_error: rms_samples(errors),
    maximum_error: maximum,
    correlation,
    lag: correlation_result.max_lag,
  })
}

///|
/// Return true when a comparison meets correlation and RMS tolerances.
pub fn comparison_passes(
  comparison : TraceComparison,
  minimum_correlation : Double,
  maximum_rms_error : Double,
) -> Bool {
  comparison.correlation >= minimum_correlation &&
  comparison.rms_error <= maximum_rms_error
}

///|
/// Return the normalized improvement between two RMS errors.
pub fn rms_improvement(before : Double, after : Double) -> Double {
  if before <= 0.0 {
    0.0
  } else {
    (before - after) / before
  }
}

///|
/// Return an error summary as CSV.
pub fn comparison_csv(comparison : TraceComparison) -> String {
  "samples,mean_error,rms_error,max_error,correlation,lag\n" +
  comparison.sample_count.to_string() +
  "," +
  comparison.mean_error.to_string() +
  "," +
  comparison.rms_error.to_string() +
  "," +
  comparison.maximum_error.to_string() +
  "," +
  comparison.correlation.to_string() +
  "," +
  comparison.lag.to_string() +
  "\n"
}

///|
/// Return the percentage of samples within an absolute error tolerance.
pub fn agreement_ratio(
  left : Array[Double],
  right : Array[Double],
  tolerance : Double,
) -> Double {
  let count = if left.length() < right.length() {
    left.length()
  } else {
    right.length()
  }
  if count == 0 {
    return 0.0
  }
  let mut accepted = 0
  for index = 0; index < count; index = index + 1 {
    if (left[index] - right[index]).abs() <= tolerance.abs() {
      accepted += 1
    }
  }
  accepted.to_double() / count.to_double()
}

///|
/// Return a trace copy after removing its mean error against a reference.
pub fn bias_correct(trace : Trace, reference : Trace) -> Trace {
  let count = if trace.npts() < reference.npts() {
    trace.npts()
  } else {
    reference.npts()
  }
  let errors = Array::new(capacity=count)
  for index = 0; index < count; index = index + 1 {
    errors.push(trace.data[index] - reference.data[index])
  }
  let bias = mean_samples(errors)
  let data = trace.data.copy()
  for index = 0; index < data.length(); index = index + 1 {
    data[index] = data[index] - bias
  }
  trace_with_data(trace, data)
}

///|
/// Return the maximum absolute difference between two sample arrays.
pub fn maximum_difference(
  left : Array[Double],
  right : Array[Double],
) -> Double {
  let count = if left.length() < right.length() {
    left.length()
  } else {
    right.length()
  }
  let mut maximum = 0.0
  for index = 0; index < count; index = index + 1 {
    let difference = (left[index] - right[index]).abs()
    if difference > maximum {
      maximum = difference
    }
  }
  maximum
}

///|
/// Return a copy of a signal with its initial sample removed as DC offset.
pub fn remove_initial_offset(
  data : Array[Double],
  samples : Int,
) -> Array[Double] {
  let count = if samples < 1 {
    1
  } else if samples > data.length() {
    data.length()
  } else {
    samples
  }
  let offset = mean_samples(slice_samples(data, 0, count))
  let result = data.copy()
  for index = 0; index < result.length(); index = index + 1 {
    result[index] = result[index] - offset
  }
  result
}

///|
/// Return a normalized absolute difference score in `[0, 1]`.
pub fn difference_score(left : Array[Double], right : Array[Double]) -> Double {
  let scale = rms_samples(left) + rms_samples(right) + 0.0000000001
  let count = if left.length() < right.length() {
    left.length()
  } else {
    right.length()
  }
  if count == 0 {
    0.0
  } else {
    1.0 - (maximum_difference(left, right) / scale).min(1.0)
  }
}
