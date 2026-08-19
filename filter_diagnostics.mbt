///|
/// Measured magnitude response of a filter at one frequency.
pub struct FilterResponsePoint {
  frequency : Double
  input_rms : Double
  output_rms : Double
  gain : Double
  gain_db : Double
}

///|
/// Measure an arbitrary trace transformation over a frequency grid.
pub fn measure_filter_response(
  transform : (Trace) -> Trace,
  sample_rate : Double,
  frequencies : Array[Double],
  cycles : Int,
) -> Array[FilterResponsePoint] {
  let result = Array::new()
  let count = if cycles < 1 { 1 } else { cycles }
  for frequency in frequencies {
    let length = (sample_rate / frequency * count.to_double()).to_int().max(16)
    let input = Trace::new(
      "response",
      "test",
      sample_rate,
      0.0,
      tapered_sine(length, frequency, sample_rate, 1.0, 0.05),
    )
    let output = transform(input)
    let input_rms = rms_samples(input.data)
    let output_rms = rms_samples(output.data)
    let gain = if input_rms == 0.0 { 0.0 } else { output_rms / input_rms }
    result.push({
      frequency,
      input_rms,
      output_rms,
      gain,
      gain_db: if gain <= 0.0 {
        -300.0
      } else {
        20.0 * @math.log10(gain)
      },
    })
  }
  result
}

///|
/// Measure the existing low-pass implementation on a frequency grid.
pub fn lowpass_response(
  sample_rate : Double,
  cutoff : Double,
  frequencies : Array[Double],
) -> Array[FilterResponsePoint] {
  measure_filter_response(
    trace => {
      lowpass(trace, cutoff, false)
      trace
    },
    sample_rate,
    frequencies,
    8,
  )
}

///|
/// Find the first frequency where measured gain falls below a dB threshold.
pub fn cutoff_from_response(
  response : Array[FilterResponsePoint],
  threshold_db : Double,
) -> Double? {
  for point in response {
    if point.gain_db <= threshold_db {
      return Some(point.frequency)
    }
  }
  None
}

///|
/// Return the maximum passband ripple in dB.
pub fn passband_ripple(
  response : Array[FilterResponsePoint],
  low_frequency : Double,
  high_frequency : Double,
) -> Double {
  let values = Array::new()
  for point in response {
    if point.frequency >= low_frequency && point.frequency <= high_frequency {
      values.push(point.gain_db)
    }
  }
  if values.is_empty() {
    0.0
  } else {
    maximum_sample(values) - minimum_sample(values)
  }
}
