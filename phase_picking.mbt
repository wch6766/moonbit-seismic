///|
/// Arrival picking methods exposed by the high-level API.
pub enum PickMethod {
  AmplitudeThreshold
  EnergyRatio
  AkaikeInformationCriterion
  KurtosisChange
}

///|
/// Return pick method constructors.
pub fn pick_by_amplitude() -> PickMethod {
  AmplitudeThreshold
}

///|
pub fn pick_by_energy_ratio() -> PickMethod {
  EnergyRatio
}

///|
pub fn pick_by_aic() -> PickMethod {
  AkaikeInformationCriterion
}

///|
pub fn pick_by_kurtosis() -> PickMethod {
  KurtosisChange
}

///|
/// A phase pick with method and confidence metadata.
pub struct ArrivalPick {
  index : Int
  time : Double
  strategy : PickMethod
  confidence : Double
  characteristic_value : Double
}

///|
/// Compute a running kurtosis characteristic function.
pub fn running_kurtosis(
  data : Array[Double],
  window_size : Int,
) -> Array[Double] {
  let result = Array::make(data.length(), 0.0)
  if data.is_empty() || window_size <= 1 {
    return result
  }
  for index = 0; index < data.length(); index = index + 1 {
    let start = if index + 1 < window_size {
      0
    } else {
      index + 1 - window_size
    }
    let values = slice_samples(data, start, index + 1)
    let mean = mean_samples(values)
    let deviation = standard_deviation(values)
    if deviation > 0.0 {
      let mut fourth = 0.0
      for value in values {
        let normalized = (value - mean) / deviation
        fourth += normalized * normalized * normalized * normalized
      }
      result[index] = fourth / values.length().to_double() - 3.0
    }
  }
  result
}

///|
/// Compute the AIC change-point characteristic.
pub fn aic_characteristic(
  data : Array[Double],
  minimum_window : Int,
) -> Array[Double] {
  let result = Array::make(data.length(), 0.0)
  let minimum = if minimum_window < 2 { 2 } else { minimum_window }
  if data.length() < 2 * minimum {
    return result
  }
  for split = minimum; split < data.length() - minimum; split = split + 1 {
    let left = slice_samples(data, 0, split)
    let right = slice_samples(data, split, data.length())
    let left_variance = variance_samples(left).max(1.0e-20)
    let right_variance = variance_samples(right).max(1.0e-20)
    result[split] = split.to_double() * @math.ln(left_variance) +
      (data.length() - split).to_double() * @math.ln(right_variance)
  }
  result
}

///|
/// Pick an arrival from a selected characteristic function.
pub fn pick_characteristic(
  characteristic : Array[Double],
  trace : Trace,
  minimum_value : Double,
  strategy : PickMethod,
) -> ArrivalPick? {
  if characteristic.is_empty() || trace.sampling_rate <= 0.0 {
    return None
  }
  let mut best = 0
  let mut value = characteristic[0]
  for index = 1; index < characteristic.length(); index = index + 1 {
    let candidate = if strategy is AkaikeInformationCriterion {
      -characteristic[index]
    } else {
      characteristic[index]
    }
    let current = if strategy is AkaikeInformationCriterion {
      -value
    } else {
      value
    }
    if candidate > current {
      best = index
      value = characteristic[index]
    }
  }
  let measured = if strategy is AkaikeInformationCriterion {
    -value
  } else {
    value
  }
  if measured < minimum_value {
    None
  } else {
    Some({
      index: best,
      time: sample_to_time(trace, best),
      strategy,
      confidence: if minimum_value == 0.0 {
        1.0
      } else {
        (measured / minimum_value).min(1.0)
      },
      characteristic_value: measured,
    })
  }
}

///|
/// Pick the first amplitude threshold crossing after a pre-event window.
pub fn pick_amplitude_arrival(
  trace : Trace,
  threshold : Double,
  quiet_samples : Int,
) -> ArrivalPick? {
  match first_arrival_index(trace.data, threshold, quiet_samples) {
    None => None
    Some(index) =>
      Some({
        index,
        time: sample_to_time(trace, index),
        strategy: AmplitudeThreshold,
        confidence: 1.0,
        characteristic_value: trace.data[index].abs(),
      })
  }
}

///|
/// Pick an arrival with the selected characteristic function.
pub fn pick_arrival(
  trace : Trace,
  strategy : PickMethod,
  threshold : Double,
  window_size : Int,
) -> ArrivalPick? {
  match strategy {
    AmplitudeThreshold => pick_amplitude_arrival(trace, threshold, window_size)
    EnergyRatio =>
      pick_characteristic(
        energy_ratio(trace.data, window_size, window_size * 4),
        trace,
        threshold,
        strategy,
      )
    AkaikeInformationCriterion =>
      pick_characteristic(
        aic_characteristic(trace.data, window_size),
        trace,
        threshold,
        strategy,
      )
    KurtosisChange =>
      pick_characteristic(
        running_kurtosis(trace.data, window_size),
        trace,
        threshold,
        strategy,
      )
  }
}

///|
/// Pick all methods and return the available consensus candidates.
pub fn consensus_picks(
  trace : Trace,
  threshold : Double,
  window_size : Int,
) -> Array[ArrivalPick] {
  let result = Array::new()
  for
    strategy in [
      AmplitudeThreshold,
      EnergyRatio,
      AkaikeInformationCriterion,
      KurtosisChange,
    ] {
    match pick_arrival(trace, strategy, threshold, window_size) {
      Some(pick) => result.push(pick)
      None => ()
    }
  }
  result
}
