///|
/// Configuration for production-style event detection.
pub struct EventDetectionConfig {
  sta_samples : Int
  lta_samples : Int
  onset_threshold : Double
  offset_threshold : Double
  minimum_samples : Int
  merge_gap_samples : Int
  dt : Double
  start_time : Double
}

///|
/// A detected event with timing and signal quality measurements.
pub struct EventCandidate {
  on_index : Int
  off_index : Int
  peak_index : Int
  on_time : Double
  off_time : Double
  peak_time : Double
  peak_ratio : Double
  energy : Double
  signal_to_noise : Double
}

///|
/// Return a conservative default detector configuration.
pub fn default_event_config(
  dt : Double,
  start_time : Double,
) -> EventDetectionConfig {
  {
    sta_samples: 5,
    lta_samples: 20,
    onset_threshold: 2.0,
    offset_threshold: 1.2,
    minimum_samples: 3,
    merge_gap_samples: 2,
    dt,
    start_time,
  }
}

///|
/// Compute the short-to-long energy ratio with an O(N) prefix sum.
pub fn energy_ratio(
  data : Array[Double],
  short_samples : Int,
  long_samples : Int,
) -> Array[Double] {
  let result = Array::make(data.length(), 1.0)
  if data.is_empty() || short_samples <= 0 || long_samples <= short_samples {
    return result
  }
  let prefix = Array::make(data.length() + 1, 0.0)
  for index = 0; index < data.length(); index = index + 1 {
    prefix[index + 1] = prefix[index] + data[index] * data[index]
  }
  for index = long_samples - 1; index < data.length(); index = index + 1 {
    let long_start = index + 1 - long_samples
    let short_start = index + 1 - short_samples
    let long_energy = prefix[index + 1] - prefix[long_start]
    let short_energy = prefix[index + 1] - prefix[short_start]
    let long_mean = long_energy / long_samples.to_double()
    let short_mean = short_energy / short_samples.to_double()
    result[index] = if long_mean <= 0.0 { 1.0 } else { short_mean / long_mean }
  }
  result
}

///|
/// Compute a recursive STA/LTA ratio for streaming-friendly detection.
pub fn recursive_sta_lta(
  data : Array[Double],
  sta_samples : Int,
  lta_samples : Int,
) -> Array[Double] {
  let result = Array::make(data.length(), 1.0)
  if data.is_empty() || sta_samples <= 0 || lta_samples <= sta_samples {
    return result
  }
  let sta_decay = @math.exp(-1.0 / sta_samples.to_double())
  let lta_decay = @math.exp(-1.0 / lta_samples.to_double())
  let mut sta = 0.0
  let mut lta = 0.0
  for index = 0; index < data.length(); index = index + 1 {
    let amplitude = data[index].abs()
    sta = sta_decay * sta + (1.0 - sta_decay) * amplitude
    lta = lta_decay * lta + (1.0 - lta_decay) * amplitude
    result[index] = if lta <= 0.0 { 1.0 } else { sta / lta }
  }
  result
}

///|
/// Compute a short-time characteristic function from first differences.
pub fn difference_characteristic(data : Array[Double]) -> Array[Double] {
  let result = Array::make(data.length(), 0.0)
  if data.length() < 2 {
    return result
  }
  for index = 1; index < data.length(); index = index + 1 {
    result[index] = (data[index] - data[index - 1]).abs()
  }
  result
}

///|
/// Return the first index crossing a threshold after a quiet pre-event window.
pub fn first_arrival_index(
  data : Array[Double],
  threshold : Double,
  quiet_samples : Int,
) -> Int? {
  if data.is_empty() {
    return None
  }
  let quiet = if quiet_samples < 1 { 1 } else { quiet_samples }
  let baseline_data = slice_samples(data, 0, quiet)
  let baseline = mean_samples(baseline_data) +
    3.0 * standard_deviation(baseline_data)
  let effective = if threshold > baseline { threshold } else { baseline }
  for index = quiet; index < data.length(); index = index + 1 {
    if data[index].abs() >= effective {
      return Some(index)
    }
  }
  None
}

///|
/// Find the peak index and value within an inclusive event range.
fn event_peak(data : Array[Double], start : Int, end : Int) -> (Int, Double) {
  let first = if start < 0 { 0 } else { start }
  let last = if end >= data.length() { data.length() - 1 } else { end }
  if data.is_empty() || first > last {
    return (0, 0.0)
  }
  let mut peak_index = first
  let mut peak = data[first].abs()
  for index = first + 1; index <= last; index = index + 1 {
    if data[index].abs() > peak {
      peak = data[index].abs()
      peak_index = index
    }
  }
  (peak_index, peak)
}

///|
/// Merge adjacent trigger windows separated by at most `gap_samples`.
pub fn merge_triggers(
  triggers : Array[Trigger],
  gap_samples : Int,
) -> Array[Trigger] {
  let result = Array::new()
  if triggers.is_empty() {
    return result
  }
  let gap = if gap_samples < 0 { 0 } else { gap_samples }
  let mut current = triggers[0]
  for index = 1; index < triggers.length(); index = index + 1 {
    let next = triggers[index]
    if next.on_index - current.off_index <= gap + 1 {
      current = {
        on_index: current.on_index,
        off_index: if next.off_index > current.off_index {
          next.off_index
        } else {
          current.off_index
        },
        on_time: current.on_time,
        off_time: if next.off_time > current.off_time {
          next.off_time
        } else {
          current.off_time
        },
      }
    } else {
      result.push(current)
      current = next
    }
  }
  result.push(current)
  result
}

///|
/// Convert a trigger into a feature-rich event candidate.
fn candidate_from_trigger(
  trigger : Trigger,
  data : Array[Double],
  ratio : Array[Double],
  dt : Double,
) -> EventCandidate {
  let (peak_index, peak) = event_peak(data, trigger.on_index, trigger.off_index)
  let values = slice_samples(data, trigger.on_index, trigger.off_index + 1)
  let noise_end = if trigger.on_index < 20 { trigger.on_index } else { 20 }
  let noise = slice_samples(data, 0, noise_end)
  let noise_rms = rms_samples(noise)
  let snr = if noise_rms <= 0.0 {
    0.0
  } else {
    20.0 * @math.log10(peak / noise_rms)
  }
  let peak_ratio = if peak_index < ratio.length() {
    ratio[peak_index]
  } else {
    0.0
  }
  {
    on_index: trigger.on_index,
    off_index: trigger.off_index,
    peak_index,
    on_time: trigger.on_time,
    off_time: trigger.off_time,
    peak_time: trigger.on_time +
    (peak_index - trigger.on_index).to_double() * dt,
    peak_ratio,
    energy: signal_energy(values),
    signal_to_noise: snr,
  }
}

///|
/// Detect, merge, and characterize event windows.
pub fn detect_events(
  data : Array[Double],
  config : EventDetectionConfig,
) -> Array[EventCandidate] {
  let ratio = classic_sta_lta(data, config.sta_samples, config.lta_samples)
  let raw = trigger_onset(
    ratio,
    config.onset_threshold,
    config.offset_threshold,
    config.dt,
    config.start_time,
  )
  let merged = merge_triggers(raw, config.merge_gap_samples)
  let result = Array::new()
  for trigger in merged {
    if trigger.off_index - trigger.on_index + 1 >= config.minimum_samples {
      result.push(candidate_from_trigger(trigger, data, ratio, config.dt))
    }
  }
  result
}

///|
/// Estimate the energy in a physical time interval.
pub fn event_energy(
  trace : Trace,
  start_time : Double,
  end_time : Double,
) -> Double {
  if trace.sampling_rate <= 0.0 || end_time <= start_time {
    return 0.0
  }
  let first = ((start_time - trace.start_time) * trace.sampling_rate).to_int()
  let last = ((end_time - trace.start_time) * trace.sampling_rate).to_int() + 1
  signal_energy(slice_samples(trace.data, first, last)) / trace.sampling_rate
}

///|
/// Return events whose signal-to-noise estimate passes a minimum in dB.
pub fn filter_events_by_snr(
  events : Array[EventCandidate],
  minimum_db : Double,
) -> Array[EventCandidate] {
  events.filter(event => event.signal_to_noise >= minimum_db)
}

///|
/// Return the largest candidate by energy.
pub fn strongest_event(events : Array[EventCandidate]) -> EventCandidate? {
  if events.is_empty() {
    return None
  }
  let mut best = events[0]
  for index = 1; index < events.length(); index = index + 1 {
    if events[index].energy > best.energy {
      best = events[index]
    }
  }
  Some(best)
}

///|
/// Compute the duration of an event in seconds.
pub fn event_duration(event : EventCandidate) -> Double {
  event.off_time - event.on_time
}
