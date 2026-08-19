///|
/// A station-local trigger with an explicit station identifier.
pub struct StationTrigger {
  station : String
  trigger : Trigger
  peak_amplitude : Double
  quality_score : Double
}

///|
/// A network-associated event hypothesis.
pub struct NetworkEvent {
  origin_time : Double
  station_count : Int
  duration : Double
  magnitude : Double
  stations : Array[String]
  confidence : Double
}

///|
/// Convert local detector results into station-tagged triggers.
pub fn station_triggers(
  station : String,
  trace : Trace,
  config : EventDetectionConfig,
) -> Array[StationTrigger] {
  let ratio = classic_sta_lta(
    trace.data,
    config.sta_samples,
    config.lta_samples,
  )
  let triggers = trigger_onset(
    ratio,
    config.onset_threshold,
    config.offset_threshold,
    config.dt,
    config.start_time,
  )
  let result = Array::new()
  let quality = assess_quality(trace, 1.0e6, 0.0000001)
  for trigger in triggers {
    let (_, peak) = event_peak(trace.data, trigger.on_index, trigger.off_index)
    result.push({
      station,
      trigger,
      peak_amplitude: peak,
      quality_score: quality.score,
    })
  }
  result
}

///|
/// Associate triggers whose onset times fall within a tolerance window.
pub fn associate_station_triggers(
  triggers : Array[StationTrigger],
  tolerance_seconds : Double,
) -> Array[NetworkEvent] {
  let result = Array::new()
  let tolerance = if tolerance_seconds < 0.0 { 0.0 } else { tolerance_seconds }
  let used = Array::make(triggers.length(), false)
  for seed = 0; seed < triggers.length(); seed = seed + 1 {
    if used[seed] {
      continue
    }
    used[seed] = true
    let origin = triggers[seed].trigger.on_time
    let stations = Array::new()
    stations.push(triggers[seed].station)
    let mut total_peak = triggers[seed].peak_amplitude
    let mut total_quality = triggers[seed].quality_score
    let mut earliest = triggers[seed].trigger.on_time
    let mut latest = triggers[seed].trigger.off_time
    for candidate = seed + 1
        candidate < triggers.length()
        candidate = candidate + 1 {
      if !used[candidate] &&
        (triggers[candidate].trigger.on_time - origin).abs() <= tolerance {
        used[candidate] = true
        stations.push(triggers[candidate].station)
        total_peak += triggers[candidate].peak_amplitude
        total_quality += triggers[candidate].quality_score
        if triggers[candidate].trigger.on_time < earliest {
          earliest = triggers[candidate].trigger.on_time
        }
        if triggers[candidate].trigger.off_time > latest {
          latest = triggers[candidate].trigger.off_time
        }
      }
    }
    let count = stations.length()
    result.push({
      origin_time: earliest,
      station_count: count,
      duration: latest - earliest,
      magnitude: if count == 0 {
        0.0
      } else {
        @math.log10(total_peak / count.to_double() + 1.0)
      },
      stations,
      confidence: if count == 0 {
        0.0
      } else {
        total_quality / count.to_double()
      },
    })
  }
  result
}

///|
/// Return only network events seen by at least a station quorum.
pub fn quorum_events(
  events : Array[NetworkEvent],
  minimum_stations : Int,
) -> Array[NetworkEvent] {
  let quorum = if minimum_stations < 1 { 1 } else { minimum_stations }
  events.filter(event => event.station_count >= quorum)
}

///|
/// Return the most confident network event.
pub fn strongest_network_event(events : Array[NetworkEvent]) -> NetworkEvent? {
  if events.is_empty() {
    return None
  }
  let mut best = events[0]
  for index = 1; index < events.length(); index = index + 1 {
    if events[index].confidence > best.confidence {
      best = events[index]
    }
  }
  Some(best)
}

///|
/// Estimate a network magnitude as the median station local magnitude.
pub fn network_local_magnitude(
  traces : Array[Trace],
  distance_km : Double,
) -> Result[Double, String] {
  let values = Array::new()
  for trace in traces {
    match local_magnitude(wood_anderson_peak(trace), distance_km) {
      Ok(value) => values.push(value)
      Err(_) => ()
    }
  }
  if values.is_empty() {
    Err("No trace produced a valid local magnitude")
  } else {
    Ok(median_samples(values))
  }
}
