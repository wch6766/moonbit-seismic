///|
/// Wave phase selector for travel-time and ray diagnostics.
pub enum WavePhase {
  PWave
  SWave
}

///|
/// Return the P-wave phase tag.
pub fn wave_phase_p() -> WavePhase {
  PWave
}

///|
/// Return the S-wave phase tag.
pub fn wave_phase_s() -> WavePhase {
  SWave
}

///|
/// A single sample in a travel-time curve.
pub struct TravelTimeSample {
  distance : Double
  source_depth : Double
  time : Double
  phase : WavePhase
}

///|
/// A sampled travel-time curve for plotting or phase picking.
pub struct TravelTimeCurve {
  samples : Array[TravelTimeSample]
  phase : WavePhase
  minimum_time : Double
  maximum_time : Double
}

///|
/// Validate a one-dimensional velocity model.
pub fn validate_velocity_model(model : VelocityModel) -> Result[Unit, String] {
  if model.layers.is_empty() {
    return Err("Velocity model must contain at least one layer")
  }
  for index = 0; index < model.layers.length(); index = index + 1 {
    let layer = model.layers[index]
    if layer.depth < 0.0 ||
      layer.vp <= 0.0 ||
      layer.vs <= 0.0 ||
      layer.vs > layer.vp {
      return Err("Velocity layers require depth >= 0 and 0 < Vs <= Vp")
    }
    if index > 0 && layer.depth <= model.layers[index - 1].depth {
      return Err("Velocity layer depths must be strictly increasing")
    }
  }
  Ok(())
}

///|
/// Return the shallowest and deepest defined layer depths.
pub fn velocity_depth_range(model : VelocityModel) -> (Double, Double)? {
  if model.layers.is_empty() {
    return None
  }
  Some((model.layers[0].depth, model.layers[model.layers.length() - 1].depth))
}

///|
/// Return the minimum velocity over all model layers for a phase.
pub fn minimum_model_velocity(
  model : VelocityModel,
  phase : WavePhase,
) -> Double {
  if model.layers.is_empty() {
    return 0.0
  }
  let mut minimum = 1.0e300
  for layer in model.layers {
    let value = match phase {
      PWave => layer.vp
      SWave => layer.vs
    }
    if value < minimum {
      minimum = value
    }
  }
  minimum
}

///|
/// Return the maximum velocity over all model layers for a phase.
pub fn maximum_model_velocity(
  model : VelocityModel,
  phase : WavePhase,
) -> Double {
  let mut maximum = 0.0
  for layer in model.layers {
    let value = match phase {
      PWave => layer.vp
      SWave => layer.vs
    }
    if value > maximum {
      maximum = value
    }
  }
  maximum
}

///|
/// Return a phase-aware wrapper around the existing travel-time API.
pub fn travel_time_phase(
  model : VelocityModel,
  source_depth : Double,
  distance : Double,
  phase : WavePhase,
) -> Double {
  match phase {
    PWave => travel_time(model, source_depth, distance, false)
    SWave => travel_time(model, source_depth, distance, true)
  }
}

///|
/// Integrate slowness along a straight ray with a selectable step count.
pub fn travel_time_steps(
  model : VelocityModel,
  source_depth : Double,
  distance : Double,
  phase : WavePhase,
  steps : Int,
) -> Double {
  if source_depth < 0.0 || distance < 0.0 {
    return 0.0
  }
  let count = if steps < 1 { 1 } else { steps }
  let path_length = (source_depth * source_depth + distance * distance).sqrt()
  if path_length == 0.0 {
    return 0.0
  }
  let mut slowness = 0.0
  for index = 0; index <= count; index = index + 1 {
    let fraction = index.to_double() / count.to_double()
    let depth = source_depth * (1.0 - fraction)
    let velocity = get_velocity_at_depth(model, depth, phase is SWave)
    let weight = if index == 0 || index == count { 0.5 } else { 1.0 }
    slowness += weight / velocity
  }
  path_length * slowness / count.to_double()
}

///|
/// Return uniformly spaced points along a straight source-receiver ray.
pub fn ray_path_points(
  source_depth : Double,
  distance : Double,
  points : Int,
) -> Array[(Double, Double)] {
  let result = Array::new()
  let count = if points < 2 { 2 } else { points }
  for index = 0; index < count; index = index + 1 {
    let fraction = index.to_double() / (count - 1).to_double()
    result.push((distance * fraction, source_depth * (1.0 - fraction)))
  }
  result
}

///|
/// Compute a depth-dependent vertical travel time by layer integration.
pub fn vertical_travel_time(
  model : VelocityModel,
  source_depth : Double,
  phase : WavePhase,
) -> Double {
  if source_depth <= 0.0 || model.layers.is_empty() {
    return 0.0
  }
  let mut time = 0.0
  let mut top = 0.0
  for index = 0; index < model.layers.length(); index = index + 1 {
    let bottom = model.layers[index].depth
    let end = if source_depth < bottom { source_depth } else { bottom }
    if end > top {
      let thickness = end - top
      let velocity = match phase {
        PWave => model.layers[index].vp
        SWave => model.layers[index].vs
      }
      time += thickness / velocity
    }
    top = bottom
    if top >= source_depth {
      break
    }
  }
  if top < source_depth {
    let velocity = match phase {
      PWave => model.layers[model.layers.length() - 1].vp
      SWave => model.layers[model.layers.length() - 1].vs
    }
    time += (source_depth - top) / velocity
  }
  time
}

///|
/// Estimate a ray parameter from a straight-ray angle.
pub fn estimate_ray_parameter(
  model : VelocityModel,
  source_depth : Double,
  distance : Double,
  phase : WavePhase,
) -> Double {
  let velocity = get_velocity_at_depth(
    model,
    source_depth / 2.0,
    phase is SWave,
  )
  let angle = @math.atan2(distance, source_depth)
  if velocity <= 0.0 {
    0.0
  } else {
    @math.sin(angle) / velocity
  }
}

///|
/// Sample a source-depth/distance grid for a travel-time plot.
pub fn sample_travel_time_curve(
  model : VelocityModel,
  source_depth : Double,
  minimum_distance : Double,
  maximum_distance : Double,
  samples : Int,
  phase : WavePhase,
) -> TravelTimeCurve {
  let result = Array::new()
  let count = if samples < 1 { 1 } else { samples }
  for index = 0; index < count; index = index + 1 {
    let fraction = if count == 1 {
      0.0
    } else {
      index.to_double() / (count - 1).to_double()
    }
    let distance = minimum_distance +
      fraction * (maximum_distance - minimum_distance)
    result.push({
      distance,
      source_depth,
      time: travel_time_steps(model, source_depth, distance, phase, 128),
      phase,
    })
  }
  let times = Array::new(capacity=result.length())
  for sample in result {
    times.push(sample.time)
  }
  {
    samples: result,
    phase,
    minimum_time: minimum_sample(times),
    maximum_time: maximum_sample(times),
  }
}

///|
/// Return a P and S arrival pair for a source-receiver geometry.
pub fn phase_arrivals(
  model : VelocityModel,
  source_depth : Double,
  distance : Double,
) -> Array[TravelTimeSample] {
  [
    {
      distance,
      source_depth,
      time: travel_time_phase(model, source_depth, distance, PWave),
      phase: PWave,
    },
    {
      distance,
      source_depth,
      time: travel_time_phase(model, source_depth, distance, SWave),
      phase: SWave,
    },
  ]
}

///|
/// Return the P-minus-S or S-minus-P arrival separation in seconds.
pub fn phase_separation(
  model : VelocityModel,
  source_depth : Double,
  distance : Double,
) -> Double {
  travel_time_phase(model, source_depth, distance, SWave) -
  travel_time_phase(model, source_depth, distance, PWave)
}

///|
/// Estimate an epicentral distance from a P/S time separation by bisection.
pub fn distance_from_phase_separation(
  model : VelocityModel,
  source_depth : Double,
  observed_separation : Double,
  maximum_distance : Double,
) -> Double {
  if observed_separation <= 0.0 || maximum_distance <= 0.0 {
    return 0.0
  }
  let mut low = 0.0
  let mut high = maximum_distance
  for _ in 0..<48 {
    let middle = (low + high) / 2.0
    let separation = phase_separation(model, source_depth, middle)
    if separation < observed_separation {
      low = middle
    } else {
      high = middle
    }
  }
  (low + high) / 2.0
}
