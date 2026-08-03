///|
/// Represents a single flat layer in a 1D velocity model.
pub struct Layer {
  depth : Double // depth to top of layer (km)
  vp : Double // P-wave velocity (km/s)
  vs : Double // S-wave velocity (km/s)
}

///|
/// Create a new Layer.
pub fn Layer::new(depth : Double, vp : Double, vs : Double) -> Layer {
  { depth, vp, vs }
}

///|
/// Represents a 1D velocity model consisting of multiple layers.
pub struct VelocityModel {
  name : String
  layers : Array[Layer]
}

///|
/// Create a new VelocityModel.
pub fn VelocityModel::new(
  name : String,
  layers : Array[Layer],
) -> VelocityModel {
  { name, layers }
}

///|
/// Get the P or S wave velocity at a specific depth in the 1D model.
pub fn get_velocity_at_depth(
  model : VelocityModel,
  depth : Double,
  is_s_wave : Bool,
) -> Double {
  let n = model.layers.length()
  if n == 0 {
    return 1.0 // safety fallback
  }

  let mut active_layer = model.layers[0]
  for i = 1; i < n; i = i + 1 {
    if depth >= model.layers[i].depth {
      active_layer = model.layers[i]
    } else {
      break
    }
  }

  if is_s_wave {
    active_layer.vs
  } else {
    active_layer.vp
  }
}

///|
/// Estimate the travel time of P or S wave between a source and a station
/// using a path-integral slowness approximation along the straight-line ray.
pub fn travel_time(
  model : VelocityModel,
  source_depth : Double,
  distance : Double,
  is_s_wave : Bool,
) -> Double {
  let path_len = (distance * distance + source_depth * source_depth).sqrt()
  if path_len == 0.0 {
    return 0.0
  }

  // Integrate slowness over 100 steps along the straight ray path
  let steps = 100
  let mut total_slowness = 0.0
  for i = 0; i < steps; i = i + 1 {
    let frac = i.to_double() / steps.to_double()
    // depth along the path from source (frac = 0) to receiver (frac = 1, depth = 0)
    let z = source_depth * (1.0 - frac)
    let vel = get_velocity_at_depth(model, z, is_s_wave)
    total_slowness = total_slowness + 1.0 / vel
  }

  let avg_slowness = total_slowness / steps.to_double()
  path_len * avg_slowness
}

///|
/// Construct the standard IASPEI91 global velocity model.
pub fn iaspei91() -> VelocityModel {
  let layers = [
    Layer::new(0.0, 5.8, 3.36),
    Layer::new(20.0, 6.5, 3.75),
    Layer::new(35.0, 8.04, 4.47),
    Layer::new(210.0, 8.3, 4.5),
    Layer::new(410.0, 9.13, 4.9),
  ]
  VelocityModel::new("IASPEI91", layers)
}
