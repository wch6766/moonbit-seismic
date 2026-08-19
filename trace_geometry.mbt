///|
/// Geographic station location in degrees.
pub struct StationPosition {
  station : String
  latitude : Double
  longitude : Double
  elevation_km : Double
}

///|
/// A small grid-search location result.
pub struct LocationEstimate {
  latitude : Double
  longitude : Double
  depth_km : Double
  origin_time : Double
  residual : Double
}

///|
/// Construct a station position.
pub fn StationPosition::new(
  station : String,
  latitude : Double,
  longitude : Double,
  elevation_km : Double,
) -> StationPosition {
  { station, latitude, longitude, elevation_km }
}

///|
/// Approximate surface distance using the haversine formula.
pub fn surface_distance_km(
  left : StationPosition,
  right : StationPosition,
) -> Double {
  let earth_radius = 6371.0
  let dlat = (right.latitude - left.latitude) * @math.PI / 180.0
  let dlon = (right.longitude - left.longitude) * @math.PI / 180.0
  let lat1 = left.latitude * @math.PI / 180.0
  let lat2 = right.latitude * @math.PI / 180.0
  let a = @math.sin(dlat / 2.0) * @math.sin(dlat / 2.0) +
    @math.cos(lat1) *
    @math.cos(lat2) *
    @math.sin(dlon / 2.0) *
    @math.sin(dlon / 2.0)
  2.0 * earth_radius * @math.asin(a.sqrt())
}

///|
/// Return a clockwise azimuth from one station to another.
pub fn azimuth_degrees(
  left : StationPosition,
  right : StationPosition,
) -> Double {
  let lat1 = left.latitude * @math.PI / 180.0
  let lat2 = right.latitude * @math.PI / 180.0
  let dlon = (right.longitude - left.longitude) * @math.PI / 180.0
  let y = @math.sin(dlon) * @math.cos(lat2)
  let x = @math.cos(lat1) * @math.sin(lat2) -
    @math.sin(lat1) * @math.cos(lat2) * @math.cos(dlon)
  let degrees = @math.atan2(y, x) * 180.0 / @math.PI
  if degrees < 0.0 {
    degrees + 360.0
  } else {
    degrees
  }
}

///|
/// Return the RMS residual between observed and predicted arrivals.
pub fn arrival_residual(
  observed : Array[Double],
  predicted : Array[Double],
) -> Double {
  let count = if observed.length() < predicted.length() {
    observed.length()
  } else {
    predicted.length()
  }
  if count == 0 {
    return 0.0
  }
  let mut sum = 0.0
  for index = 0; index < count; index = index + 1 {
    let difference = observed[index] - predicted[index]
    sum += difference * difference
  }
  (sum / count.to_double()).sqrt()
}

///|
/// Evaluate a candidate source against station observations.
pub fn location_residual(
  model : VelocityModel,
  source_depth : Double,
  source_latitude : Double,
  source_longitude : Double,
  origin_time : Double,
  stations : Array[StationPosition],
  observed : Array[Double],
) -> Double {
  let predicted = Array::new()
  let source = StationPosition::new(
    "source", source_latitude, source_longitude, 0.0,
  )
  for station in stations {
    let distance = surface_distance_km(source, station)
    predicted.push(
      origin_time + travel_time(model, source_depth, distance, false),
    )
  }
  arrival_residual(observed, predicted)
}

///|
/// Search a rectangular epicenter grid at one depth.
pub fn grid_search_location(
  model : VelocityModel,
  source_depth : Double,
  latitude_min : Double,
  latitude_max : Double,
  longitude_min : Double,
  longitude_max : Double,
  grid_size : Int,
  stations : Array[StationPosition],
  observed : Array[Double],
) -> LocationEstimate {
  let count = if grid_size < 2 { 2 } else { grid_size }
  let mut best = {
    latitude: latitude_min,
    longitude: longitude_min,
    depth_km: source_depth,
    origin_time: minimum_sample(observed),
    residual: 1.0e300,
  }
  for lat_index = 0; lat_index < count; lat_index = lat_index + 1 {
    let latitude = latitude_min +
      (latitude_max - latitude_min) *
      lat_index.to_double() /
      (count - 1).to_double()
    for lon_index = 0; lon_index < count; lon_index = lon_index + 1 {
      let longitude = longitude_min +
        (longitude_max - longitude_min) *
        lon_index.to_double() /
        (count - 1).to_double()
      let origin = minimum_sample(observed)
      let residual = location_residual(
        model, source_depth, latitude, longitude, origin, stations, observed,
      )
      if residual < best.residual {
        best = {
          latitude,
          longitude,
          depth_km: source_depth,
          origin_time: origin,
          residual,
        }
      }
    }
  }
  best
}
