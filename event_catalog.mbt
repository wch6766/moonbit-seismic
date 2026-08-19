///|
/// A catalog entry combining an origin estimate and a network event.
pub struct CatalogEntry {
  id : Int
  origin_time : Double
  latitude : Double
  longitude : Double
  depth_km : Double
  magnitude : Double
  station_count : Int
  confidence : Double
}

///|
/// Mutable in-memory event catalog for small workflows.
pub struct EventCatalog {
  entries : Array[CatalogEntry]
}

///|
/// Create an empty catalog.
pub fn EventCatalog::new() -> EventCatalog {
  { entries: [] }
}

///|
/// Add an entry and return its assigned identifier.
pub fn EventCatalog::add(
  self : EventCatalog,
  origin_time : Double,
  latitude : Double,
  longitude : Double,
  depth_km : Double,
  magnitude : Double,
  station_count : Int,
  confidence : Double,
) -> Int {
  let id = self.entries.length()
  self.entries.push({
    id,
    origin_time,
    latitude,
    longitude,
    depth_km,
    magnitude,
    station_count,
    confidence,
  })
  id
}

///|
/// Return catalog entries in a time range.
pub fn EventCatalog::between(
  self : EventCatalog,
  start : Double,
  end : Double,
) -> Array[CatalogEntry] {
  self.entries.filter(entry => {
    entry.origin_time >= start && entry.origin_time <= end
  })
}

///|
/// Return entries meeting a magnitude threshold.
pub fn EventCatalog::above_magnitude(
  self : EventCatalog,
  minimum : Double,
) -> Array[CatalogEntry] {
  self.entries.filter(entry => entry.magnitude >= minimum)
}

///|
/// Return the average magnitude of the catalog.
pub fn EventCatalog::mean_magnitude(self : EventCatalog) -> Double {
  let values = self.entries.map(entry => entry.magnitude)
  mean_samples(values)
}

///|
/// Serialize catalog entries as a stable CSV table.
pub fn EventCatalog::to_csv(self : EventCatalog) -> String {
  let mut output = "id,origin_time,latitude,longitude,depth_km,magnitude,stations,confidence\n"
  for entry in self.entries {
    output = output +
      entry.id.to_string() +
      "," +
      entry.origin_time.to_string() +
      "," +
      entry.latitude.to_string() +
      "," +
      entry.longitude.to_string() +
      "," +
      entry.depth_km.to_string() +
      "," +
      entry.magnitude.to_string() +
      "," +
      entry.station_count.to_string() +
      "," +
      entry.confidence.to_string() +
      "\n"
  }
  output
}

///|
/// Convert network events into catalog entries at a supplied location.
pub fn catalog_from_network_events(
  events : Array[NetworkEvent],
  latitude : Double,
  longitude : Double,
  depth_km : Double,
) -> EventCatalog {
  let catalog = EventCatalog::new()
  for event in events {
    ignore(
      catalog.add(
        event.origin_time,
        latitude,
        longitude,
        depth_km,
        event.magnitude,
        event.station_count,
        event.confidence,
      ),
    )
  }
  catalog
}
