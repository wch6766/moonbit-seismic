///|
/// Parsed fields of a fixed-section miniSEED data header.
pub struct MiniSeedHeader {
  sequence : String
  station : String
  location : String
  channel : String
  network : String
  year : Int
  julian_day : Int
  hour : Int
  minute : Int
  second : Int
  fractional_seconds : Int
  sample_count : Int
  sample_rate : Double
  data_offset : Int
}

///|
/// Convert miniSEED factor/multiplier pairs to a positive sample rate.
pub fn decode_miniseed_rate(factor : Int, multiplier : Int) -> Double {
  if factor > 0 && multiplier > 0 {
    factor.to_double() * multiplier.to_double()
  } else if factor > 0 && multiplier < 0 {
    factor.to_double() / (-multiplier).to_double()
  } else if factor < 0 && multiplier > 0 {
    (-factor).to_double() / multiplier.to_double()
  } else if factor < 0 && multiplier < 0 {
    1.0 / ((-factor).to_double() * (-multiplier).to_double())
  } else {
    0.0
  }
}

///|
/// Encode a positive sample rate as a compact factor/multiplier pair.
pub fn encode_miniseed_rate(rate : Double) -> (Int, Int) {
  if rate <= 0.0 {
    return (0, 0)
  }
  let rounded = rate.round()
  if (rounded - rate).abs() < 0.0000001 && rounded <= 32767.0 {
    (rounded.to_int(), 1)
  } else {
    let factor = (rate * 10000.0).round().to_int()
    (factor, -10000)
  }
}

///|
/// Parse a miniSEED fixed-section header without decoding samples.
pub fn parse_miniseed_header(bytes : Bytes) -> Result[MiniSeedHeader, String] {
  if bytes.length() < 48 {
    return Err("miniSEED header requires at least 48 bytes")
  }
  let factor = read_int16_be(bytes, 32)
  let multiplier = read_int16_be(bytes, 34)
  let data_offset = read_uint16_be(bytes, 44)
  let sample_count = read_uint16_be(bytes, 30)
  if data_offset < 48 || data_offset > bytes.length() {
    return Err("miniSEED data offset is outside the record")
  }
  Ok({
    sequence: read_string(bytes, 0, 6),
    station: read_string(bytes, 8, 5),
    location: read_string(bytes, 13, 2),
    channel: read_string(bytes, 15, 3),
    network: read_string(bytes, 18, 2),
    year: read_uint16_be(bytes, 20),
    julian_day: read_uint16_be(bytes, 22),
    hour: bytes[24].to_int(),
    minute: bytes[25].to_int(),
    second: bytes[26].to_int(),
    fractional_seconds: read_uint16_be(bytes, 28),
    sample_count,
    sample_rate: decode_miniseed_rate(factor, multiplier),
    data_offset,
  })
}

///|
/// Validate a fixed-width uncompressed miniSEED record.
pub fn validate_miniseed_record(
  bytes : Bytes,
) -> Result[MiniSeedHeader, String] {
  let header = match parse_miniseed_header(bytes) {
    Ok(value) => value
    Err(message) => return Err(message)
  }
  let expected = header.data_offset + header.sample_count * 4
  if expected > bytes.length() {
    Err("miniSEED record is truncated for its sample count")
  } else {
    Ok(header)
  }
}

///|
/// Decode all uncompressed float or integer samples from a miniSEED record.
pub fn miniseed_samples(
  bytes : Bytes,
  is_float : Bool,
) -> Result[Array[Double], String] {
  let header = match validate_miniseed_record(bytes) {
    Ok(value) => value
    Err(message) => return Err(message)
  }
  let result = Array::new(capacity=header.sample_count)
  for index = 0; index < header.sample_count; index = index + 1 {
    let offset = header.data_offset + index * 4
    let value = if is_float {
      read_float32(bytes, offset, true).to_double()
    } else {
      read_int32(bytes, offset, true).to_double()
    }
    result.push(value)
  }
  Ok(result)
}

///|
/// Return the nominal record size rounded to a power-of-two frame.
pub fn miniseed_record_size(sample_count : Int, data_offset : Int) -> Int {
  let count = if sample_count < 0 { 0 } else { sample_count }
  let offset = if data_offset < 48 { 48 } else { data_offset }
  let needed = offset + count * 4
  if needed % 256 == 0 {
    needed
  } else {
    (needed / 256 + 1) * 256
  }
}

///|
/// Convert a Julian day and clock fields to seconds from the year start.
pub fn miniseed_seconds_from_year(
  julian_day : Int,
  hour : Int,
  minute : Int,
  second : Int,
  fractional_100us : Int,
) -> Double {
  let day = if julian_day < 1 { 1 } else { julian_day }
  (day - 1).to_double() * 86400.0 +
  hour.to_double() * 3600.0 +
  minute.to_double() * 60.0 +
  second.to_double() +
  fractional_100us.to_double() / 10000.0
}

///|
/// Convert a trace start time in seconds to miniSEED clock fields.
pub fn miniseed_clock_from_seconds(
  seconds : Double,
) -> (Int, Int, Int, Int, Int) {
  let safe = if seconds < 0.0 { 0.0 } else { seconds }
  let julian_day = (safe / 86400.0).to_int() + 1
  let remainder = safe % 86400.0
  let hour = (remainder / 3600.0).to_int()
  let minute = ((remainder - hour.to_double() * 3600.0) / 60.0).to_int()
  let second = (remainder -
  hour.to_double() * 3600.0 -
  minute.to_double() * 60.0).to_int()
  let fractional = ((remainder - remainder.to_int().to_double()) * 10000.0)
    .round()
    .to_int()
  (julian_day, hour, minute, second, fractional)
}

///|
/// Return the duration represented by a miniSEED sample count.
pub fn miniseed_duration(sample_count : Int, sample_rate : Double) -> Double {
  if sample_count <= 0 || sample_rate <= 0.0 {
    0.0
  } else {
    (sample_count - 1).to_double() / sample_rate
  }
}
