///|
/// SAC byte-order marker.
pub enum SacEndian {
  SacLittleEndian
  SacBigEndian
}

///|
/// Return the little-endian marker.
pub fn sac_little_endian() -> SacEndian {
  SacLittleEndian
}

///|
/// Return the big-endian marker.
pub fn sac_big_endian() -> SacEndian {
  SacBigEndian
}

///|
/// Parsed SAC header fields needed by downstream processing.
pub struct SacHeaderInfo {
  endian : SacEndian
  version : Int
  sample_count : Int
  delta : Double
  begin : Double
  station : String
  channel : String
  header_size : Int
}

///|
/// Detect SAC byte order from the NVHDR field.
pub fn detect_sac_endian(bytes : Bytes) -> Result[SacEndian, String] {
  if bytes.length() < 632 {
    return Err("SAC header requires at least 632 bytes")
  }
  let little = read_int32(bytes, 304, false)
  if little == 6 {
    return Ok(SacLittleEndian)
  }
  let big = read_int32(bytes, 304, true)
  if big == 6 {
    Ok(SacBigEndian)
  } else {
    Err("SAC NVHDR is not a supported version")
  }
}

///|
/// Extract a validated SAC header.
pub fn parse_sac_header(bytes : Bytes) -> Result[SacHeaderInfo, String] {
  let endian = match detect_sac_endian(bytes) {
    Ok(value) => value
    Err(message) => return Err(message)
  }
  let is_big = endian is SacBigEndian
  let version = read_int32(bytes, 304, is_big)
  let sample_count = read_int32(bytes, 316, is_big)
  if sample_count < 0 {
    return Err("SAC sample count cannot be negative")
  }
  let delta = read_float32(bytes, 0, is_big).to_double()
  let begin = read_float32(bytes, 20, is_big).to_double()
  Ok({
    endian,
    version,
    sample_count,
    delta,
    begin,
    station: read_string(bytes, 440, 8),
    channel: read_string(bytes, 608, 8),
    header_size: 632,
  })
}

///|
/// Return the physical size required by a SAC header and sample count.
pub fn sac_record_size(sample_count : Int) -> Int {
  if sample_count < 0 {
    632
  } else {
    632 + sample_count * 4
  }
}

///|
/// Validate that a byte array contains the complete SAC record.
pub fn validate_sac_record(bytes : Bytes) -> Result[Unit, String] {
  let header = match parse_sac_header(bytes) {
    Ok(value) => value
    Err(message) => return Err(message)
  }
  let expected = sac_record_size(header.sample_count)
  if bytes.length() < expected {
    Err("SAC record is shorter than its header declares")
  } else {
    Ok(())
  }
}

///|
/// Return the sample interval from a SAC record.
pub fn sac_delta(bytes : Bytes) -> Result[Double, String] {
  match parse_sac_header(bytes) {
    Ok(header) => Ok(header.delta)
    Err(message) => Err(message)
  }
}

///|
/// Return all SAC data samples as doubles.
pub fn sac_samples(bytes : Bytes) -> Result[Array[Double], String] {
  let header = match parse_sac_header(bytes) {
    Ok(value) => value
    Err(message) => return Err(message)
  }
  let expected = sac_record_size(header.sample_count)
  if bytes.length() < expected {
    return Err("SAC record is truncated")
  }
  let is_big = header.endian is SacBigEndian
  let result = Array::new(capacity=header.sample_count)
  for index = 0; index < header.sample_count; index = index + 1 {
    result.push(read_float32(bytes, 632 + index * 4, is_big).to_double())
  }
  Ok(result)
}

///|
/// Return a compact header-plus-data summary.
pub fn sac_summary(bytes : Bytes) -> Result[DescriptiveStats, String] {
  match sac_samples(bytes) {
    Ok(data) => Ok(describe_samples(data))
    Err(message) => Err(message)
  }
}

///|
/// Serialize a trace with a selected endian marker.
pub fn serialize_sac(trace : Trace, endian : SacEndian) -> Bytes {
  to_sac_bytes(trace, endian is SacBigEndian)
}
