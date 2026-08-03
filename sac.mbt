///|
/// Helper to read a 32-bit integer from bytes with specified endianness.
fn read_int32(bytes : Bytes, offset : Int, is_big_endian : Bool) -> Int {
  let b0 = bytes[offset].to_int()
  let b1 = bytes[offset + 1].to_int()
  let b2 = bytes[offset + 2].to_int()
  let b3 = bytes[offset + 3].to_int()
  if is_big_endian {
    (b0 << 24) | (b1 << 16) | (b2 << 8) | b3
  } else {
    b0 | (b1 << 8) | (b2 << 16) | (b3 << 24)
  }
}

///|
/// Helper to read a 32-bit float from bytes with specified endianness.
fn read_float32(bytes : Bytes, offset : Int, is_big_endian : Bool) -> Float {
  let bits = read_int32(bytes, offset, is_big_endian)
  Float::reinterpret_from_int(bits)
}

///|
/// Helper to write a 32-bit integer to a FixedArray[Byte] with specified endianness.
fn write_int32(
  bytes : FixedArray[Byte],
  offset : Int,
  value : Int,
  is_big_endian : Bool,
) -> Unit {
  let b0 = (value & 0xff).to_byte()
  let b1 = ((value >> 8) & 0xff).to_byte()
  let b2 = ((value >> 16) & 0xff).to_byte()
  let b3 = ((value >> 24) & 0xff).to_byte()
  if is_big_endian {
    bytes[offset] = b3
    bytes[offset + 1] = b2
    bytes[offset + 2] = b1
    bytes[offset + 3] = b0
  } else {
    bytes[offset] = b0
    bytes[offset + 1] = b1
    bytes[offset + 2] = b2
    bytes[offset + 3] = b3
  }
}

///|
/// Helper to write a 32-bit float to a FixedArray[Byte] with specified endianness.
fn write_float32(
  bytes : FixedArray[Byte],
  offset : Int,
  value : Float,
  is_big_endian : Bool,
) -> Unit {
  let bits = Float::reinterpret_as_int(value)
  write_int32(bytes, offset, bits, is_big_endian)
}

///|
/// Read an ASCII string from a fixed-length header field in Bytes.
fn read_string(bytes : Bytes, offset : Int, length : Int) -> String {
  let mut s = ""
  for i = 0; i < length; i = i + 1 {
    let b = bytes[offset + i].to_int()
    if b == 0 {
      break
    }
    if b >= 32 && b <= 126 {
      let char_opt = b.to_char()
      match char_opt {
        Some(c) => s = s + c.to_string()
        None => ()
      }
    }
  }

  // Trim spaces
  let mut start = 0
  let mut end = s.length()
  while start < end && s[start] == ' ' {
    start = start + 1
  }
  while end > start && s[end - 1] == ' ' {
    end = end - 1
  }
  s.unsafe_substring(start~, end~)
}

///|
/// Write a string to a fixed-length header field in a FixedArray[Byte], padding with spaces.
fn write_string_field(
  bytes : FixedArray[Byte],
  offset : Int,
  s : String,
  max_length : Int,
) -> Unit {
  let len = s.length()
  for i = 0; i < max_length; i = i + 1 {
    if i < len {
      let char_code = s[i].to_int()
      bytes[offset + i] = char_code.to_byte()
    } else {
      bytes[offset + i] = b'\x20' // space
    }
  }
}

///|
/// Parse a SAC binary file from a byte array.
pub fn parse_sac(bytes : Bytes) -> Result[Trace, String] {
  if bytes.length() < 632 {
    return Err("Byte array too small for SAC header")
  }

  // Detect endianness using nvhdr (offset 304, integer index 6)
  let nvhdr_le = read_int32(bytes, 304, false)
  let (is_big_endian, version) = if nvhdr_le == 6 {
    (false, 6)
  } else {
    let nvhdr_be = read_int32(bytes, 304, true)
    if nvhdr_be == 6 {
      (true, 6)
    } else {
      (false, nvhdr_le)
    }
  }

  if version != 6 {
    return Err("Unsupported SAC version: \{version}")
  }

  let npts = read_int32(bytes, 316, is_big_endian) // integer index 9
  if npts < 0 {
    return Err("Invalid npts in SAC header: \{npts}")
  }

  let expected_size = 632 + npts * 4
  if bytes.length() < expected_size {
    return Err(
      "Byte array too small for data: expected \{expected_size}, got \{bytes.length()}",
    )
  }

  let delta = read_float32(bytes, 0, is_big_endian).to_double()
  let b = read_float32(bytes, 20, is_big_endian).to_double()

  let station = read_string(bytes, 440, 8)
  let channel = read_string(bytes, 608, 8)

  // Read data points
  let data = Array::new(capacity=npts)
  for i = 0; i < npts; i = i + 1 {
    let val = read_float32(bytes, 632 + i * 4, is_big_endian).to_double()
    data.push(val)
  }

  let sampling_rate = if delta <= 0.0 { 1.0 } else { 1.0 / delta }

  Ok(Trace::new(station, channel, sampling_rate, b, data))
}

///|
/// Serialize a Trace into SAC binary bytes.
pub fn to_sac_bytes(trace : Trace, is_big_endian : Bool) -> Bytes {
  let npts = trace.npts()
  let size = 632 + npts * 4
  let bytes = FixedArray::make(size, b'\x00')

  // Fill float headers with standard undefined values (-12345.0)
  for i = 0; i < 70; i = i + 1 {
    write_float32(bytes, i * 4, -12345.0, is_big_endian)
  }

  // Fill integer headers with -12345
  for i = 0; i < 40; i = i + 1 {
    write_int32(bytes, 280 + i * 4, -12345, is_big_endian)
  }

  // Fill string headers with spaces (ASCII 32)
  for i = 0; i < 192; i = i + 1 {
    bytes[440 + i] = b'\x20'
  }

  // Set specific fields
  let delta = Float::from_double(trace.delta())
  write_float32(bytes, 0, delta, is_big_endian)

  let b = Float::from_double(trace.start_time)
  write_float32(bytes, 20, b, is_big_endian)

  // nvhdr (version) = 6
  write_int32(bytes, 304, 6, is_big_endian)

  // npts
  write_int32(bytes, 316, npts, is_big_endian)

  // iftype = 1 (Time series)
  write_int32(bytes, 340, 1, is_big_endian)

  // station
  write_string_field(bytes, 440, trace.station, 8)
  // channel
  write_string_field(bytes, 608, trace.channel, 8)

  // Write data
  for i = 0; i < npts; i = i + 1 {
    let val = Float::from_double(trace.data[i])
    write_float32(bytes, 632 + i * 4, val, is_big_endian)
  }

  bytes.unsafe_reinterpret_as_bytes()
}
