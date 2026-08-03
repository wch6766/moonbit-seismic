///|
/// Helper to read a 16-bit unsigned integer in big-endian byte order.
fn read_uint16_be(bytes : Bytes, offset : Int) -> Int {
  let b0 = bytes[offset].to_int()
  let b1 = bytes[offset + 1].to_int()
  (b0 << 8) | b1
}

///|
/// Helper to read a 16-bit signed integer in big-endian byte order.
fn read_int16_be(bytes : Bytes, offset : Int) -> Int {
  let u = read_uint16_be(bytes, offset)
  if (u & 0x8000) != 0 {
    u - 65536
  } else {
    u
  }
}

///|
/// Parse a single miniSEED data record.
/// Supports uncompressed 32-bit float (is_float = true) or uncompressed 32-bit integer (is_float = false) data.
pub fn parse_miniseed_record(
  bytes : Bytes,
  is_float : Bool,
) -> Result[Trace, String] {
  if bytes.length() < 48 {
    return Err("Bytes too small for miniSEED FSDH")
  }

  // miniSEED headers are big-endian by definition
  let station = read_string(bytes, 8, 5)
  let channel = read_string(bytes, 15, 3)

  let _year = read_uint16_be(bytes, 20)
  let jday = read_uint16_be(bytes, 22)
  let hour = bytes[24].to_int()
  let min = bytes[25].to_int()
  let sec = bytes[26].to_int()
  let msec = read_uint16_be(bytes, 28) / 10 // 0.0001s to ms

  let num_samples = read_uint16_be(bytes, 30)
  let rate_factor = read_int16_be(bytes, 32)
  let rate_multiplier = read_int16_be(bytes, 34)

  let data_offset = read_uint16_be(bytes, 44)

  if data_offset + num_samples * 4 > bytes.length() {
    return Err("Record size too small for specified sample count and offset")
  }

  // Compute sampling rate from factor (F) and multiplier (M)
  let sampling_rate = if rate_factor > 0 && rate_multiplier > 0 {
    rate_factor.to_double() * rate_multiplier.to_double()
  } else if rate_factor > 0 && rate_multiplier < 0 {
    rate_factor.to_double() / (-rate_multiplier).to_double()
  } else if rate_factor < 0 && rate_multiplier > 0 {
    (-rate_factor).to_double() / rate_multiplier.to_double()
  } else if rate_factor < 0 && rate_multiplier < 0 {
    1.0 / ((-rate_factor).to_double() * (-rate_multiplier).to_double())
  } else {
    1.0
  }

  let data = Array::new(capacity=num_samples)
  for i = 0; i < num_samples; i = i + 1 {
    let offset = data_offset + i * 4
    let val = if is_float {
      read_float32(bytes, offset, true).to_double()
    } else {
      read_int32(bytes, offset, true).to_double()
    }
    data.push(val)
  }

  // Approximate start time in seconds relative to the beginning of the year
  let start_time = (jday * 24 * 3600 + hour * 3600 + min * 60 + sec).to_double() +
    msec.to_double() / 1000.0

  Ok(Trace::new(station, channel, sampling_rate, start_time, data))
}

///|
/// Serialize a Trace into miniSEED binary bytes.
pub fn to_miniseed_bytes(
  trace : Trace,
  sequence_number : Int,
  is_float : Bool,
) -> Bytes {
  let npts = trace.npts()
  let data_offset = 64
  let required_size = data_offset + npts * 4
  let record_size = if required_size % 256 == 0 {
    required_size
  } else {
    (required_size / 256 + 1) * 256
  }

  let bytes = FixedArray::make(record_size, b'\x00')

  // Sequence number (6 ASCII chars, e.g. "000001")
  let seq_str = sequence_number.to_string()
  let seq_len = seq_str.length()
  for i = 0; i < 6; i = i + 1 {
    if i < 6 - seq_len {
      bytes[i] = b'\x30' // '0'
    } else {
      let idx = i - (6 - seq_len)
      bytes[i] = seq_str[idx].to_int().to_byte()
    }
  }

  bytes[6] = b'\x44' // 'D'
  bytes[7] = b'\x20' // space

  write_string_field_bytes(bytes, 8, trace.station, 5)
  write_string_field_bytes(bytes, 15, trace.channel, 3)
  write_string_field_bytes(bytes, 13, "", 2)
  write_string_field_bytes(bytes, 18, "XX", 2)

  // Derive start time parameters
  let total_sec = trace.start_time.to_int()
  let msec = ((trace.start_time - total_sec.to_double()) * 10000.0).to_int()
  let jday = total_sec / (24 * 3600)
  let rem = total_sec % (24 * 3600)
  let hour = rem / 3600
  let rem2 = rem % 3600
  let min = rem2 / 60
  let sec = rem2 % 60

  write_uint16_be_bytes(bytes, 20, 2026) // default year
  write_uint16_be_bytes(bytes, 22, jday + 1)
  bytes[24] = hour.to_byte()
  bytes[25] = min.to_byte()
  bytes[26] = sec.to_byte()
  bytes[27] = b'\x00'
  write_uint16_be_bytes(bytes, 28, msec)

  write_uint16_be_bytes(bytes, 30, npts)

  // Format sample rate
  let rate_int = trace.sampling_rate.to_int()
  if rate_int.to_double() == trace.sampling_rate {
    write_int16_be_bytes(bytes, 32, rate_int)
    write_int16_be_bytes(bytes, 34, 1)
  } else {
    let factor = (trace.sampling_rate * 100.0).to_int()
    write_int16_be_bytes(bytes, 32, factor)
    write_int16_be_bytes(bytes, 34, -100)
  }

  write_uint16_be_bytes(bytes, 44, data_offset)

  // Write data
  for i = 0; i < npts; i = i + 1 {
    let offset = data_offset + i * 4
    if is_float {
      let val = Float::from_double(trace.data[i])
      let bits = Float::reinterpret_as_int(val)
      write_int32_bytes(bytes, offset, bits, true)
    } else {
      let bits = trace.data[i].to_int()
      write_int32_bytes(bytes, offset, bits, true)
    }
  }

  bytes.unsafe_reinterpret_as_bytes()
}

///|
fn write_string_field_bytes(
  bytes : FixedArray[Byte],
  offset : Int,
  s : String,
  max_length : Int,
) -> Unit {
  let len = s.length()
  for i = 0; i < max_length; i = i + 1 {
    if i < len {
      bytes[offset + i] = s[i].to_int().to_byte()
    } else {
      bytes[offset + i] = b'\x20'
    }
  }
}

///|
fn write_uint16_be_bytes(
  bytes : FixedArray[Byte],
  offset : Int,
  value : Int,
) -> Unit {
  bytes[offset] = (value >> 8).to_byte()
  bytes[offset + 1] = (value & 0xff).to_byte()
}

///|
fn write_int16_be_bytes(
  bytes : FixedArray[Byte],
  offset : Int,
  value : Int,
) -> Unit {
  write_uint16_be_bytes(bytes, offset, value)
}

///|
fn write_int32_bytes(
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
