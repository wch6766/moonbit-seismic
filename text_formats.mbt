///|
/// A portable text representation for exchanging waveform samples.
pub struct TextWaveform {
  station : String
  channel : String
  sample_rate : Double
  start_time : Double
  data : Array[Double]
}

///|
/// Construct a text waveform value.
pub fn TextWaveform::new(
  station : String,
  channel : String,
  sample_rate : Double,
  start_time : Double,
  data : Array[Double],
) -> TextWaveform {
  { station, channel, sample_rate, start_time, data }
}

///|
/// Parse a single numeric token with a useful column error.
fn parse_text_number(
  token : String,
  line_number : Int,
) -> Result[Double, String] {
  Ok(@strconv.from_str(token.to_string_view())) catch {
    _ => Err("Invalid waveform number at line \{line_number}: \{token}")
  }
}

///|
/// Parse whitespace-separated samples, ignoring blank and comment lines.
pub fn parse_ascii_samples(text : String) -> Result[Array[Double], String] {
  let result = Array::new()
  let mut line_number = 0
  for line_view in text.split("\n") {
    line_number = line_number + 1
    let line = line_view.trim().to_owned()
    if line == "" || line.has_prefix("#") {
      continue
    }
    for token_view in line.split(" ") {
      let token = token_view.trim().to_owned()
      if token != "" && token != "\r" {
        match parse_text_number(token, line_number) {
          Ok(value) => result.push(value)
          Err(message) => return Err(message)
        }
      }
    }
  }
  Ok(result)
}

///|
/// Parse comma-separated samples with an optional header row.
pub fn parse_csv_samples(
  text : String,
  has_header : Bool,
) -> Result[Array[Double], String] {
  let result = Array::new()
  let mut line_number = 0
  for line_view in text.split("\n") {
    line_number = line_number + 1
    if has_header && line_number == 1 {
      continue
    }
    let line = line_view.trim().to_owned()
    if line == "" || line.has_prefix("#") {
      continue
    }
    let fields = line.split(",")
    for field_view in fields {
      let token = field_view.trim().to_owned()
      match parse_text_number(token, line_number) {
        Ok(value) => result.push(value)
        Err(message) => return Err(message)
      }
    }
  }
  Ok(result)
}

///|
/// Serialize samples as one value per line.
pub fn samples_to_ascii(
  data : Array[Double],
  precision_label : String,
) -> String {
  let mut output = "# moonbit-seismic waveform ".add(precision_label).add("\n")
  for value in data {
    output = output + value.to_string() + "\n"
  }
  output
}

///|
/// Serialize samples as a CSV column with an explicit header.
pub fn samples_to_csv(data : Array[Double], header : String) -> String {
  let mut output = header + "\n"
  for value in data {
    output = output + value.to_string() + "\n"
  }
  output
}

///|
/// Parse a metadata-prefixed CSV waveform.
///
/// The first four rows are `station`, `channel`, `sample_rate`, and
/// `start_time`; remaining rows contain one sample each.
pub fn parse_text_waveform(text : String) -> Result[TextWaveform, String] {
  let lines = Array::new()
  for line_view in text.split("\n") {
    let line = line_view.trim().to_owned()
    if line != "" && !line.has_prefix("#") {
      lines.push(line)
    }
  }
  if lines.length() < 4 {
    return Err("Text waveform requires four metadata lines")
  }
  let station = lines[0]
  let channel = lines[1]
  let sample_rate = match parse_text_number(lines[2], 3) {
    Ok(value) => value
    Err(message) => return Err(message)
  }
  let start_time = match parse_text_number(lines[3], 4) {
    Ok(value) => value
    Err(message) => return Err(message)
  }
  let data = Array::new()
  for index = 4; index < lines.length(); index = index + 1 {
    match parse_text_number(lines[index], index + 1) {
      Ok(value) => data.push(value)
      Err(message) => return Err(message)
    }
  }
  Ok({ station, channel, sample_rate, start_time, data })
}

///|
/// Serialize a waveform in the portable metadata-prefixed text format.
pub fn text_waveform_to_string(waveform : TextWaveform) -> String {
  let mut output = waveform.station + "\n" + waveform.channel + "\n"
  output = output + waveform.sample_rate.to_string() + "\n"
  output = output + waveform.start_time.to_string() + "\n"
  for value in waveform.data {
    output = output + value.to_string() + "\n"
  }
  output
}

///|
/// Convert a trace to the text exchange type.
pub fn trace_to_text_waveform(trace : Trace) -> TextWaveform {
  TextWaveform::new(
    trace.station,
    trace.channel,
    trace.sampling_rate,
    trace.start_time,
    trace.data.copy(),
  )
}

///|
/// Convert a text exchange value to a trace.
pub fn text_waveform_to_trace(waveform : TextWaveform) -> Trace {
  Trace::new(
    waveform.station,
    waveform.channel,
    waveform.sample_rate,
    waveform.start_time,
    waveform.data,
  )
}

///|
/// Serialize a trace as a two-column CSV of time and amplitude.
pub fn trace_to_csv(trace : Trace) -> String {
  let mut output = "time,amplitude\n"
  for index = 0; index < trace.npts(); index = index + 1 {
    let time = trace.start_time + index.to_double() * trace.delta()
    output = output +
      time.to_string() +
      "," +
      trace.data[index].to_string() +
      "\n"
  }
  output
}

///|
/// Parse a two-column time/amplitude CSV and infer sample rate.
pub fn parse_time_amplitude_csv(
  text : String,
  station : String,
  channel : String,
) -> Result[Trace, String] {
  let times = Array::new()
  let values = Array::new()
  let mut line_number = 0
  for line_view in text.split("\n") {
    line_number = line_number + 1
    let line = line_view.trim().to_owned()
    if line == "" || line.has_prefix("#") || line.has_prefix("time,") {
      continue
    }
    let fields = Array::new()
    for field_view in line.split(",") {
      fields.push(field_view.trim().to_owned())
    }
    if fields.length() < 2 {
      return Err("Expected time and amplitude at line {line_number}")
    }
    let time = match parse_text_number(fields[0], line_number) {
      Ok(value) => value
      Err(message) => return Err(message)
    }
    let value = match parse_text_number(fields[1], line_number) {
      Ok(value) => value
      Err(message) => return Err(message)
    }
    times.push(time)
    values.push(value)
  }
  if times.length() < 2 {
    return Err("Time/amplitude CSV requires at least two rows")
  }
  let dt = times[1] - times[0]
  if dt <= 0.0 {
    return Err("CSV times must be strictly increasing")
  }
  for index = 2; index < times.length(); index = index + 1 {
    if (times[index] - times[index - 1] - dt).abs() > dt * 0.01 {
      return Err("CSV time steps must be uniform")
    }
  }
  Ok(Trace::new(station, channel, 1.0 / dt, times[0], values))
}
