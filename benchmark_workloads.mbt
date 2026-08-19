///|
/// Reproducible benchmark workload description.
pub struct BenchmarkWorkload {
  name : String
  sample_count : Int
  sample_rate : Double
  event_count : Int
  noise_amplitude : Double
}

///|
/// Work completed by a benchmark run.
pub struct BenchmarkSummary {
  workload : String
  iterations : Int
  samples_processed : Int
  events_found : Int
  checksum : Double
  dominant_frequency : Double
}

///|
/// Construct a benchmark workload.
pub fn BenchmarkWorkload::new(
  name : String,
  sample_count : Int,
  sample_rate : Double,
  event_count : Int,
  noise_amplitude : Double,
) -> BenchmarkWorkload {
  { name, sample_count, sample_rate, event_count, noise_amplitude }
}

///|
/// Return the standard acceptance benchmark workload.
pub fn acceptance_workload() -> BenchmarkWorkload {
  {
    name: "synthetic-4k-event-trace",
    sample_count: 4096,
    sample_rate: 100.0,
    event_count: 3,
    noise_amplitude: 0.02,
  }
}

///|
/// Generate a deterministic benchmark trace.
pub fn benchmark_trace(workload : BenchmarkWorkload) -> Trace {
  let data = deterministic_noise(
    workload.sample_count,
    workload.noise_amplitude,
    1337,
  )
  for event_index = 0
      event_index < workload.event_count
      event_index = event_index + 1 {
    add_synthetic_event(
      data,
      workload.sample_rate,
      SyntheticEvent::new(
        600 + event_index * 900,
        1.0 + event_index.to_double() * 0.2,
        4.0 + event_index.to_double(),
        160,
        0.0,
      ),
    )
  }
  Trace::new("BENCH", "BHZ", workload.sample_rate, 0.0, data)
}

///|
/// Run the standard workload repeatedly and return a non-optimizable summary.
pub fn run_benchmark_workload(
  workload : BenchmarkWorkload,
  iterations : Int,
) -> BenchmarkSummary {
  let count = if iterations < 1 { 1 } else { iterations }
  let mut samples = 0
  let mut events = 0
  let mut checksum = 0.0
  let mut dominant = 0.0
  for _ in 0..<count {
    let trace = benchmark_trace(workload)
    let config = default_processing_config()
    let result = analyze_trace(trace, config)
    match result {
      Ok(value) => {
        samples += value.processed.npts()
        events += value.events.length()
        checksum += value.processed.mean() + value.quality.score
        dominant = value.dominant_frequency
      }
      Err(_) => checksum = checksum - 1.0
    }
  }
  {
    workload: workload.name,
    iterations: count,
    samples_processed: samples,
    events_found: events,
    checksum,
    dominant_frequency: dominant,
  }
}

///|
/// Return a stable benchmark result line for documentation and CI logs.
pub fn benchmark_summary_text(summary : BenchmarkSummary) -> String {
  "workload=" +
  summary.workload +
  " iterations=" +
  summary.iterations.to_string() +
  " samples=" +
  summary.samples_processed.to_string() +
  " events=" +
  summary.events_found.to_string() +
  " checksum=" +
  summary.checksum.to_string() +
  " dominant_hz=" +
  summary.dominant_frequency.to_string()
}
