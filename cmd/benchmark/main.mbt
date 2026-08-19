///|
/// Run the deterministic benchmark workload used in the documentation.
fn main {
  let workload = @seismic.acceptance_workload()
  let summary = @seismic.run_benchmark_workload(workload, 1)
  println(@seismic.benchmark_summary_text(summary))
}
