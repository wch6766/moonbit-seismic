# Benchmarks

This file records a reproducible local run of the checked-in benchmark command.
It is intended to catch accidental changes in workload output and to provide a
transparent reference point for future optimization work.

## Workload

The command is implemented in `benchmark_workloads.mbt` and executed by
`cmd/benchmark`:

- 4,096 samples at 100 Hz;
- three deterministic synthetic event insertions;
- noise amplitude `0.02`;
- fixed noise seed `1337`;
- one pipeline iteration per command invocation.

The detector currently reports five event candidates for this workload. That is
the observed result of the current detector configuration, not an assumption
that the number of inserted events must equal the number of candidates.

## Observed output

Command:

```text
moon run cmd/benchmark
```

Toolchain: Moon `0.1.20260807`, MoonC `0.10.7+bc794d341`, Windows local
runner, recorded 2026-08-19.

```text
workload=synthetic-4k-event-trace iterations=1 samples=4096 events=5 checksum=1.0002339599117342 dominant_hz=0.2197265625
```

The checksum, processed sample count, event count, and dominant frequency are
stable-value checks for the deterministic workload.

## CLI wall-clock sample

Five consecutive `moon run cmd/benchmark` invocations from the same warm local
workspace measured the complete CLI process, including Moon's cached build
lookup and process startup:

| Run | Wall time |
| ---: | ---: |
| 1 | 696.43 ms |
| 2 | 716.64 ms |
| 3 | 950.82 ms |
| 4 | 858.18 ms |
| 5 | 731.41 ms |

Median: **731.41 ms**. These are machine-specific CLI timings rather than a
portable performance guarantee; rerun the command on the target hardware when
comparing changes.
