# Reference-engine case bundles

Native input bundles executed by the controlled v2 reference service
(`run_reference_calculation`). Each directory holds `case.json` (declared
metadata; not independently verified physics) plus the engine-native inputs.
Raw outputs, logs and provenance are archived by the service under the
experiment state directory and referenced by run ID from validation records.

- `calibration/`: workload-envelope cases used to measure reference-engine
  throughput during performance calibration (decision 0039). Their physics
  outputs may later be reused as evaluation evidence only with an explicit
  role record.
