# gperftools CPU profiler reference

This reference is loaded on demand by the `profile-and-optimize` skill
when interpreting profile output or modifying profiling settings. The
critical fact for Lancet2 work is that **most of the upstream gperftools
documentation describes runtime knobs that Lancet2 deliberately bypasses** —
every relevant setting is hardcoded in source. Edit and rebuild to change
anything.


Lancet2 statically links gperftools when the build is configured with `LANCET_PROFILE_MODE=ON` (the `pixi run configure-profile` task sets this, requires `RelWithDebInfo`). The full upstream documentation is at [gperftools.github.io/gperftools/cpuprofile.html](https://gperftools.github.io/gperftools/cpuprofile.html); the salient points for working on Lancet2 are below. Note that **most of the upstream documentation describes runtime knobs (env vars, signals, `LD_PRELOAD`) that Lancet2 deliberately bypasses** — every relevant setting is hardcoded in the source. To change anything, edit and rebuild.

### How profiling starts and stops

The library provides three upstream mechanisms — `LD_PRELOAD`, the `CPUPROFILE` environment variable, and explicit `ProfilerStart()`/`ProfilerStop()` API calls. Lancet2 uses only the third, in `src/lancet/cli/pipeline_runner.cpp`:

```cpp
#ifdef LANCET_PROFILE_MODE
  setenv("CPUPROFILE_PER_THREAD_TIMERS", "1", 1);
  setenv("CPUPROFILE_FREQUENCY", "250", 1);
  auto const timestamp = absl::FormatTime("%Y%m%d%ET%H%M%S", absl::Now(), absl::LocalTimeZone());
  auto const fname = fmt::format("Lancet.cpu_profile.{}.bin", timestamp);
  ProfilerStart(fname.c_str());
#endif
```

`ProfilerStop()` and `ProfilerFlush()` are called from `src/lancet/core/pipeline_executor.cpp` after the pipeline finishes. The implications for users:

- **The output filename is fixed.** `ProfilerStart` takes the filename directly, so `CPUPROFILE` is ignored. The file always lands in the current working directory as `Lancet.cpu_profile.<TIMESTAMP>.bin`. To control where the file goes, control the working directory of the Lancet2 invocation.
- **The two `setenv` calls happen inside the binary, before `ProfilerStart`.** This is the only way these settings take effect, because gperftools reads its env vars during static initialization (before `main()` runs). Setting `CPUPROFILE_PER_THREAD_TIMERS` or `CPUPROFILE_FREQUENCY` from your shell before invoking Lancet2 has no effect.
- **The stack-unwinder is compile-time-selected.** `--enable-frame-pointers` in `cmake/configure_gperftools.sh` decides this; no env var or runtime setting can change it.

To change the sampling frequency (default 250 Hz in the project), the per-thread timer setting, or any other gperftools knob, edit `pipeline_runner.cpp` and rebuild via `pixi run build-profile`. Do not paste `export CPUPROFILE_*=...` recipes from upstream gperftools docs into Lancet2 work; they will silently do nothing.

### What gets sampled

Each sample captures the call stack at the moment SIGPROF fires. The sampler attributes the time slice to whichever function is currently executing on each thread — under threading, total profiled time can exceed wall-clock time (it's CPU-time across all cores).

The project's `CPUPROFILE_PER_THREAD_TIMERS=1` choice ensures every thread has its own SIGPROF timer. Without that, all threads share one timer; samples get biased toward whichever thread happens to be running when SIGPROF fires. For Lancet2's multi-threaded pipeline, per-thread timers are the right setting; that's why it's hardcoded.

### Wall-time vs CPU-time profiling

By default gperftools uses `ITIMER_PROF` (CPU time only — sleeps and I/O waits are excluded from samples). Setting `CPUPROFILE_REALTIME=1` switches to `ITIMER_REAL` (wall-clock — sleeps and I/O get sampled too).

Lancet2 does not currently set this, and the existing analyses are CPU-time analyses. If you specifically need wall-time profiling (you suspect I/O wait is the bottleneck and you want to see it represented in the profile), the change is in `pipeline_runner.cpp`: add `setenv("CPUPROFILE_REALTIME", "1", 1);` alongside the other two `setenv` calls and rebuild. Be aware that this will mix CPU-active and CPU-blocked time in the same profile, which makes hot-function analysis harder; revert before merging unless I/O is genuinely the topic.

### Caveats from the upstream docs

- If the program exits via signal, the profile may be incomplete or empty. Make sure the pipeline reaches its normal completion (`ProfilerStop()` flushes the file).
- If a library was compiled without sufficient symbol info, samples may be charged to the last symbol the profiler could resolve. Lancet2's `RelWithDebInfo` profile build avoids this for project code; vendored libraries (HTSlib, SPOA) may still show up as opaque symbols.
- Profiling on one machine and analyzing on another with different shared libraries produces confusing output. Symbolize on the same host you profiled on.

### Common Lancet2-specific pitfalls

- **Forgetting to use the profiling build.** A Debug or plain Release build won't have the gperftools instrumentation linked in, so `Lancet.cpu_profile.*.bin` won't be produced. Always `pixi run build-profile` first.
- **Looking for the file in the wrong place.** The file lands in the *current working directory of the Lancet2 invocation*, not in the project root, not in the build directory, and not in `/tmp` (unless you `cd /tmp` first). Steps 2 and 7 of this skill use a `cd /tmp/lancet2-profile && ./Lancet2 ...` pattern for exactly this reason.
- **Trying to set `CPUPROFILE` from the shell.** It's silently ignored. The filename is hardcoded to `Lancet.cpu_profile.<TIMESTAMP>.bin`.

