#!/usr/bin/env python3
"""Prepare exact rational blocks for the native beta-series backend."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import pathlib
import re
import selectors
import shutil
import signal
import subprocess
import sys
import time


VERSION = "DirectRationalSeries-1.0"
active_process: subprocess.Popen | None = None


def sha(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: pathlib.Path, value: object) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    os.replace(temporary, path)


def stop(signum, _frame):
    global active_process
    if active_process is not None and active_process.poll() is None:
        active_process.terminate()
    raise SystemExit(128 + signum)


def log(text: str) -> None:
    print(time.strftime("[%Y-%m-%dT%H:%M:%S] ") + text, flush=True)


def invoke(command: list[str], output: pathlib.Path, seconds: int) -> None:
    global active_process

    started = time.monotonic()
    with output.open("a", encoding="utf-8") as handle:
        handle.write("COMMAND " + repr(command) + "\n")
        handle.flush()

        active_process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        assert active_process.stdout is not None

        selector = selectors.DefaultSelector()
        selector.register(active_process.stdout, selectors.EVENT_READ)

        try:
            while selector.get_map():
                if time.monotonic() - started > seconds:
                    active_process.kill()
                    active_process.wait()
                    raise RuntimeError(
                        f"Native operation exceeded {seconds} seconds: "
                        f"{command[0]}; see {output}"
                    )

                for key, _ in selector.select(1):
                    chunk = os.read(key.fd, 65536)
                    if not chunk:
                        selector.unregister(key.fileobj)
                        continue

                    text = chunk.decode("utf-8", errors="replace")
                    handle.write(text)
                    handle.flush()
                    sys.stdout.write(text)
                    sys.stdout.flush()

            code = active_process.wait()
            if code:
                raise RuntimeError(
                    f"Native operation exited with code {code}; see {output}"
                )
        finally:
            selector.close()
            if active_process.poll() is None:
                active_process.kill()
                active_process.wait()
            active_process = None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run", type=pathlib.Path)
    parser.add_argument("--block-seconds", type=int, default=3600)
    args = parser.parse_args()

    run = args.run.resolve()
    source = run / "source"
    checkpoints = run / "checkpoints"

    if not 300 <= args.block_seconds <= 86400:
        raise RuntimeError("block-seconds must be between 300 and 86400")

    native = checkpoints / "native"
    native.mkdir(exist_ok=True)

    backend_sources = [
        "rational_dag.cpp",
        "jet_dag.cpp",
        "reference_dag.cpp",
        "exact_rf.hpp",
        "check_fixed_reference.cpp",
    ]
    signature = hashlib.sha256(
        "".join(name + sha(source / name) for name in backend_sources).encode()
    ).hexdigest()

    blocks = []
    reused = 0

    for projector in ["P1", "P2"]:
        folder = run / "collected" / projector
        index = (folder / "blocks.tsv").read_text().splitlines()
        if not index or index[0] != "block\tnf_power\tmonomial":
            raise RuntimeError("unexpected block index header")

        for line in index[1:]:
            fields = line.split("\t")
            if len(fields) != 3:
                raise RuntimeError("malformed block index")

            block, nf_power = int(fields[0]), int(fields[1])
            if block < 0 or nf_power not in (0, 1, 2):
                raise RuntimeError("invalid block identifier")

            block_source = folder / f"rational_{block:05d}.wl"
            source_hash = sha(block_source)
            cached = (
                checkpoints
                / "reduced"
                / projector
                / f"block_{block:05d}.wl"
            )

            if cached.exists():
                text = cached.read_text()
                match = re.search(
                    '"SourceBlockHash"\\s*->\\s*"([0-9a-f]{64})"',
                    text,
                )
                if not match or match.group(1) != source_hash:
                    raise RuntimeError(
                        "source hash mismatch for cached block " + str(cached)
                    )
                if "GGAA-Unified-RealMPL22-0.1" not in text:
                    raise RuntimeError("unexpected checkpoint version")
                reused += 1
                continue

            blocks.append(
                (projector, block, nf_power, block_source, source_hash)
            )

    # Put the two historically troublesome blocks first, so failures show up early.
    blocks.sort(
        key=lambda value: (
            0
            if value[0] == "P1" and value[1] == 1438
            else 1
            if value[0] == "P2" and value[1] == 6017
            else 2,
            value[0],
            value[1],
        )
    )

    log(
        f"Reusing {reused} checked Wolfram block files; "
        f"{len(blocks)} original blocks need native rational series."
    )

    records = []
    for number, (projector, block, nf_power, block_source, source_hash) in enumerate(
        blocks, 1
    ):
        output = native / projector / f"block_{block:05d}"
        output.mkdir(parents=True, exist_ok=True)
        audit_path = output / "audit.json"

        if audit_path.exists():
            audit = json.loads(audit_path.read_text())
            if not (
                audit.get("passed") is True
                and audit.get("source_sha256") == source_hash
                and audit.get("code_signature") == signature
            ):
                raise RuntimeError("stale native checkpoint " + str(output))

            for record in audit["files"].values():
                path = output / record["name"]
                if not path.is_file() or sha(path) != record["sha256"]:
                    raise RuntimeError("changed native output " + str(path))

            records.append(audit)
            if number % 100 == 0:
                log(f"Reusing native rational series {number}/{len(blocks)}")
            continue

        if shutil.disk_usage(run).free < 512 * 1024**2:
            raise RuntimeError(
                "Less than 512 MiB free disk remains; completed checkpoints are preserved."
            )

        label = (
            f"{projector} block {block:05d} "
            f"({block_source.stat().st_size:,} source bytes), "
            f"{number}/{len(blocks)}"
        )
        log("START " + label)
        (run / "current_stage.txt").write_text(
            "NATIVE_EXACT_RATIONAL_SERIES\n" + label + "\n"
        )

        started = time.monotonic()
        prefix = output / "arithmetic"
        dag = output / "arithmetic.dag.tsv"
        logfile = output / "native.log"
        partial_path = output / "partial.json"

        if partial_path.exists():
            partial = json.loads(partial_path.read_text())
            if (
                partial.get("source_sha256") != source_hash
                or partial.get("code_signature") != signature
            ):
                raise RuntimeError("stale partial checkpoint " + str(partial_path))
        else:
            partial = {
                "source_sha256": source_hash,
                "code_signature": signature,
                "steps": {},
                "files": {},
            }

        phase_seconds = args.block_seconds

        def phase(key, command, result):
            if key in partial["steps"]:
                record = partial["steps"][key]
                if not result.exists() or sha(result) != record["sha256"]:
                    raise RuntimeError("partial output changed: " + str(result))
                return

            phase_started = time.monotonic()
            invoke(command, logfile, phase_seconds)
            if not result.is_file() or result.stat().st_size == 0:
                raise RuntimeError("native phase wrote no output " + key)

            record = {
                "name": result.name,
                "sha256": sha(result),
                "seconds": time.monotonic() - phase_started,
            }
            partial["steps"][key] = record
            if key != "dag":
                partial["files"][key] = record
            write_json(partial_path, partial)

        phase(
            "dag",
            [str(source / "rational_dag"), str(block_source), str(prefix)],
            dag,
        )

        with dag.open() as handle:
            header = handle.readline().strip().split("\t")
            root_line = handle.readline().split()

        if (
            len(header) != 4
            or header[0] != "roots"
            or int(header[1]) != 1
            or len(root_line) != 1
            or int(header[3]) not in (0, 1)
        ):
            raise RuntimeError("unexpected DAG header/mass-degree certificate")

        for order in (4, 6):
            result = output / f"jet{order}.wl"
            phase(
                f"jet{order}",
                [str(source / "jet_dag"), str(dag), str(result), str(order)],
                result,
            )

        references = []
        beta_points = [
            (1, 1000),
            (1, 500),
            (1, 200),
            (1, 100),
            (1, 50),
            (1, 20),
            (1, 10),
            (1, 5),
            (3, 10),
        ]

        for numerator, denominator in beta_points:
            key = f"ref_{numerator}_{denominator}"
            result = output / f"{key}.wl"
            phase(
                key,
                [
                    str(source / "reference_dag"),
                    str(dag),
                    str(result),
                    str(numerator),
                    str(denominator),
                ],
                result,
            )

            check_key = key + "_check"
            if check_key not in partial["steps"]:
                check_started = time.monotonic()
                invoke(
                    [
                        str(source / "check_fixed_reference"),
                        str(block_source),
                        str(result),
                        str(numerator),
                        str(denominator),
                    ],
                    logfile,
                    phase_seconds,
                )
                partial["steps"][check_key] = {
                    "passed": True,
                    "modular_comparisons": 6,
                    "seconds": time.monotonic() - check_started,
                }
                write_json(partial_path, partial)

            references.append(
                {
                    "Numerator": numerator,
                    "Denominator": denominator,
                    "File": result.name,
                    "SHA256": sha(result),
                }
            )

        audit = {
            "method": VERSION,
            "passed": True,
            "projector": projector,
            "block": block,
            "nf_power": nf_power,
            "source_sha256": source_hash,
            "code_signature": signature,
            "source_bytes": block_source.stat().st_size,
            "mass_degree": int(header[3]),
            "files": partial["files"],
            "references": references,
            "modular_comparisons": 54,
            "seconds_this_invocation": time.monotonic() - started,
            "total_measured_seconds": sum(
                value.get("seconds", 0) for value in partial["steps"].values()
            ),
            "scope": (
                "Exact beta Laurent coefficients of the rational source divided by M2 "
                "through orders4 and6. Unexpanded references are separate exact rational "
                "functions of cos(theta), obtained by substituting each beta in the "
                "unchanged arithmetic DAG. Modular checks test those references; no "
                "fit/reconstruction from sampled numerical data is used for series "
                "coefficients. This is not full amplitude validation."
            ),
        }
        write_json(audit_path, audit)
        records.append(audit)

        if projector == "P1" and block == 1438:
            write_json(run / "results/failed_block_recovery_audit.json", audit)

        dag.unlink()
        write_json(
            run / "results/native_reduction_progress.json",
            {
                "reused_wolfram_blocks": reused,
                "ready_native_blocks": len(records),
                "remaining_native_blocks": len(blocks) - len(records),
                "last": audit,
            },
        )
        log(
            f"COMPLETED {projector} block {block:05d}: two exact truncated series; "
            "nine exact unexpanded references; 54 modular comparisons passed"
        )

    report = {
        "method": VERSION,
        "reused_wolfram_blocks": reused,
        "native_blocks": len(records),
        "all_required_blocks_ready": True,
        "code_signature": signature,
        "records": records,
        "scope": (
            "This completes rational preparation only. The Wolfram construction and "
            "numerical validation must still finish."
        ),
    }
    write_json(run / "results/native_reduction_report.json", report)
    (run / "NATIVE_RATIONAL_PREPARATION_COMPLETE.txt").write_text(
        "All required rational inputs are available; not an amplitude-completion marker.\n"
    )
    log("ALL REQUIRED RATIONAL BLOCKS ARE READY. Starting Wolfram stages next.")
    return 0


for caught_signal in (signal.SIGTERM, signal.SIGINT):
    signal.signal(caught_signal, stop)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        log("STOP: " + str(exc))
        raise SystemExit(1)
