"""Run figure renderers that use the bundled compact source data.

The default ``standard`` profile covers the standard Python panels. Manifold and
persistent-homology panels are provided through the ``topology`` profile.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class Task:
    name: str
    script: str
    profiles: tuple[str, ...]
    description: str
    environment: str
    arguments: tuple[str, ...] = ()
    depends_on: tuple[str, ...] = ()
    optional: bool = False

    @property
    def path(self) -> Path:
        return REPO_ROOT / self.script


TASKS = (
    Task(
        "figure_2i",
        "workflows/encoding_examples/render_figure_2i.py",
        ("standard", "all"),
        "Representative encoding-fit traces",
        "standard",
    ),
    Task(
        "figure_2d",
        "workflows/example_activity/render_figure_2d.py",
        ("standard", "all"),
        "Representative neural, platform, and right-toe traces",
        "standard",
    ),
    Task(
        "figure_2c",
        "workflows/neuron_density/render_figure_2c.py",
        ("standard", "all"),
        "Recorded-neuron locations and density",
        "standard",
    ),
    Task(
        "figure_3a_e",
        "workflows/decoding/render_figure_3a_e.py",
        ("standard", "all"),
        "Regional population-decoding panels",
        "standard",
    ),
    Task(
        "extended_data_7c",
        "workflows/encoding_performance/render_extended_data_figure_7c.py",
        ("standard", "all"),
        "Encoding-performance distributions",
        "standard",
    ),
    Task(
        "extended_data_9",
        "workflows/sex_indicator/render_extended_data_figure_9.py",
        ("standard", "all"),
        "Sex and indicator sensitivity panels",
        "standard",
    ),
    Task(
        "extended_data_10_neuron_count",
        "workflows/decoding/build_paper_ready_assets.py",
        ("standard", "all"),
        "Neuron-count sensitivity panels and saved statistics",
        "standard",
    ),
    Task(
        "extended_data_10_empirical_null",
        "workflows/decoding/recenter_empirical_null_figures.py",
        ("standard", "all"),
        "Recentered shuffled-target panels",
        "standard",
        depends_on=("extended_data_10_neuron_count",),
    ),
    Task(
        "figure_3h_and_extended_data_10k",
        "workflows/frontoparietal_coordination/step5_paper_ready_figures.py",
        ("topology", "all"),
        "Frontal-parietal error coordination panels",
        "topology",
    ),
    Task(
        "figure_3f_g_i",
        "workflows/frontoparietal_coordination/render_figure_3fgi.py",
        ("standard", "all"),
        "Frontal and parietal confusion matrices and joint decoding-error excess",
        "standard",
    ),
    Task(
        "figure_3l",
        "workflows/gcarp/render_figure_3l.py",
        ("standard", "all"),
        "GCARP contributions for movement, tilt, and posture",
        "standard",
    ),
    Task(
        "figures_4_6",
        "workflows/manifold_topology/render_figures_4_6.py",
        ("topology", "all"),
        "Computational panels of manuscript Figures 4, 5, and 6",
        "topology",
    ),
    Task(
        "extended_data_11cd",
        "workflows/manifold_topology/render_extended_data_figure_11_cd.py",
        ("topology", "all"),
        "Embedding-dimension controls rendered from bundled summary arrays",
        "topology",
        optional=True,
    ),
)


def selected_tasks(profile: str, requested: set[str], skipped: set[str]) -> list[Task]:
    by_name = {task.name: task for task in TASKS}
    unknown = (requested | skipped) - set(by_name)
    if unknown:
        raise ValueError(f"Unknown task name(s): {', '.join(sorted(unknown))}")

    initial = [
        task
        for task in TASKS
        if (not requested and profile in task.profiles) or task.name in requested
    ]
    names = {task.name for task in initial}
    changed = True
    while changed:
        changed = False
        for task in tuple(initial):
            for dependency in task.depends_on:
                if dependency not in names:
                    initial.insert(0, by_name[dependency])
                    names.add(dependency)
                    changed = True

    skipped_with_dependents = set(skipped)
    changed = True
    while changed:
        changed = False
        for task in initial:
            if any(dep in skipped_with_dependents for dep in task.depends_on):
                if task.name not in skipped_with_dependents:
                    skipped_with_dependents.add(task.name)
                    changed = True
    return [task for task in initial if task.name not in skipped_with_dependents]


def command_for(task: Task, interpreter: str | None = None) -> list[str]:
    arguments = [str(REPO_ROOT / arg) if arg.startswith("outputs/") else arg for arg in task.arguments]
    return [interpreter or sys.executable, str(task.path), *arguments]


def matplotlib_version(interpreter: str) -> str | None:
    result = subprocess.run(
        [interpreter, "-c", "import matplotlib; print(matplotlib.__version__)"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    return result.stdout.strip() if result.returncode == 0 else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=("standard", "topology", "all"), default="standard")
    parser.add_argument("--task", action="append", default=[], help="Run only this named task. Repeatable.")
    parser.add_argument("--skip", action="append", default=[], help="Skip this task and its dependents. Repeatable.")
    parser.add_argument("--list", action="store_true", help="List tasks without running them.")
    parser.add_argument("--continue-on-error", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--standard-python", default=sys.executable)
    parser.add_argument("--topology-python", default=sys.executable)
    parser.add_argument(
        "--ignore-version",
        action="store_true",
        help="Run despite a Matplotlib version mismatch. Figure appearance may vary.",
    )
    args = parser.parse_args()

    try:
        tasks = selected_tasks(args.profile, set(args.task), set(args.skip))
    except ValueError as exc:
        parser.error(str(exc))

    if args.list:
        for task in TASKS:
            profiles = ",".join(task.profiles)
            availability = "available" if task.path.is_file() else "missing"
            print(
                f"{task.name:44} [{profiles:13}] {task.environment:8} "
                f"{availability:9} {task.description}"
            )
        return 0

    env = os.environ.copy()
    source_path = str(REPO_ROOT / "src")
    env["PYTHONPATH"] = source_path + os.pathsep + env.get("PYTHONPATH", "")
    failures: list[str] = []
    completed: set[str] = set()
    interpreters = {"standard": args.standard_python, "topology": args.topology_python}
    expected_versions = {"standard": "3.8.4", "topology": "3.10.5"}
    checked_environments: dict[tuple[str, str], bool] = {}

    for task in tasks:
        if any(dep not in completed for dep in task.depends_on):
            print(f"SKIP {task.name}: a dependency did not complete")
            failures.append(task.name)
            continue
        if not task.path.is_file():
            label = "optional task" if task.optional else "required task"
            print(f"SKIP {task.name}: {label} script is absent: {task.script}")
            if not task.optional:
                failures.append(task.name)
            continue
        interpreter = interpreters[task.environment]
        environment_key = (task.environment, interpreter)
        if environment_key not in checked_environments and not args.dry_run:
            observed = matplotlib_version(interpreter)
            expected = expected_versions[task.environment]
            version_ok = observed == expected
            checked_environments[environment_key] = version_ok
            if not version_ok:
                message = (
                    f"{task.environment} profile requires Matplotlib {expected}, "
                    f"but {interpreter} reports {observed or 'unavailable'}"
                )
                if args.ignore_version:
                    print(f"NOTE: {message}")
                else:
                    print(f"FAIL {task.name}: {message}")
                    failures.append(task.name)
                    if not args.continue_on_error:
                        break
                    continue
        elif environment_key in checked_environments and not checked_environments[environment_key] and not args.ignore_version:
            print(f"FAIL {task.name}: incompatible {task.environment} environment")
            failures.append(task.name)
            if not args.continue_on_error:
                break
            continue
        command = command_for(task, interpreter)
        print(f"RUN  {task.name}: {' '.join(command)}")
        if args.dry_run:
            completed.add(task.name)
            continue
        result = subprocess.run(command, cwd=REPO_ROOT, env=env, check=False)
        if result.returncode == 0:
            completed.add(task.name)
            print(f"PASS {task.name}")
        else:
            failures.append(task.name)
            print(f"FAIL {task.name}: exit code {result.returncode}")
            if not args.continue_on_error:
                break

    if failures:
        print("Failed or blocked tasks: " + ", ".join(failures))
        return 1
    verb = "Planned" if args.dry_run else "Completed"
    print(f"{verb} {len(completed)} task(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
