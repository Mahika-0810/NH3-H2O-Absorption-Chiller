"""Run the chiller calculations, studies, checks, and plots."""

import subprocess
import sys
from pathlib import Path

folder = Path(__file__).resolve().parent
results = folder / "results"


def run_file(name, *args):
    """Run one project script and stop if it reports an error."""
    print(f"\nRunning {name}", *args)
    result = subprocess.run(
        [sys.executable, str(folder / name), *args],
        cwd=folder,
    )
    if result.returncode != 0:
        raise SystemExit(f"There was an error in {name}")


def organize_results():
    """Group generated files so the main results folder is easy to navigate."""
    import shutil

    plot_folder = results / "plots"
    detail_folder = results / "details"
    plot_folder.mkdir(exist_ok=True)
    detail_folder.mkdir(exist_ok=True)

    main_tables = {
        "results.json",
        "model_check_summary.csv",
        "study_operating_point.csv",
        "study_state_points.csv",
        "study_component_balances.csv",
        "study_temperature_sweep.csv",
    }

    for item in list(results.iterdir()):
        if not item.is_file():
            continue
        if item.suffix.lower() == ".png":
            destination = plot_folder / item.name
        elif item.suffix.lower() == ".csv" and item.name not in main_tables:
            destination = detail_folder / item.name
        else:
            continue
        if destination.exists():
            destination.unlink()
        shutil.move(str(item), str(destination))

    readme = results / "RESULTS_README.txt"
    readme.write_text(
        "MAIN FILES\n"
        "- study_operating_point.csv: summary of the chosen operating point\n"
        "- study_state_points.csv: states 1-10 and their calculated properties\n"
        "- study_component_balances.csv: mass and energy balance checks\n"
        "- study_temperature_sweep.csv: performance across condenser temperatures\n"
        "- model_check_summary.csv: project-level checks\n"
        "- results.json: complete structured run summary\n\n"
        "PLOTS\n"
        "Open the plots folder to view all generated figures.\n\n"
        "DETAILS\n"
        "The details folder contains additional CSVs from seasonal, condenser,\n"
        "air-handler, and sensitivity analyses. They are kept for traceability.\n",
        encoding="utf-8",
    )


def run_property_tests():
    """Run the equilibrium and thermodynamic-property regression tests."""
    print("\nRunning equilibrium property tests")
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "tests/test_equilibrium.py"],
        cwd=folder,
    )
    if result.returncode != 0:
        raise SystemExit("The equilibrium property tests failed")


def main():
    results.mkdir(exist_ok=True)
    files = [
        ("run_cycle.py",),
        ("study.py",),
        ("cycle_study.py",),
        ("study.py", "--sensitivity"),
        ("check_balances.py",),
    ]

    for command in files:
        run_file(*command)
    run_property_tests()

    try:
        import pandas as pd
        table = pd.read_csv(results / "model_check_summary.csv")
    except Exception as error:
        raise SystemExit(f"Could not read the check results: {error}")

    if "Result" not in table.columns or (table["Result"] == "FAIL").any():
        print(table)
        raise SystemExit("Some checks failed. See results/model_check_summary.csv")

    organize_results()
    print("\nDone! Main tables are in results/, plots are in results/plots/, and supporting tables are in results/details/.")


if __name__ == "__main__":
    main()
