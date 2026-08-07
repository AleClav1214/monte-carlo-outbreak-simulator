"""
Generate docs/evidence_tables.md from data/parameters/pathogens.yaml.

Deliberately a generated artifact, not hand-maintained: the YAML file is the
single source of truth (see data/schemas.py module docstring), and this
script is how it gets rendered into the human-readable "structured evidence
table" the project brief asks for. Re-run this after any edit to
pathogens.yaml rather than hand-editing docs/evidence_tables.md.

Run with: python scripts/generate_evidence_tables.py
"""

from __future__ import annotations

from pathlib import Path

from outbreak_simulator.data import load_all

OUTPUT_PATH = Path(__file__).resolve().parent.parent / "docs" / "evidence_tables.md"

QUALITY_LABEL = {
    "meta_analysis": "Meta-analysis",
    "multi_study": "Multi-study",
    "single_study": "Single study",
    "expert_consensus": "Expert consensus",
    "modeled_estimate": "Modeled estimate",
    "low_confidence": "⚠️ LOW CONFIDENCE",
}


def render() -> str:
    pathogens, scenarios = load_all()
    lines = [
        "# Evidence Tables",
        "",
        "**Generated from `data/parameters/pathogens.yaml` — do not hand-edit this file.**",
        "Re-run `python scripts/generate_evidence_tables.py` after changing the YAML.",
        "",
        "See `data/parameters/pathogens.yaml`'s file header for the full literature-review",
        "methodology statement (targeted search, not a PRISMA systematic review).",
        "",
    ]

    for pid in sorted(pathogens.keys()):
        p = pathogens[pid]
        lines.append(f"## {p.display_name}")
        lines.append("")
        lines.append(f"*{p.pathogen_class}. Last reviewed: {p.last_reviewed}.*")
        lines.append("")
        lines.append("| Parameter | Value | Range (approx. 95% interval) | Unit | Evidence quality | Source |")
        lines.append("|---|---|---|---|---|---|")
        for param_name in sorted(p.parameters.keys()):
            param = p.parameters[param_name]
            rng = f"[{param.low:.3g}, {param.high:.3g}]" if param.low is not None else "—"
            quality = QUALITY_LABEL.get(param.evidence_quality.value, param.evidence_quality.value)
            source_short = param.source if len(param.source) < 90 else param.source[:87] + "..."
            lines.append(f"| {param.display_name} | {param.point_estimate:.3g} | {rng} | {param.unit} | {quality} | {source_short} |")
        lines.append("")
        lines.append("**Justifications and notes:**")
        lines.append("")
        for param_name in sorted(p.parameters.keys()):
            param = p.parameters[param_name]
            lines.append(f"- **{param.display_name}:** {param.justification.strip()}")
            if param.notes:
                lines.append(f"  *Notes:* {param.notes.strip()}")
        lines.append("")

    lines.append("## Scenario-level setting parameters")
    lines.append("")
    lines.append("| Scenario | Pathogen | Population | Initial cases | Contact multiplier | Real-world benchmark(s) |")
    lines.append("|---|---|---|---|---|---|")
    for sid in sorted(scenarios.keys()):
        s = scenarios[sid]
        benchmarks = "; ".join(f"{o.attack_rate:.1%}" for o in s.observed_outcomes if o.attack_rate is not None) or "none"
        lines.append(
            f"| {s.display_name} | {s.pathogen_id} | {s.population.population_size} | "
            f"{s.population.initial_cases} | {s.population.contact_rate_multiplier}x | {benchmarks} |"
        )
    lines.append("")

    return "\n".join(lines)


if __name__ == "__main__":
    content = render()
    OUTPUT_PATH.write_text(content)
    print(f"Wrote {OUTPUT_PATH} ({len(content)} characters)")
