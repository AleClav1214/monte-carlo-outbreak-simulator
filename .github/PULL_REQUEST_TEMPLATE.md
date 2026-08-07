## What does this PR do?



## Checklist

- [ ] `pytest tests/` passes
- [ ] `ruff check src/ tests/` passes
- [ ] If I changed `data/parameters/pathogens.yaml` or `data/scenarios/scenarios.yaml`,
      I ran `python scripts/generate_evidence_tables.py` and committed the
      regenerated `docs/evidence_tables.md`
- [ ] If I added a parameter, it has a real `source`, a specific `justification`,
      and an honest `evidence_quality` (see CONTRIBUTING.md)
- [ ] If I added a scenario, it has at least one `observed_outcomes` entry
      (or, if none exists in the literature, I said so explicitly in `limitations`)
- [ ] If I fixed a bug, I added a test that would have failed before the fix

## Related issue(s)

Closes #

## Notes for reviewers

Anything that needs particular scrutiny (a new statistical method, a
parameter you're less confident about, a design tradeoff you want feedback on).
