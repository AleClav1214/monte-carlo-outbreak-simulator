---
name: Bug report
about: Something in the code doesn't work as documented
title: "[BUG] "
labels: bug
---

**Describe the bug**
A clear description of what's wrong.

**To reproduce**
```python
# Minimal code that reproduces the issue
```

**Expected behavior**
What you expected to happen.

**Environment**
- Python version:
- outbreak-simulator version / commit:
- OS:
- Output of `pip list | grep -E "numpy|scipy|pydantic|pandas"`:

**Additional context**
If this involves a specific scenario or pathogen's output looking wrong,
please include the `scenario_id`/`pathogen_id` and the seed you used —
this project's Monte Carlo runs are meant to be exactly reproducible
(see docs/reproducibility.md), so a seed lets us reproduce your exact result.
