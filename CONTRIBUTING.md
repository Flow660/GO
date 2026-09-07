# Contributing

Thank you for your interest in this project. It accompanies a scientific
publication, so the priority is reproducibility of the published analysis.

## Reporting issues

Please open a GitHub issue describing:

- what you expected to happen and what happened instead;
- the versions of Python, Streamlit and pandas you are using;
- a minimal input workbook (with any confidential data removed) that reproduces
  the problem.

## Pull requests

1. Fork the repository and create a branch for your change.
2. Keep changes focused; behavioural changes to the processing logic should be
   explained in the pull request description, because they may affect the
   published results.
3. Add or update the tests under `tests/` and make sure they pass:

   ```bash
   pip install -r requirements.txt pytest
   pytest
   ```

4. Follow the existing code style (PEP 8, four-space indentation, English
   comments and docstrings).

## License of contributions

By contributing you agree that your contributions are licensed under the
BSD 3-Clause License that covers this project.
