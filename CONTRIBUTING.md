# Contributing

Thank you for your interest in improving ISTSD-Profiler.

## Reporting issues

Please use GitHub Issues for bug reports, installation problems, and feature
requests. Include:

- the ISTSD-Profiler version or commit hash
- the operating system and environment manager used
- the command that was run
- the relevant error message or log excerpt
- a small input example, if possible

## Development setup

The recommended setup is the conda or mamba environment in `environment.yml`:

```bash
mamba env create -f environment.yml
conda activate istsd-profiler
```

Before opening a pull request, please run at least:

```bash
python -m py_compile ISTSD-Profiler.py
```

If your change affects pipeline output, also run the bundled example data or a
small equivalent test case and describe the result in the pull request.

## Pull requests

- Keep changes focused on one topic.
- Update `README.md` when command-line usage, dependencies, or outputs change.
- Do not commit generated runtime outputs, indexes, logs, or temporary files.
- Do not include unpublished third-party data unless redistribution is allowed.

