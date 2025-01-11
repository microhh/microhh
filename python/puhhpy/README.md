# PUHHPY
### Python MicroHH open-boundary Pre-processor 4 You.

Yes, worst name ever. Suggestions are welcome.

### Usage
Either add the `puhhpy` package location to your `PYTHONPATH`:

    export PYTHONPATH="${PYTHONPATH}:/path/to/microhh/python/puhhpy"

Or specify the path using `sys`, before importing `puhppy`:

    import sys
    sys.path.append('/path/to/microhh/python/puhhpy')

Now `puhhpy` should be available as an import, e.g.

    from puhhpy.spatial import Domain
    from puhhpy.spatial import Projection
