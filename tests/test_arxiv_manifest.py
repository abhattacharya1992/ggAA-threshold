"""Check ancillary selection without producing the submission archive."""
import importlib.util
import pathlib
import tempfile
import unittest
from unittest.mock import patch


SCRIPT = pathlib.Path(__file__).resolve().parents[1] / 'scripts' / 'make_arxiv.py'
SPEC = importlib.util.spec_from_file_location('make_arxiv', SCRIPT)
make_arxiv = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(make_arxiv)


class AncillaryManifest(unittest.TestCase):
    def test_excludes_manuscript_snapshots_and_keeps_scientific_data(self):
        kept = {'data/born/original/BornSqr.log', 'src/python/model.py', 'results/table.json'}
        excluded = {
            'manuscript/main.pdf',
            'manuscript_arXiv_v1/main.pdf',
            'manuscript_arXiv_v1/figures/result.pdf',
            'history/pre_v2_snapshot/main.pdf',
            'history/pre_v2_snapshot/main.bbl',
            'build/previous.tar.gz',
            'results/exact_amplitude_cache.json',
            'run.log',
            'notes.tex',
            'package.zip',
            'data/archive.zip',
        }
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            for name in kept | excluded:
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text('fixture\n')
            with patch.object(make_arxiv, 'ROOT', root):
                selected = {path.as_posix() for path in make_arxiv.ancillary_files()}
        self.assertEqual(selected, kept)


if __name__ == '__main__':
    unittest.main()
