import importlib.util
from pathlib import Path


MODULE_PATH = (
    Path(__file__).resolve().parents[2]
    / "infrastructure"
    / "host_runner"
    / "host_runner.py"
)

spec = importlib.util.spec_from_file_location("ionmc_host_runner", MODULE_PATH)
host_runner = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(host_runner)


def test_bounded_tail_short_text():
    text, truncated, total = host_runner.bounded_tail("hello", 10)

    assert text == "hello"
    assert truncated is False
    assert total == 5


def test_bounded_tail_exact_limit():
    text, truncated, total = host_runner.bounded_tail("12345", 5)

    assert text == "12345"
    assert truncated is False
    assert total == 5


def test_bounded_tail_long_text_returns_tail():
    text, truncated, total = host_runner.bounded_tail("abcdefghij", 4)

    assert text == "ghij"
    assert truncated is True
    assert total == 10


def test_inline_output_limit_is_bounded():
    assert host_runner.MAX_INLINE_OUTPUT_CHARS > 0
    assert host_runner.MAX_INLINE_OUTPUT_CHARS <= 100000


def test_linked_worktree_git_metadata_is_bound_read_only(tmp_path):
    worktree = tmp_path / "worktrees" / "v2-001"
    worktree.mkdir(parents=True)
    common = tmp_path / "main" / ".git"
    common.mkdir(parents=True)
    assert host_runner.git_dir_binds(worktree, common) == [
        "--ro-bind",
        str(common.resolve()),
        str(common.resolve()),
    ]


def test_primary_checkout_git_metadata_needs_no_extra_bind(tmp_path):
    inside = tmp_path / ".git"
    inside.mkdir()
    assert host_runner.git_dir_binds(tmp_path, inside) == []
    assert host_runner.git_dir_binds(tmp_path, None) == []
