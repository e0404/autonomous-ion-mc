import json
from pathlib import Path

import pytest

from infrastructure.experiment_v2 import data, notify, reference
from infrastructure.experiment_v2.common import file_hash, inventory, safe_path
from infrastructure.experiment_v2.policy import direct_github_command


def test_runtime_hash_detects_change_and_external_symlink(tmp_path):
    f = tmp_path / "bin"
    f.write_bytes(b"version1")
    before = reference.runtime_inventory(tmp_path)
    f.write_bytes(b"version2")
    assert reference.runtime_inventory(tmp_path) != before
    (tmp_path / "escape").symlink_to("/etc/passwd")
    with pytest.raises(ValueError):
        reference.runtime_inventory(tmp_path)


@pytest.mark.parametrize("path", ["../secret", "/etc/passwd", "a/../../secret", ""])
def test_paths_reject_escape(tmp_path, path):
    with pytest.raises(ValueError):
        safe_path(tmp_path, path)


def test_input_symlink_rejected(tmp_path):
    (tmp_path / "link").symlink_to("/etc/passwd")
    with pytest.raises(ValueError):
        inventory(tmp_path)


def test_controlled_argv_has_no_shell_network_home_or_inherited_env(
    tmp_path, monkeypatch
):
    monkeypatch.setattr(reference.shutil, "which", lambda _: "/usr/bin/bwrap")
    config = {
        "engine": "topas",
        "executable": "/opt/reference/topas/bin/topas",
        "mounts": [{"source": "/runtime", "destination": "/opt/reference/topas"}],
    }
    argv = reference.command(
        config, {"input": "input.txt", "arguments": ["literal;no-shell"]}, tmp_path
    )
    assert "--unshare-all" in argv and "--clearenv" in argv
    assert argv[-3:] == [
        "/opt/reference/topas/bin/topas",
        "input.txt",
        "literal;no-shell",
    ]
    assert "/home" not in argv and "/var/run/docker.sock" not in argv
    assert "--share-net" not in argv


def test_permission_notification_is_redacted_and_not_approval(tmp_path):
    calls = []

    def sender(*args):
        calls.append(args)
        return {"status": "sent"}

    payload = {
        "hook_event_name": "PermissionRequest",
        "tool_input": {"command": "SECRET_TOKEN"},
        "session_id": "secret",
    }
    assert notify.handle(payload, sender=sender, state=tmp_path) == "sent"
    assert "SECRET_TOKEN" not in str(calls)
    events = (tmp_path / "events.jsonl").read_text()
    assert "SECRET_TOKEN" not in events and "secret" not in events
    assert "permission_prompt" in events and "notification_result" in events
    assert "intervention_resolved" not in events


def test_failed_notification_and_idle_prompt(tmp_path):
    assert (
        notify.handle(
            {"hook_event_name": "Notification", "notification_type": "idle_prompt"},
            state=tmp_path,
        )
        is None
    )
    result = notify.handle(
        {"hook_event_name": "Notification", "notification_type": "permission_prompt"},
        sender=lambda *args: {"status": "failed", "error": "SECRET_URL"},
        state=tmp_path,
    )
    assert result == "failed"
    assert "SECRET_URL" not in (tmp_path / "events.jsonl").read_text()


@pytest.mark.parametrize(
    "command",
    [
        "gh pr checks",
        "/usr/bin/gh pr checks",
        "env GH_PAGER=cat gh run list",
        "true && gh pr merge",
    ],
)
def test_direct_ci_guard(command):
    assert direct_github_command(command)


def test_data_use_retains_calibration_history(tmp_path):
    obj = tmp_path / "data/objects"
    obj.mkdir(parents=True)
    f = obj / "temp"
    f.write_bytes(b"fixture dataset")
    checksum = file_hash(f)
    f.rename(obj / checksum)
    data.assign_role(checksum, "calibration", "fit model", sha="a" * 40, state=tmp_path)
    data.assign_role(
        checksum,
        "evaluation",
        "correlated check, not held out",
        sha="a" * 40,
        state=tmp_path,
    )
    roles = [
        json.loads(p.read_text())["role"]
        for p in (tmp_path / "data/uses").glob("*.json")
    ]
    assert set(roles) == {"calibration", "evaluation"}
    (obj / checksum).write_bytes(b"corrupt")
    with pytest.raises(ValueError):
        data.assign_role(checksum, "evaluation", "why", sha="a" * 40, state=tmp_path)


@pytest.mark.parametrize(
    "url",
    [
        "file:///etc/passwd",
        "http://example.org/data",
        "https://user:pass@example.org/data",
        "https://localhost/data",
    ],
)
def test_data_rejects_nonpublic_or_authenticated_urls(url):
    with pytest.raises(ValueError):
        data.public_url(url)


def test_bounded_artifact_and_traversal(tmp_path):
    run = tmp_path / "references/REF-test"
    run.mkdir(parents=True)
    (run / "log").write_bytes(b"x" * 20000)
    result = reference.read_artifact("REF-test", "log", limit=100, state=tmp_path)
    assert result["bytes"] == 100 and result["total_bytes"] == 20000
    with pytest.raises(ValueError):
        reference.read_artifact("REF-test", "../../secret", state=tmp_path)


def test_process_failure_missing_outputs_and_timeout_are_retained(
    tmp_path, monkeypatch
):
    root = tmp_path / "repo"
    root.mkdir()
    case = root / "case"
    case.mkdir()
    (case / "input").write_text("fixture")
    envelope = {
        k: "fixture"
        for k in (
            "physics",
            "source",
            "geometry",
            "materials",
            "cuts",
            "scorers",
            "units",
            "rationale",
        )
    }
    envelope.update(
        input="input", histories=1, seeds=[1], expected_outputs=["out"], arguments=[]
    )
    (case / "case.json").write_text(json.dumps(envelope))
    monkeypatch.setattr(reference, "exact_state", lambda root: "a" * 40)
    monkeypatch.setattr(reference, "require_committed_inputs", lambda *args: None)
    monkeypatch.setattr(reference, "load_engine", lambda *args: {"engine": "topas"})
    monkeypatch.setattr(reference, "command", lambda *args: ["/usr/bin/false"])
    result = reference.run_reference(root, "case", "topas", state=tmp_path / "state")
    assert not result["succeeded"] and result["exit_code"] != 0
    assert result["missing_outputs"] == ["out"]
    monkeypatch.setattr(reference, "command", lambda *args: ["/usr/bin/sleep", "5"])
    result = reference.run_reference(
        root, "case", "topas", state=tmp_path / "state", timeout=1
    )
    assert result["timed_out"] and not result["succeeded"]
    monkeypatch.setattr(reference, "command", lambda *args: ["/usr/bin/true"])
    result = reference.run_reference(root, "case", "topas", state=tmp_path / "state")
    assert result["exit_code"] == 0 and not result["succeeded"]
    assert len(list((tmp_path / "state/references").glob("*/result.json"))) == 3


def test_native_internal_output_alias_is_allowed_but_escape_is_not(tmp_path):
    (tmp_path / "score").mkdir()
    (tmp_path / "score/dose").write_bytes(b"dose data")
    (tmp_path / "dose").symlink_to("score/dose")
    assert inventory(tmp_path, allow_internal_links=True)["dose"]["bytes"] == 9
    with pytest.raises(ValueError):
        inventory(tmp_path)
    (tmp_path / "outside").symlink_to("/etc/passwd")
    with pytest.raises(ValueError):
        inventory(tmp_path, allow_internal_links=True)


def test_notification_transport_to_local_receiver(monkeypatch):
    import threading
    from http.server import BaseHTTPRequestHandler, HTTPServer

    from infrastructure.interventions.request_intervention import send_notification

    received = []

    class Receiver(BaseHTTPRequestHandler):
        def do_POST(self):
            received.append(self.rfile.read(int(self.headers["Content-Length"])))
            self.send_response(200)
            self.end_headers()

        def log_message(self, *args):
            pass

    monkeypatch.setenv("no_proxy", "127.0.0.1,localhost")
    server = HTTPServer(("127.0.0.1", 0), Receiver)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        monkeypatch.setenv(
            "IONMC_NTFY_URL", f"http://127.0.0.1:{server.server_port}/fixture"
        )
        result = send_notification(
            "TEST", "setup", "permission-approval", "Generic fixture notification"
        )
        assert result["status"] == "sent" and b"TEST" in received[0]
    finally:
        server.shutdown()
        server.server_close()
        thread.join()


def test_v2_intervention_requires_real_investigation_record():
    from infrastructure.experiment_v2.intervention import validate_context

    with pytest.raises(ValueError):
        validate_context(None)
    with pytest.raises(ValueError):
        validate_context({"blocked_claim": "fragmentation"})
    value = {
        "blocked_claim": "fragmentation evidence",
        "budget": "bounded fixture",
        "alternatives": "no accessible substitute",
        "attempts": [
            {
                "approach": "public source",
                "outcome": "unavailable",
                "artifact": "fixture-log",
            }
        ],
    }
    assert validate_context(value) == value


def test_ignored_reference_inputs_cannot_claim_committed_provenance(monkeypatch):
    monkeypatch.setattr(reference, "git", lambda *args: "case/input.txt")
    reference.require_committed_inputs(Path("/repo"), "case", {"input.txt": {}})
    with pytest.raises(ValueError):
        reference.require_committed_inputs(
            Path("/repo"), "case", {"input.txt": {}, "ignored.txt": {}}
        )
