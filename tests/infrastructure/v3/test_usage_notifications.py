import json
import os
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path

import pytest

from infrastructure.experiment_v3 import notify
from infrastructure.experiment_v3 import usage_notifications as usage


def payload(session="fixture-session", message="You've hit your weekly limit SECRET"):
    return {
        "hook_event_name": "StopFailure",
        "error": "rate_limit",
        "session_id": session,
        "last_assistant_message": message,
        "error_details": "SECRET_URL",
        "transcript_path": "/secret/transcript",
    }


def test_alert_cooldown_sessions_and_redaction(tmp_path, monkeypatch):
    now = [10000]
    monkeypatch.setattr(usage.time, "time", lambda: now[0])
    sent = []

    def sender(*args):
        sent.append(args)
        return {"status": "sent"}

    assert notify.handle(payload(), sender=sender, state=tmp_path) == "sent"
    assert notify.handle(payload(), sender=sender, state=tmp_path) == "suppressed"
    assert notify.handle(payload("another"), sender=sender, state=tmp_path) == "sent"
    now[0] += usage.SUCCESS_COOLDOWN
    assert notify.handle(payload(), sender=sender, state=tmp_path) == "sent"
    assert len(sent) == 3
    records = (tmp_path / "events.jsonl").read_text()
    for secret in ["SECRET", "fixture-session", "/secret/transcript"]:
        assert secret not in str(sent) and secret not in records
    assert "experiment-v3" in sent[0][-1] and "weekly" in sent[0][-1]
    assert records.count('"event": "runtime_limit"') == 4
    assert '"status": "suppressed"' in records
    assert "intervention_requested" not in records


@pytest.mark.parametrize("failure", ["failed", "not-configured", "raises"])
def test_failed_delivery_retries_after_short_cooldown(tmp_path, monkeypatch, failure):
    now = [10000]
    monkeypatch.setattr(usage.time, "time", lambda: now[0])
    attempts = []

    def sender(*args):
        attempts.append(args)
        if len(attempts) > 1:
            return {"status": "sent"}
        if failure == "raises":
            raise RuntimeError("SECRET transport exception")
        return {"status": failure}

    first = notify.handle(payload(), sender=sender, state=tmp_path)
    assert first == ("failed" if failure == "raises" else failure)
    assert notify.handle(payload(), sender=sender, state=tmp_path) == "suppressed"
    now[0] += usage.FAILURE_COOLDOWN
    assert notify.handle(payload(), sender=sender, state=tmp_path) == "sent"
    assert len(attempts) == 2
    assert "SECRET" not in (tmp_path / "events.jsonl").read_text()


def test_concurrent_hooks_send_only_one_alert(tmp_path):
    entered, release = threading.Event(), threading.Event()
    calls = []

    def sender(*args):
        calls.append(args)
        entered.set()
        assert release.wait(5)
        return {"status": "sent"}

    with ThreadPoolExecutor(max_workers=2) as pool:
        first = pool.submit(notify.handle, payload(), sender=sender, state=tmp_path)
        try:
            assert entered.wait(5)
            assert (
                notify.handle(payload(), sender=sender, state=tmp_path) == "suppressed"
            )
        finally:
            release.set()
        assert first.result(timeout=5) == "sent"
    assert len(calls) == 1


def test_other_failures_do_not_trigger_usage_alert(tmp_path):
    def sender(*args):
        pytest.fail("Unexpected notification")

    for error in ["server_error", "authentication_failed", "max_output_tokens"]:
        value = payload()
        value["error"] = error
        assert notify.handle(value, sender=sender, state=tmp_path) is None
    assert not (tmp_path / "events.jsonl").exists()


@pytest.mark.parametrize(
    "text,expected",
    [
        ("You've hit your session limit", "session"),
        ("You've hit your Opus limit", "opus"),
        ("You've hit your Fable limit", "fable"),
        ("You've hit your Sonnet limit", "sonnet"),
        ("429 Too Many Requests", "unspecified"),
    ],
)
def test_category_is_best_effort_fixed_label(text, expected):
    assert usage.limit_kind(payload(message=text)) == expected


def test_cooldown_state_does_not_follow_symlink(tmp_path):
    directory = tmp_path / "notification-cooldowns"
    target = tmp_path / "outside"
    target.mkdir()
    directory.symlink_to(target, target_is_directory=True)
    with pytest.raises(OSError):
        notify.handle(
            payload(), sender=lambda *a: pytest.fail("Unexpected send"), state=tmp_path
        )
    assert not list(target.iterdir())


def test_configured_hook_subprocess_delivers_once_to_local_receiver(tmp_path):
    received = []

    class Receiver(BaseHTTPRequestHandler):
        def do_POST(self):
            received.append(self.rfile.read(int(self.headers["Content-Length"])))
            self.send_response(200)
            self.end_headers()

        def log_message(self, *args):
            pass

    server = HTTPServer(("127.0.0.1", 0), Receiver)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    root = Path(__file__).resolve().parents[3]
    settings = json.loads((root / ".claude/settings.json").read_text())
    hooks = [
        x for x in settings["hooks"]["StopFailure"] if x.get("matcher") == "rate_limit"
    ]
    assert len(hooks) == 1
    assert (
        hooks[0]["hooks"][0]["command"]
        == "python3 -m infrastructure.experiment_v3.notify"
    )
    template = json.loads((root / "experiment/v3/claude-settings.json").read_text())
    assert template["hooks"] == settings["hooks"]
    env = dict(
        os.environ,
        IONMC_NTFY_URL=f"http://127.0.0.1:{server.server_port}/test",
        IONMC_V3_EVENT_DIR=str(tmp_path),
        no_proxy="127.0.0.1,localhost",
    )
    try:
        for _ in range(2):
            result = subprocess.run(
                [sys.executable, "-m", "infrastructure.experiment_v3.notify"],
                input=json.dumps(payload()),
                text=True,
                capture_output=True,
                env=env,
                cwd=root,
                timeout=10,
                check=True,
            )
            assert not result.stdout and not result.stderr
        assert len(received) == 1
        assert b"usage-limit" in received[0] and b"experiment-v3" in received[0]
        assert b"SECRET" not in received[0]
    finally:
        server.shutdown()
        server.server_close()
        thread.join()
