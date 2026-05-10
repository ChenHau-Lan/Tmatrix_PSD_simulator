"""Small stdlib web server for the local T-matrix PSD explorer."""

from __future__ import annotations

from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
import os
from pathlib import Path
import traceback
from urllib.parse import urlparse

from .local_tmatrix import HYDROMETEORS, compute


ROOT = Path(__file__).resolve().parents[2]
FRONTEND = ROOT / "frontend"


class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(FRONTEND), **kwargs)

    def log_message(self, fmt: str, *args) -> None:
        print(fmt % args)

    def _send_json(self, status: int, payload: dict) -> None:
        body = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:
        path = urlparse(self.path).path
        if path == "/api/health":
            self._send_json(
                200,
                {
                    "status": "ok",
                    "engine": "local Fortran py_Tmatrix_Mueller",
                    "hydrometeors": {
                        key: {"label": value.label, "verified": value.verified}
                        for key, value in HYDROMETEORS.items()
                    },
                },
            )
            return
        if path == "/":
            self.path = "/index.html"
        return super().do_GET()

    def do_POST(self) -> None:
        path = urlparse(self.path).path
        if path != "/api/compute":
            self._send_json(404, {"error": "not found"})
            return

        try:
            length = int(self.headers.get("Content-Length", "0"))
            req = json.loads(self.rfile.read(length) or b"{}")
            result = compute(
                dm=float(req.get("dm", 2.0)),
                log_nw=float(req.get("log_nw", 3.5)),
                mu=float(req.get("mu", 0.0)),
                band=str(req.get("band", "S")),
                hydro_key=str(req.get("hydrometeor", "rain")),
                canting_std_deg=float(req.get("canting_std_deg", 0.0)),
                d_step=float(req.get("d_step", 0.1)),
            )
            self._send_json(200, result)
        except Exception as exc:
            traceback.print_exc()
            self._send_json(500, {"error": str(exc)})


def run(host: str | None = None, port: int | None = None) -> None:
    host = host or os.environ.get("HOST", "127.0.0.1")
    port = port or int(os.environ.get("PORT", "7860"))
    server = ThreadingHTTPServer((host, port), Handler)
    print(f"Serving local T-matrix PSD explorer at http://{host}:{port}")
    server.serve_forever()


if __name__ == "__main__":
    run()
