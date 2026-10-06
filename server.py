"""Web server for the System of Linear Equations Calculator.

Serves the page in static/ and answers POST /api/solve by running the MATLAB
solver code (matlab/solver_core.m) in GNU Octave. The browser does no math:
every number on the page comes back from the .m code.

Run locally:  python3 server.py      then open http://localhost:7860
"""
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
import threading
import time
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

APP_DIR = os.path.dirname(os.path.abspath(__file__))
STATIC_DIR = os.path.join(APP_DIR, "static")
MATLAB_DIR = os.path.join(APP_DIR, "matlab")
OCTAVE = shutil.which("octave-cli") or shutil.which("octave") or "octave-cli"
OCTAVE_ARGS = ["--norc", "--quiet", "--no-history", "--no-window-system"]
PORT = int(os.environ.get("PORT", "7860"))
TIMEOUT_S = int(os.environ.get("SOLVE_TIMEOUT", "60"))
MAX_BODY = 64 * 1024
# Free hosting has little memory, so run at most two Octave processes at once.
SLOTS = threading.BoundedSemaphore(int(os.environ.get("MAX_OCTAVE", "2")))

METHODS = {"doolittle", "crout", "cholesky", "cramer", "gaussjordan", "seidel", "jacobi"}
CRITERIA = {"abs", "rel", "repeat"}
NAME_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]{0,11}$")


def octave_version():
    try:
        out = subprocess.run([OCTAVE, "--version"], capture_output=True, text=True, timeout=20).stdout
        m = re.search(r"version\s+([\d.]+)", out)
        return f"GNU Octave {m.group(1)}" if m else "GNU Octave"
    except Exception:  # noqa: BLE001
        return "GNU Octave (not found)"


ENGINE = octave_version()


class BadRequest(ValueError):
    pass


def number(v, what):
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
        raise BadRequest(f"{what} must be a finite number.")
    return float(v)


def validate(req):
    if not isinstance(req, dict):
        raise BadRequest("The request must be a JSON object.")
    method = req.get("method")
    if method not in METHODS:
        raise BadRequest("Unknown method.")
    A, b = req.get("A"), req.get("b")
    if not isinstance(A, list) or not 3 <= len(A) <= 10:
        raise BadRequest("Enter between 3 and 10 equations.")
    n = len(A)
    if any(not isinstance(r, list) or len(r) != n for r in A):
        raise BadRequest("The coefficient matrix must be square.")
    if not isinstance(b, list) or len(b) != n:
        raise BadRequest(f"There must be {n} constants.")
    A = [[number(v, "Every coefficient") for v in row] for row in A]
    b = [number(v, "Every constant") for v in b]
    names = req.get("names")
    if not isinstance(names, list) or len(names) != n or not all(isinstance(s, str) and NAME_RE.match(s) for s in names):
        raise BadRequest(f"Give {n} names for the unknowns (letters, digits or _, starting with a letter).")
    if len(set(names)) != n:
        raise BadRequest("The names of the unknowns must all be different.")
    dp = req.get("dp", 4)
    if not isinstance(dp, int) or isinstance(dp, bool) or not 0 <= dp <= 10:
        raise BadRequest("Decimal places must be a whole number from 0 to 10.")
    clean = {"method": method, "A": A, "b": b, "names": names, "dp": dp}
    if method in ("seidel", "jacobi"):
        x0 = req.get("x0", [0] * n)
        if not isinstance(x0, list) or len(x0) != n:
            raise BadRequest(f"The initial guess needs {n} numbers.")
        crit = req.get("crit", "abs")
        if crit not in CRITERIA:
            raise BadRequest("Unknown stopping criterion.")
        tol = number(req.get("tol", 1e-4), "The tolerance")
        if tol <= 0:
            raise BadRequest("The tolerance must be greater than 0.")
        max_it = req.get("maxIt", 100)
        if not isinstance(max_it, int) or isinstance(max_it, bool) or not 1 <= max_it <= 5000:
            raise BadRequest("Maximum iterations must be a whole number from 1 to 5000.")
        clean.update(x0=[number(v, "Every initial guess") for v in x0], crit=crit, tol=tol,
                     maxIt=max_it, reorder=bool(req.get("reorder", True)))
    return clean


def run_octave(req):
    with tempfile.TemporaryDirectory(prefix="lss-") as wd:
        with open(os.path.join(wd, "request.json"), "w") as f:
            json.dump(req, f)
        started = time.time()
        with SLOTS:
            proc = subprocess.run(
                [OCTAVE, *OCTAVE_ARGS, os.path.join(MATLAB_DIR, "web_solve.m"), wd],
                capture_output=True, text=True, timeout=TIMEOUT_S, cwd=wd,
            )
        ms = int((time.time() - started) * 1000)
        result_path = os.path.join(wd, "result.json")
        if not os.path.exists(result_path):
            detail = (proc.stderr or proc.stdout or "").strip().splitlines()[-3:]
            raise RuntimeError("The MATLAB code did not finish. " + " ".join(detail))
        with open(result_path) as f:
            result = json.load(f)
        with open(os.path.join(wd, "output.txt")) as f:
            text = f.read()
    if isinstance(result.get("x"), str):
        result["x"] = [result["x"]]
    result.update(text=text, engine=ENGINE, ms=ms)
    return result


class Handler(SimpleHTTPRequestHandler):
    server_version = "LinearSystemsCalculator/1.0"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=STATIC_DIR, **kwargs)

    def send_json(self, code, payload):
        body = json.dumps(payload).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        if self.path.split("?")[0] == "/api/health":
            return self.send_json(200, {"ok": True, "engine": ENGINE})
        return super().do_GET()

    def do_POST(self):
        if self.path.split("?")[0] != "/api/solve":
            return self.send_json(404, {"ok": False, "error": "Not found."})
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if length <= 0 or length > MAX_BODY:
                raise BadRequest("The request is empty or too large.")
            req = validate(json.loads(self.rfile.read(length)))
            return self.send_json(200, run_octave(req))
        except (BadRequest, json.JSONDecodeError) as e:
            return self.send_json(400, {"ok": False, "status": "invalid", "error": str(e)})
        except subprocess.TimeoutExpired:
            return self.send_json(504, {"ok": False, "status": "timeout",
                                        "error": f"The MATLAB code took longer than {TIMEOUT_S} s. Try fewer iterations."})
        except Exception as e:  # noqa: BLE001
            return self.send_json(500, {"ok": False, "status": "server", "error": str(e)})

    def end_headers(self):
        self.send_header("X-Content-Type-Options", "nosniff")
        super().end_headers()

    def log_message(self, fmt, *args):
        print("%s - %s" % (self.address_string(), fmt % args), flush=True)


if __name__ == "__main__":
    print(f"Serving on port {PORT} using {ENGINE} ({OCTAVE})", flush=True)
    ThreadingHTTPServer(("0.0.0.0", PORT), Handler).serve_forever()
