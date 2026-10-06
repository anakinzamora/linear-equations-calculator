# System of Linear Equations Calculator: web page + GNU Octave running the MATLAB solver code.
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    HOME=/tmp \
    PORT=7860 \
    PYTHONUNBUFFERED=1

RUN apt-get update \
 && apt-get install -y --no-install-recommends octave python3 \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY server.py ./
COPY matlab ./matlab
COPY static ./static

# Fail the build early if Octave cannot load the solver code.
RUN octave-cli --norc --quiet --no-history --eval "source('/app/matlab/solver_core.m'); disp(fnum(2.00005, 4))"

EXPOSE 7860
CMD ["python3", "server.py"]
