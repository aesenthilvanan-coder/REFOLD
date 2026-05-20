FROM python:3.11-slim

WORKDIR /refold

RUN apt-get update && apt-get install -y \
    gcc g++ make wget git \
    libopenbabel-dev openbabel \
    dssp fpocket \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY requirements-ml.txt .
RUN pip install --no-cache-dir torch==2.3.0 --index-url https://download.pytorch.org/whl/cpu \
    && pip install --no-cache-dir fair-esm==2.0.0 torch-geometric==2.5.3

COPY . .
RUN pip install --no-cache-dir -e .

RUN mkdir -p data/raw data/processed data/results checkpoints logs

ENTRYPOINT ["python", "-m", "refold"]
CMD ["--help"]
