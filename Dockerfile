FROM python:3.9-slim

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV PORT=8080

WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY polygen_scripts ./polygen_scripts

WORKDIR /app/polygen_scripts
EXPOSE 8080

CMD exec gunicorn --bind "0.0.0.0:${PORT}" --workers 1 --threads 8 --timeout 120 polygen:app
