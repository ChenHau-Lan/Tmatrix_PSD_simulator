FROM python:3.11-slim

WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends gfortran \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY engine/py_Tmatrix_Mueller ./py_Tmatrix_Mueller
COPY tmatrix_web ./tmatrix_web

RUN cd py_Tmatrix_Mueller/fortran_tm \
    && gfortran -fdefault-real-8 py_tmat.f90 -o tmat_py.exe \
    && gfortran -fdefault-real-8 py_musingle.f90 -o mueller_py.exe \
    && chmod 755 tmat_py.exe mueller_py.exe

ENV HOST=0.0.0.0
ENV PORT=7860

EXPOSE 7860

CMD ["python", "-m", "tmatrix_web.backend.app.main"]
