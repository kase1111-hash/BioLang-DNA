@echo off
REM download_genome.bat — Download C. elegans WBcel235 reference genome from NCBI
REM Windows equivalent of download_genome.sh
REM
REM Requires: curl (built into Windows 10+)
REM Usage: download_genome.bat [--no-decompress]

setlocal enabledelayedexpansion

REM --- Configuration ---
set "ACCESSION=GCF_000002985.6"
set "ASSEMBLY=WBcel235"
set "BASE_URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/%ACCESSION%_%ASSEMBLY%"

REM Resolve data directory relative to this script
set "SCRIPT_DIR=%~dp0"
set "DATA_DIR=%SCRIPT_DIR%data"

set "DECOMPRESS=1"
if "%~1"=="--no-decompress" set "DECOMPRESS=0"

REM --- Main ---
echo [%TIME%] === C. elegans Genome Acquisition ===
echo [%TIME%] Assembly: %ASSEMBLY% ^(%ACCESSION%^)
echo [%TIME%] Target:   %DATA_DIR%

if not exist "%DATA_DIR%" mkdir "%DATA_DIR%"
cd /d "%DATA_DIR%"

set "FAILED=0"

REM 1. Full genome FASTA
call :download_with_retry "%BASE_URL%/%ACCESSION%_%ASSEMBLY%_genomic.fna.gz" "genome.fna.gz"
if !errorlevel! neq 0 set /a FAILED+=1

REM 2. Gene annotations (GFF3)
call :download_with_retry "%BASE_URL%/%ACCESSION%_%ASSEMBLY%_genomic.gff.gz" "annotations.gff3.gz"
if !errorlevel! neq 0 set /a FAILED+=1

REM 3. Predicted proteome
call :download_with_retry "%BASE_URL%/%ACCESSION%_%ASSEMBLY%_protein.faa.gz" "proteome.faa.gz"
if !errorlevel! neq 0 set /a FAILED+=1

REM 4. Coding sequences (CDS)
call :download_with_retry "%BASE_URL%/%ACCESSION%_%ASSEMBLY%_cds_from_genomic.fna.gz" "cds.fna.gz"
if !errorlevel! neq 0 set /a FAILED+=1

if !FAILED! gtr 0 (
    echo.
    echo [%TIME%] WARNING: !FAILED! download^(s^) failed.
    echo [%TIME%] Download manually from: %BASE_URL%/
    exit /b 1
)

REM --- Decompress ---
if "%DECOMPRESS%"=="0" goto :skip_decompress

echo.
echo [%TIME%] --- Decompressing ---
call :decompress_file "genome.fna.gz" "genome.fna"
call :decompress_file "annotations.gff3.gz" "annotations.gff3"
call :decompress_file "proteome.faa.gz" "proteome.faa"
call :decompress_file "cds.fna.gz" "cds.fna"

:skip_decompress

REM --- Verification ---
echo.
echo [%TIME%] --- Verification ---
call :verify_file "genome.fna"
call :verify_file "annotations.gff3"
call :verify_file "proteome.faa"
call :verify_file "cds.fna"

echo.
echo [%TIME%] === Acquisition complete ===
exit /b 0


REM ============================================================
REM Subroutine: download_with_retry URL OUTPUT
REM ============================================================
:download_with_retry
set "DL_URL=%~1"
set "DL_OUTPUT=%~2"
set "DL_DELAY=2"

set "DL_ATTEMPT=1"
:download_loop
if !DL_ATTEMPT! gtr 4 goto :download_failed

echo [%TIME%] Downloading %~2 ^(attempt !DL_ATTEMPT!/4^)...
curl -sfL -o "%DL_OUTPUT%" "%DL_URL%"
if !errorlevel! equ 0 (
    echo [%TIME%]   OK Downloaded %~2
    exit /b 0
)

if !DL_ATTEMPT! lss 4 (
    echo [%TIME%]   Retrying in !DL_DELAY!s...
    timeout /t !DL_DELAY! /nobreak >nul
    set /a DL_DELAY*=2
)
set /a DL_ATTEMPT+=1
goto :download_loop

:download_failed
echo [%TIME%]   FAILED to download %~2 after 4 attempts
exit /b 1


REM ============================================================
REM Subroutine: decompress_file GZFILE OUTFILE
REM ============================================================
:decompress_file
if not exist "%~1" exit /b 0

REM Try tar first (Windows 10 1803+)
tar -xzf "%~1" 2>nul
if exist "%~2" goto :decompress_ok

REM Fallback: PowerShell gzip decompression
powershell -NoProfile -Command ^
    "$fs = [IO.File]::OpenRead('%~1'); $gz = [IO.Compression.GzipStream]::new($fs, [IO.Compression.CompressionMode]::Decompress); $os = [IO.File]::Create('%~2'); $gz.CopyTo($os); $os.Close(); $gz.Close(); $fs.Close()" 2>nul

if not exist "%~2" (
    echo [%TIME%]   WARN Could not decompress %~1
    exit /b 1
)

:decompress_ok
del "%~1"
echo [%TIME%]   Decompressed %~1
exit /b 0


REM ============================================================
REM Subroutine: verify_file FILENAME
REM ============================================================
:verify_file
if exist "%~1" (
    for %%s in ("%~1") do echo [%TIME%]   OK %~1 -- %%~zs bytes
    exit /b 0
)
if exist "%~1.gz" (
    for %%s in ("%~1.gz") do echo [%TIME%]   OK %~1.gz -- %%~zs bytes ^(compressed^)
    exit /b 0
)
echo [%TIME%]   MISSING %~1
exit /b 1
