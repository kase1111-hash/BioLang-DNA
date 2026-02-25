@echo off
REM download_genome.bat — Download C. elegans WBcel235 reference genome from NCBI
REM Windows equivalent of download_genome.sh
REM
REM Requires: curl (built into Windows 10+) and tar (for gzip decompression)
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
echo [%TIME%] Assembly: %ASSEMBLY% (%ACCESSION%)
echo [%TIME%] Target:   %DATA_DIR%

if not exist "%DATA_DIR%" mkdir "%DATA_DIR%"
cd /d "%DATA_DIR%"

set "FAILED=0"

REM 1. Full genome FASTA
call :download_with_retry "%BASE_URL%/%ACCESSION%_%ASSEMBLY%_genomic.fna.gz" "genome.fna.gz"
if errorlevel 1 set /a FAILED+=1

REM 2. Gene annotations (GFF3)
call :download_with_retry "%BASE_URL%/%ACCESSION%_%ASSEMBLY%_genomic.gff.gz" "annotations.gff3.gz"
if errorlevel 1 set /a FAILED+=1

REM 3. Predicted proteome
call :download_with_retry "%BASE_URL%/%ACCESSION%_%ASSEMBLY%_protein.faa.gz" "proteome.faa.gz"
if errorlevel 1 set /a FAILED+=1

REM 4. Coding sequences (CDS)
call :download_with_retry "%BASE_URL%/%ACCESSION%_%ASSEMBLY%_cds_from_genomic.fna.gz" "cds.fna.gz"
if errorlevel 1 set /a FAILED+=1

if %FAILED% gtr 0 (
    echo.
    echo [%TIME%] WARNING: %FAILED% download(s) failed.
    echo [%TIME%] Download manually from: %BASE_URL%/
    exit /b 1
)

REM --- Decompress ---
if "%DECOMPRESS%"=="1" (
    echo.
    echo [%TIME%] --- Decompressing ---
    for %%f in (genome.fna.gz annotations.gff3.gz proteome.faa.gz cds.fna.gz) do (
        if exist "%%f" (
            tar -xzf "%%f" 2>nul
            if !errorlevel! neq 0 (
                REM Fallback: try PowerShell for decompression
                powershell -Command "& { $input = '%DATA_DIR%\%%f'; $output = '%DATA_DIR%\%%~nf'; $stream = [System.IO.File]::OpenRead($input); $gzip = New-Object System.IO.Compression.GzipStream($stream, [System.IO.Compression.CompressionMode]::Decompress); $outStream = [System.IO.File]::Create($output); $gzip.CopyTo($outStream); $outStream.Close(); $gzip.Close(); $stream.Close() }" 2>nul
            )
            if exist "%%~nf" (
                del "%%f"
                echo [%TIME%]   Decompressed %%f
            )
        )
    )
)

REM --- Verification ---
echo.
echo [%TIME%] --- Verification ---
for %%f in (genome.fna annotations.gff3 proteome.faa cds.fna) do (
    if exist "%%f" (
        for %%s in ("%%f") do echo [%TIME%]   OK %%f — %%~zs bytes
    ) else if exist "%%f.gz" (
        for %%s in ("%%f.gz") do echo [%TIME%]   OK %%f.gz — %%~zs bytes ^(compressed^)
    ) else (
        echo [%TIME%]   MISSING %%f
    )
)

echo.
echo [%TIME%] === Acquisition complete ===
exit /b 0

REM ============================================================
REM Subroutine: download_with_retry URL OUTPUT
REM ============================================================
:download_with_retry
set "URL=%~1"
set "OUTPUT=%~2"
set "MAX_RETRIES=4"
set "DELAY=2"

for /l %%a in (1,1,%MAX_RETRIES%) do (
    echo [%TIME%] Downloading %~2 (attempt %%a/%MAX_RETRIES%)...
    curl -sfL -o "%OUTPUT%" "%URL%"
    if !errorlevel! equ 0 (
        for %%s in ("%OUTPUT%") do echo [%TIME%]   OK Downloaded %~2 ^(%%~zs bytes^)
        exit /b 0
    )
    if %%a lss %MAX_RETRIES% (
        echo [%TIME%]   Retrying in !DELAY!s...
        timeout /t !DELAY! /nobreak >nul
        set /a DELAY*=2
    )
)

echo [%TIME%]   FAILED to download %~2 after %MAX_RETRIES% attempts
exit /b 1
