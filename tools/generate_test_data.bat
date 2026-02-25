@echo off
REM generate_test_data.bat — Generate synthetic C. elegans-like test data (Windows)
REM
REM Creates a small but realistic genome, GFF3, and CDS FASTA for testing
REM the frequency_analyzer.py pipeline when real genome data is unavailable.
REM
REM Requires: Python 3
REM Usage: generate_test_data.bat [--outdir DIR]

setlocal enabledelayedexpansion

REM Resolve tools directory relative to this script
set "SCRIPT_DIR=%~dp0"

REM Resolve default output directory to a clean absolute path (no "..")
pushd "%SCRIPT_DIR%.."
set "PROJECT_DIR=%CD%"
popd
set "OUTDIR=%PROJECT_DIR%\organisms\c-elegans\data"

REM Parse arguments
:parse_args
if "%~1"=="" goto :run
if "%~1"=="--outdir" (
    set "OUTDIR=%~2"
    shift
    shift
    goto :parse_args
)
if "%~1"=="-h" goto :usage
if "%~1"=="--help" goto :usage
echo Unknown argument: %~1
goto :usage

:usage
echo Usage: generate_test_data.bat [--outdir DIR]
echo.
echo Generates synthetic C. elegans-like test data for the frequency_analyzer pipeline.
echo.
echo Options:
echo   --outdir DIR   Output directory (default: organisms\c-elegans\data)
echo   -h, --help     Show this help message
exit /b 0

:run
echo === GeneLang Test Data Generator ===
echo Output directory: %OUTDIR%
echo.

REM --- Find a working Python 3 interpreter ---
REM Try each candidate and verify it actually runs.
set "PYTHON="

REM 1) Try 'py -3' (Windows Python Launcher — most reliable)
py -3 --version >nul 2>&1
if !errorlevel! equ 0 (
    set "PYTHON=py -3"
    goto :found_python
)

REM 2) Try 'python3'
python3 --version >nul 2>&1
if !errorlevel! equ 0 (
    set "PYTHON=python3"
    goto :found_python
)

REM 3) Try 'python'
python --version >nul 2>&1
if !errorlevel! equ 0 (
    set "PYTHON=python"
    goto :found_python
)

REM Nothing works
echo ERROR: Python 3 not found.
echo.
echo Tried: py -3, python3, python — none of them work.
echo.
echo Please install Python 3 from https://www.python.org/downloads/
echo and make sure to check "Add Python to PATH" during installation.
exit /b 1

:found_python
echo Using: !PYTHON!
echo.

REM Run the Python script
!PYTHON! "%SCRIPT_DIR%generate_test_data.py" --outdir "%OUTDIR%"
if !errorlevel! neq 0 (
    echo.
    echo ERROR: Test data generation failed.
    exit /b 1
)

echo.
echo === Test data generation complete ===
echo.
echo To run the fingerprinting pipeline against this data:
echo.
echo   !PYTHON! tools\frequency_analyzer.py all --genome "%OUTDIR%\genome.fna" --gff "%OUTDIR%\annotations.gff3" --cds "%OUTDIR%\cds.fna" --outdir organisms\c-elegans\fingerprint
exit /b 0
