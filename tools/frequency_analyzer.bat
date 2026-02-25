@echo off
REM frequency_analyzer.bat — Run the GeneLang frequency-analysis pipeline (Windows)
REM
REM Convenience wrapper that locates a working Python 3 interpreter and
REM calls frequency_analyzer.py.  When invoked with no arguments it runs
REM the full "all" analysis against the default C. elegans test data.
REM
REM Requires: Python 3 + packages listed in requirements.txt
REM
REM Usage:
REM   frequency_analyzer.bat                         (run all, default paths)
REM   frequency_analyzer.bat all   [OPTIONS]
REM   frequency_analyzer.bat codon [OPTIONS]
REM   frequency_analyzer.bat kmer  [OPTIONS]
REM   frequency_analyzer.bat density [OPTIONS]
REM   frequency_analyzer.bat repeats [OPTIONS]
REM   frequency_analyzer.bat -h

setlocal enabledelayedexpansion

REM Resolve directories to clean absolute paths (no "..")
set "SCRIPT_DIR=%~dp0"
pushd "%SCRIPT_DIR%.."
set "PROJECT_DIR=%CD%"
popd

set "DATA_DIR=%PROJECT_DIR%\organisms\c-elegans\data"
set "OUT_DIR=%PROJECT_DIR%\organisms\c-elegans\fingerprint"

REM --- Find a working Python 3 interpreter ---
set "PYTHON="

REM 1) Try 'py -3' (Windows Python Launcher)
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

REM --- Check for required Python packages ---
set "DEPS_OK=1"
for %%P in (numpy scipy pandas matplotlib seaborn) do (
    !PYTHON! -c "import %%P" >nul 2>&1
    if !errorlevel! neq 0 (
        set "DEPS_OK=0"
    )
)

if "!DEPS_OK!"=="0" (
    echo Missing required Python packages. Installing from requirements.txt ...
    echo.
    !PYTHON! -m pip install -r "%PROJECT_DIR%\requirements.txt"
    if !errorlevel! neq 0 (
        echo.
        echo ERROR: Failed to install dependencies.
        echo Please run manually:  !PYTHON! -m pip install -r requirements.txt
        exit /b 1
    )
    echo.
)

REM --- If the user passed arguments, forward them as-is ---
if not "%~1"=="" (
    echo === GeneLang Frequency Analyzer ===
    echo Using: !PYTHON!
    echo.
    !PYTHON! "%SCRIPT_DIR%frequency_analyzer.py" %*
    if !errorlevel! neq 0 (
        echo.
        echo ERROR: Frequency analysis failed.
        exit /b 1
    )
    exit /b 0
)

REM --- No arguments: run the full "all" analysis with default paths ---
echo === GeneLang Frequency Analyzer ===
echo Using: !PYTHON!
echo.
echo Command : all (codon + kmer + density + repeats)
echo Genome  : %DATA_DIR%\genome.fna
echo GFF3    : %DATA_DIR%\annotations.gff3
echo CDS     : %DATA_DIR%\cds.fna
echo Output  : %OUT_DIR%
echo.

REM Verify input files exist
set "MISSING="
if not exist "%DATA_DIR%\genome.fna"       set "MISSING=1"
if not exist "%DATA_DIR%\annotations.gff3" set "MISSING=1"
if not exist "%DATA_DIR%\cds.fna"          set "MISSING=1"

if defined MISSING (
    echo ERROR: One or more input files not found in %DATA_DIR%
    echo.
    echo Run generate_test_data.bat first to create the test data, or
    echo pass explicit paths:
    echo.
    echo   frequency_analyzer.bat all --genome PATH --gff PATH --cds PATH --outdir PATH
    exit /b 1
)

!PYTHON! "%SCRIPT_DIR%frequency_analyzer.py" all --genome "%DATA_DIR%\genome.fna" --gff "%DATA_DIR%\annotations.gff3" --cds "%DATA_DIR%\cds.fna" --outdir "%OUT_DIR%"
if !errorlevel! neq 0 (
    echo.
    echo ERROR: Frequency analysis failed.
    exit /b 1
)

echo.
echo === Frequency analysis complete ===
echo Results written to: %OUT_DIR%
exit /b 0
