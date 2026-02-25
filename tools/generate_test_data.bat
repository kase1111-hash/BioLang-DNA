@echo off
REM generate_test_data.bat — Generate synthetic C. elegans-like test data (Windows)
REM
REM Creates a small but realistic genome, GFF3, and CDS FASTA for testing
REM the frequency_analyzer.py pipeline when real genome data is unavailable.
REM
REM Requires: Python 3 in PATH
REM Usage: generate_test_data.bat [--outdir DIR]

setlocal enabledelayedexpansion

REM Resolve tools directory relative to this script
set "SCRIPT_DIR=%~dp0"

REM Default output directory
set "OUTDIR=%SCRIPT_DIR%..\organisms\c-elegans\data"

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

REM Check Python is available
where python >nul 2>&1
if !errorlevel! neq 0 (
    echo ERROR: Python not found in PATH.
    echo Please install Python 3 and ensure it is on your PATH.
    exit /b 1
)

REM Run the Python script
python "%SCRIPT_DIR%generate_test_data.py" --outdir "%OUTDIR%"
if !errorlevel! neq 0 (
    echo.
    echo ERROR: Test data generation failed.
    exit /b 1
)

echo.
echo === Test data generation complete ===
echo.
echo To run the fingerprinting pipeline against this data:
echo   python tools\frequency_analyzer.py all ^
echo     --genome "%OUTDIR%\genome.fna" ^
echo     --gff "%OUTDIR%\annotations.gff3" ^
echo     --cds "%OUTDIR%\cds.fna" ^
echo     --outdir organisms\c-elegans\fingerprint
exit /b 0
