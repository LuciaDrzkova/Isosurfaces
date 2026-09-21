@echo off
rem Rebuilds what changed, then starts the program.
rem Windows: double-click this file.  PowerShell: .\run.bat
rem Anything after the name is passed on, e.g. .\run.bat assets\sphere.obj -f sphere

cd /d "%~dp0"

where cmake >nul 2>nul
if errorlevel 1 (
    echo CMake is not installed or not on the PATH.
    echo Install it from https://cmake.org/download/ and choose "Add CMake to the system PATH".
    echo Then close and reopen this window and run this file again.
    pause
    exit /b 1
)

cmake -S . -B build || goto :failed
cmake --build build --config Release --parallel || goto :failed

build\Release\Isosurfaces.exe %*
exit /b

:failed
echo.
echo Build failed - see the messages above.
pause
exit /b 1
