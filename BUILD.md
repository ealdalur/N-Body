# Building N-Body Simulation

This document describes how to install all required dependencies and build the N-Body simulation project on Windows and Linux.

## Dependencies

- **CMake** 3.15 or later
- **C++17** compatible compiler (MSVC 2019+, GCC 9+, Clang 10+)
- **SDL3** (window creation and input)
- **OpenGL** (rendering)
- **FFmpeg** with libx264 (video recording)
  - libavcodec
  - libavformat
  - libavutil
  - libswscale

---

## Windows (vcpkg + MSVC)

### 1. Prerequisites

- **Visual Studio 2019 or later** with the "Desktop development with C++" workload
- **Git** (for cloning vcpkg)
- **CMake** (included with Visual Studio, or install separately)

### 2. Install and configure vcpkg

Open a terminal (Developer Command Prompt, PowerShell, or Git Bash):

```bat
cd %USERPROFILE%\source\repos
git clone https://github.com/microsoft/vcpkg.git
cd vcpkg
.\bootstrap-vcpkg.bat
```

Optionally, add vcpkg to your PATH or set the `VCPKG_ROOT` environment variable:

```bat
set VCPKG_ROOT=%USERPROFILE%\source\repos\vcpkg
```

### 3. Install packages

```bat
vcpkg install sdl3:x64-windows
vcpkg install ffmpeg[x264]:x64-windows
```

OpenGL is provided by the Windows SDK and does not need a vcpkg package.

### 4. Configure with CMake

The project includes a CMake preset that points to the vcpkg toolchain file. From the project root:

```bat
cmake --preset MSVC
```

This uses the `MSVC` preset defined in `CMakePresets.json`, which sets the toolchain file to:
```
%USERPROFILE%\source\repos\vcpkg\scripts\buildsystems\vcpkg.cmake
```

If your vcpkg is installed in a different location, either:
- Edit `CMakePresets.json` and update the `CMAKE_TOOLCHAIN_FILE` path, or
- Pass it manually:

```bat
cmake -B build -DCMAKE_TOOLCHAIN_FILE="C:\path\to\vcpkg\scripts\buildsystems\vcpkg.cmake"
```

### 5. Build

```bat
cmake --build build --config Release
```

Or use the preset:

```bat
cmake --build --preset MSVC_Release
```

The compiled executable will be at `build\Release\N-Body.exe`.

### 6. Run

FFmpeg DLLs are automatically copied to the output directory by vcpkg's app-local deployment. Run from the project root so the program can find/write data files:

```bat
cd %PROJECT_ROOT%
build\Release\N-Body.exe
```

---

## Linux (Ubuntu/Debian)

### 1. Install build tools

```bash
sudo apt update
sudo apt install -y build-essential cmake git
```

### 2. Install dependencies

```bash
sudo apt install -y \
    libsdl3-dev \
    libgl-dev \
    libavcodec-dev \
    libavformat-dev \
    libavutil-dev \
    libswscale-dev \
    libx264-dev
```

> **Note:** SDL3 may not be available in the default Ubuntu repositories for older releases. If `libsdl3-dev` is not found, you can build SDL3 from source:
>
> ```bash
> git clone https://github.com/libsdl-org/SDL.git -b main
> cd SDL
> cmake -B build -DCMAKE_BUILD_TYPE=Release
> cmake --build build -j$(nproc)
> sudo cmake --install build
> cd ..
> ```

### 3. Configure with CMake

From the project root:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
```

Or using a preset:

```bash
cmake --preset gcc_Release
```

### 4. Build

```bash
cmake --build build -j$(nproc)
```

Or using a preset:

```bash
cmake --build --preset gcc_Release
```

The compiled executable will be at `build/N-Body`.

### 5. Run

```bash
cd /path/to/project
./build/N-Body
```

---

## CMake Presets Reference

| Preset Name   | Platform | Description                |
|---------------|----------|----------------------------|
| MSVC          | Windows  | Configure for MSVC x64     |
| MSVC_Release  | Windows  | Build Release with MSVC    |
| MSVC_Debug    | Windows  | Build Debug with MSVC      |
| gcc_Release   | Linux    | Configure + build Release  |
| gcc_Debug     | Linux    | Configure + build Debug    |

---

## Video Recording

To enable video recording, set the `RECORD_VIDEO` constant to `true` in `Simulation.h`:

```cpp
const bool RECORD_VIDEO = true;
```

When enabled, the simulation writes `output.mp4` (H.264, 30 fps) to the working directory. The file is finalized when the program exits normally (close the window or press Escape).

---

## Troubleshooting

### Windows: "FFMPEG not found"
Ensure you installed FFmpeg with the x264 feature flag:
```bat
vcpkg install ffmpeg[x264]:x64-windows
```
And verify the toolchain file is correctly passed to CMake.

### Windows: Missing DLLs at runtime
Run from Visual Studio or ensure you are running the executable from the build output directory where vcpkg copies the DLLs.

### Linux: "libsdl3-dev not found"
SDL3 is relatively new. Either build from source (see instructions above) or add a PPA that provides SDL3 packages.

### Linux: "avcodec.h not found"
Install the FFmpeg development libraries:
```bash
sudo apt install libavcodec-dev libavformat-dev libavutil-dev libswscale-dev
```

### Encoder error: "H.264 encoder not found"
The FFmpeg installation must include libx264 support. On Linux, ensure `libx264-dev` is installed. On Windows, use the `ffmpeg[x264]` vcpkg feature.
