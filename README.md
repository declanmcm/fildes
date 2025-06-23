# Fildes: Interactive Digital Filter Designer Plugin

**Author:** Declan McMahon  
**License:** MIT  
**Version:** 1.0.0  

## Overview

**Fildes** is a versatile VST plugin for interactive digital filter creation and manipulation. It bridges the gap between audio developers and musicians by providing an intuitive and educational interface for designing digital filters. Unlike typical EQs or filter plugins, Fildes exposes the mathematical underpinnings of filters — including pole-zero plots, transfer functions, and frequency/phase responses — and allows real-time manipulation while processing audio.

## Features

- **Real-time audio processing** with dynamic filter coefficient updates.
- **Interactive pole-zero plot**: Drag poles and zeros to shape your filter in real-time.
- **Transfer function editor**: Enter numerator/denominator coefficients directly.
- **RBJ cookbook filter generator**: Quickly create standard biquad filters (low pass, high pass, band pass, notch, shelving, peak, all pass).
- **Visualizations**:
  - Frequency response (magnitude + phase)
  - Nyquist plot
  - Pole-zero diagram
- **Amplitude limiting tool** to keep filter gain within desired bounds.
- **Educational help system** explaining each UI component and DSP concept.

## Build Requirements

- **C++11 or later**
- **DISTRHO Plugin Framework (DPF)**
- **Cairo (for UI rendering)**
- **Eigen (for root-finding and matrix operations)**
- **CMake** (recommended for build configuration)
- **VST3 / CLAP compatible DAW for testing** (e.g., Carla, Audacity)

## Build Instructions

Install [DPF](https://github.com/DISTRHO/DPF), then navigate into the plugins directory.

```bash
git clone <this-repo>
cd <this-repo>
# Ensure DPF is set up correctly
make DPF_PATH=../DPF DPF_PLUGIN_HAS_UI=true
```

Alternatively, use your preferred CMake-based build setup if you’ve adapted it.

## Usage

- Load the plugin in your DAW.
- Enter transfer function coefficients or use the filter generator.
- Drag poles/zeros to shape your filter and hear the result instantly.
- Adjust gain or set amplitude limits to keep the filter within bounds.
- Use the help buttons (?) to learn about different parts of the interface.

## Files

- `Fildes.cpp`: Audio processing logic, handling filter state and transfer function updates.
- `FildesUI.cpp`: Custom user interface, interactive plots, and control widgets.
- `EditableText.hpp`: Component for text-based parameter editing.
- `Thesis_proposal.pdf`: Detailed proposal outlining project goals, background, and methodology.
- `Marksheet_thesis.pdf`: Assessment criteria and guidelines.

## Goals

Fildes aims to:
- Enable musicians to design and modify filters intuitively.
- Provide educational insights into digital filter design.
- Serve both as a creative tool and a learning platform.

## References

- [Steinberg VST SDK](https://steinbergmedia.github.io/vst3_dev_portal/pages/index.html)  
- [MathWorks Filter Designer](https://www.mathworks.com/help/signal/ug/introduction-to-filter-designer.html)  
- [Max/MSP](https://docs.cycling74.com/max8)  
- [Pure Data](https://puredata.info/)  

## License

MIT License — see source files for details.
